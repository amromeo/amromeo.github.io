"""Build the review tables and time estimate from completed small-N and grid results."""
import csv
import json
import math
import statistics as st
import sys
from collections import Counter, defaultdict
from pathlib import Path

root, small_archive, large_archive = map(Path, sys.argv[1:4])
out = root / "tables"
assert (root / "SCORED.txt").exists()


def read(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def save(name, rows):
    with (out/name).open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def val(x):
    return float(x) if x not in ("", "NA", "NaN") else math.nan


def average(values):
    v = [x for x in values if math.isfinite(x)]
    return st.mean(v) if v else math.nan


def sf(x):
    return f"{val(x):.3f}"


def table(headers, rows):
    return "\n".join(["| "+" | ".join(headers)+" |","| "+" | ".join("---" for h in headers)+" |"]
                     + ["| "+" | ".join(map(str, row))+" |" for row in rows])


labels = dict(anchoredRI="Improved baseline", armB="Arm B", frozenLocal="Original frozen local fit",
              archivedC2="Archived C2", archivedGrid="Archived grid fit", anchorsOnly="Anchors only", refineR="refineR")
times = read(out/"timings.csv")
assert len(times)==len({r["dataset_id"] for r in times})==3780
b0 = min(val(x["began"]) for x in times)
b1 = max(val(x["ended"]) for x in times)
runtime = []
for stage in ["small","grid","combined"]:
    a = [r for r in times if stage=="combined" or r["stage"]==stage]
    fresh = [r for r in a if r["cached"]=="FALSE"]
    elapsed = max(val(r["ended"]) for r in a)-min(val(r["began"]) for r in a)
    runtime.append(dict(stage=stage,datasets=len(a),fresh=len(fresh),cached=len(a)-len(fresh),
                        elapsed_seconds=elapsed,fit_cpu_seconds=sum(val(r["cpu"]) for r in fresh),
                        sum_worker_fit_seconds=sum(val(r["fit_seconds"]) for r in fresh)))
save("runtime.csv",runtime)
# Estimate the next stage without touching its routine data or running a large fit.
def archive_seconds(where,n):
    result=defaultdict(list)
    for f in where.glob("custom-*.csv"):
        for r in read(f):
            if r["algo"]=="anchoredRI" and int(r["N"])==n:
                t=val(r["runtime"])
                if math.isfinite(t) and t>0: result[r["Analyte"]].append(t)
    return result

old5=archive_seconds(small_archive,5000)
old50=archive_seconds(large_archive,50000)
current5=defaultdict(list)
for r in times:
    if r["stage"]=="small" and int(r["N"])==5000 and r["cached"]=="FALSE":
        current5[r["analyte"]].append(val(r["cpu"]))
# All small-N Hb and all small-N Ca at some settings can be cached. Use the fresh Hb replay clock
# via its checkpoint timestamps only if no 5k measurement is available for that analyte.
# An independent fallback uses the ratios of older and newer routines among the measured analytes.
ratios=[average(current5[a])/average(old5[a]) for a in old5 if current5[a]]
multiplier=st.median(ratios)
previous=read(out/"anchoredRI-comparator.csv")
cached50=Counter(r["Analyte"] for r in previous if r["stage"]=="50000")
projection=[]
for a in sorted(old50):
    assert len(old50[a])==144
    direct=bool(current5[a])
    t5=average(current5[a]) if direct else average(old5[a])*multiplier
    factor=average(old50[a])/average(old5[a])
    factor2=st.median(old50[a])/st.median(old5[a])
    n=144-cached50[a]
    projection.append(dict(analyte=a,remaining=n,archived_large_fits=len(old50[a]),
        fresh_small_5000=len(current5[a]),small_5000_cpu_proxy=t5,archived_50k_to_5k_ratio=factor,
        source="measured at N=5000" if direct else "cross-analyte calibration",
        all_144_from_scratch_seconds=144*t5*factor,projected_cpu_seconds=n*t5*factor,median_ratio_sensitivity_seconds=n*t5*factor2))
save("large-runtime-projection.csv",projection)
assert sum(q["remaining"] for q in projection)==1296
eff=runtime[-1]["fit_cpu_seconds"]/runtime[-1]["elapsed_seconds"]
eta=sum(q["projected_cpu_seconds"] for q in projection)/eff
sensitivity=sum(q["median_ratio_sensitivity_seconds"] for q in projection)/eff
cold=sum(q["all_144_from_scratch_seconds"] for q in projection)/eff
# Audit runner-up gaps and endpoint movement, preserving all candidate BICs.
cands=defaultdict(list)
for r in read(out/"candidates.csv"): cands[r["dataset_id"]].append(r)
gaps=[]
for key,rows in cands.items():
    for label,cap in [("anchoredRI",3),("armB",5)]:
        eligible=[r for r in rows if int(r["k"])<=cap and math.isfinite(val(r["bic"]))]
        order=sorted(eligible,key=lambda r:(val(r["bic"]),{"identity":0,"log":1,"shifted_log":2}[r["transform"]],int(r["k"])))
        first,second=order[:2]
        ds=[abs(val(first[c])-val(second[c])) for c in ("lower","upper")]
        finite=[d for d in ds if math.isfinite(d)]
        width=val(first["upper"])-(val(first["lower"]) if math.isfinite(val(first["lower"])) else 0)
        gap_in_width=max(finite)/width if finite and width>0 else math.nan
        gaps.append(dict(dataset_id=key,stage="small" if key.startswith("r-") else "grid",method=label,
                         first=first["transform"],k=first["k"],second=second["transform"],second_k=second["k"],
                         bic_gap=val(second["bic"])-val(first["bic"]),lower_difference=ds[0],
                         upper_difference=ds[1],max_limit_difference=max(finite) if finite else math.nan,
                         max_difference_as_fraction_of_width=gap_in_width))
save("runner-up-gaps.csv",gaps)
near_ties=[]
for group in ("small","grid"):
    for method in ("anchoredRI","armB"):
        eligible=[r for r in gaps if r["stage"]==group and r["method"]==method]
        near=[r for r in eligible if r["bic_gap"]<=2]
        sizes=[r["max_difference_as_fraction_of_width"] for r in near if math.isfinite(r["max_difference_as_fraction_of_width"])]
        near_ties.append(dict(stage=group,method=method,attempts=len(eligible),within_two=len(near),
                       fraction_within_two=len(near)/len(eligible),finite_differences=len(sizes),
                       median_relative_difference=st.median(sizes) if sizes else math.nan,
                       largest_relative_difference=max(sizes) if sizes else math.nan))
save("runner-up-summary.csv",near_ties)
ri=read(out/"ribench-summary.csv");grid=read(out/"grid-summary.csv")
criteria=read(out/"small-success-criteria.csv");ep=read(out/"winning-optimizer-exits.csv")
assert len(previous)==3924 and len({r["dataset_id"] for r in previous})==3924
assert Counter(r["stage"] for r in previous)==dict(small=2880,grid=900,**{"50000":144})
pairs=read(out/"ribench-pairs.csv");gridpairs=read(out/"grid-pairs.csv")
def stage_table(rows,metric):
    selected=[r for r in rows if r["group"]=="all" and r["filter"]=="before"]
    after={r["method"]:r for r in rows if r["group"]=="all" and r["filter"]=="after"}
    return table(["Method","Before","After","Reported / attempted","Declined","Technical"],[
        [labels[r["method"]],sf(r[metric]),sf(after[r["method"]][metric]),
         f'{r["reported"]}/{r["attempts"]}',f'{r["declined"]} ({100*val(r["decline_rate"]):.1f}%)',r["technical"]] for r in selected])
report=["# Task E accuracy at small N and on the simulation grid","",
    "The improved three-component baseline and five-component Arm B were evaluated on 2,880 RIbench datasets and 900 grid datasets. "
    "Each fit used 40 screened anchors before any RIbench filtering. Arm A stopped during reliability testing.","",
    "The 50-dataset replay passed for 49/50 improved baseline fits and 49/50 Arm B fits. Two distinct datasets missed. The prespecified rule therefore required Actions fits for every dataset. Index 4298 changed from three to two shifted-log components for the baseline, and Index 3033 changed from four to five for Arm B. [Replay results](runtime-replay.csv)","",
    "All 120 fitting batches succeeded. The workflow took 59 minutes 50 seconds including setup, queues and the replay decision. The two Ubuntu 24.04 runner image revisions fitted 2,062 and 1,718 datasets, respectively (20261116.133.1 and 20261116.376.1). Both used R 4.5.0, RIbench 1.0.2 and the pinned fitting source. The frozen local baseline was fitted in the earlier local runtime. Numerical runtime effects can contribute to differences from this comparator. [Per-dataset provenance](runtime-provenance.csv) · [GitHub Actions run](https://github.com/amromeo/amromeo.github.io/actions/runs/45257620744)","",
    "## Small-N RIbench","",stage_table(ri,"capped_mean"),"",
    "Values are native RIbench mean absolute z deviations after its cutoff. The native cutoff omits dataset errors above 5. "
    "The uncapped mean, number included by the cutoff, extreme endpoints and counts of implausible fits are reported separately. "
    "After-safeguard means describe the intervals retained by each method. The archived anchors-only and refineR methods use their original reporting behavior.","",
    table(["Small-N criterion","Improved score","Frozen baseline threshold","Decision"],
          [[labels[r["method"]],sf(r["small_native_capped"]),sf(r["threshold"]),"Passed" if r["passed"]=="TRUE" else "Missed"] for r in criteria]),"",
    table(["Method","Filter","Scored","Uncapped mean absolute z","Mean absolute z > 5"],
          [[labels[r["method"]],r["filter"],r["scored"],sf(r["mean_abs_z"]),r["implausible"]] for r in ri if r["group"]=="all"]),"",
    "### Skewed sets with lambda = 0.432","",
    stage_table([dict(r,group="all") for r in ri if r["group"]=="lambda=0.432"],"capped_mean"),"",
    "## Simulation grid","",stage_table(grid,"mean_width_error"),"",
    "Grid errors are the larger absolute endpoint error as a percentage of the true reference interval width. Replicates 1 to 20 were used in development. The paired comparisons use the same routine values and anchors in every method.","",
    "## Change from the original fits","",
    table(["Dataset group","Comparator","Filter","Common pairs","Mean change","Wins / ties / losses"],
     [[group,labels[r["comparator"]],r["filter"],r["complete"],sf(r["mean_difference"]),f'{r["wins"]} / {r["ties"]} / {r["losses"]}']
      for group,items in [("Small N, uncapped mean absolute z",pairs),("Grid, width-error percentage points",gridpairs)] for r in items
      if r["arm"]=="anchoredRI" and r["comparator"] in ("archivedC2","archivedGrid","frozenLocal") and r.get("stratum",r.get("group"))=="all"]),"",
    "A negative paired change favors the improved baseline. Paired after-safeguard results include only datasets on which both methods reported an interval.",""]
report.extend(['Arm B tied the improved baseline on all 900 grid datasets. Cases 2, 19 and 20 crossed the two-MCSE threshold because the paired fits were nearly identical. Their mean increases over the frozen local fit were 0.000374, 0.000143 and 0.000176 percentage points of interval width, respectively. All three were accepted as numerically negligible. [Reviewed scenarios](grid-review.csv)', ""])
report+=["The distribution and pathological-fraction table reports paired differences, wins, ties, losses and the common reporting denominator before and after safeguarding. [RIbench paired comparisons](ribench-pairs.csv) · [Grid comparisons and MCSE](grid-pairs.csv)","",
    "## BIC near ties","",
    table(["Stage","Method","Runner-up within 2 BIC","Median largest endpoint difference in near ties (% of interval width)"],
       [[r["stage"],labels[r["method"]],f'{r["within_two"]}/{r["attempts"]} ({100*r["fraction_within_two"]:.1f}%)',
         f'{100*r["median_relative_difference"]:.2f}%'] for r in near_ties]),"",
    "One-sided CRP fits use the zero-to-upper reporting span. [Near-tie summary](runner-up-summary.csv)","",
    "## Winning optimizer exits","",
    table(["Stage","Method","Iteration-limit wins","Line-search-error wins"],
            [[stage,labels[m],sum(r["stage"]==stage and r["method"]==m and r["code"]=="1" for r in ep),
             sum(r["stage"]==stage and r["method"]==m and r["code"]=="52" for r in ep)]
             for stage in ("small","grid") for m in ("anchoredRI","armB")]),"",
    "The exit log includes the source start of the winning candidate. A logged exit leaves the prespecified estimator and safeguard decision intact.","",
    "## Runtime and N = 50,000 estimate","",
    table(["Stage","Completed","Fresh","Cached","Wall time (min)","Fresh fitting CPU (hours)"],
      [[r["stage"],r["datasets"],r["fresh"],r["cached"],f'{r["elapsed_seconds"]/60:.1f}',f'{r["fit_cpu_seconds"]/3600:.2f}'] for r in runtime]),"",
    f"The two stages ran concurrently on GitHub Actions with four workers per batch. Wall times span the first to last accepted fit and include queue gaps. Extraction and scoring took {int((root/'scoring-seconds.txt').read_text())/60:.1f} additional minutes. The projected fitting time for the remaining 1,296 large datasets is {eta/3600:.2f} hours at observed throughput "
    f"({sensitivity/3600:.2f} hours using median archived timing ratios). Repeating all 1,440 from scratch would take approximately {cold/3600:.2f} hours. This uses fresh small-N fitting times and the earlier C2 50,000-to-5,000 timing ratio for each analyte. "
    "A cross-analyte correction supplies the 5,000 estimate where all corresponding small fits were cached. This is a planning estimate for fitting. Scoring and publication add time.","",
    "The comparator CSV now includes 2,880 small-N datasets, 900 grid datasets, and the 144 pre-existing large-N Hb datasets. The N = 50,000 comparison was approved after this review.","",
    "[Comparator CSV](../../anchorRefine-comparator/anchoredRI-comparator.csv) · [Summary and paired tables](.)"]
(out/"REPORT.md").write_text("\n".join(report)+"\n")
(out/"completion.json").write_text(json.dumps(dict(fitted=3780,comparator_rows=3924,elapsed_seconds=b1-b0,
             projected_large_wall_seconds=eta,large_authorized=True, large_authorized_at_small_stage=False, large_release_date="2026-11-23"),indent=2)+"\n")
