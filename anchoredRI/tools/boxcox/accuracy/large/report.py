"""Assemble the large-N report and merge verified baseline rows into the comparator."""
import csv
import json
import math
import os
import statistics
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

root = Path(sys.argv[1])
out = root / "tables"
here = Path(__file__).resolve().parent
old_file = here.parents[1] / "anchorRefine-comparator" / "anchoredRI-comparator.csv"
assert (root / "SCORED.txt").is_file()

def read(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))

def write(path, rows):
    with path.open("w", newline="") as stream:
        w = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

def number(x):
    return float(x) if x not in (None, "", "NA", "NaN") else math.nan

def value(x, digits=4):
    v = number(x)
    return f"{v:.{digits}f}" if math.isfinite(v) else "NA"

def proportion(a, b):
    return f"{100 * a / b:.1f}%" if b else "NA"

def table(headers, rows):
    return "\n".join("| " + " | ".join(map(str, row)) + " |" for row in (headers, ["---"] * len(headers), *rows))

def merge_comparator():
    prior, fresh = read(old_file), read(out / "baseline-fits.csv")
    by_id = {r["dataset_id"]: r for r in fresh}
    assert len(prior) == len({r["dataset_id"] for r in prior})
    assert len(fresh) == len(by_id) == 1440
    assert list(prior[0]) == list(fresh[0])
    assert all(r["stage"] == "50000" and r["N"] == "50000" and r["method"] == "anchoredRI" for r in fresh)
    counts = Counter(r["stage"] for r in prior)
    assert counts["small"] == 2880 and counts["grid"] == 900 and counts["50000"] in (144, 1440)
    for r in prior:
        if r["stage"] == "50000":
            new = by_id[r["dataset_id"]]
            assert all(r[k] == new[k] for k in ("Index", "routine_hash", "anchor_hash", "ribench_hashcode"))
    old_lines, new_lines = old_file.read_bytes().splitlines(keepends=True), (out / "baseline-fits.csv").read_bytes().splitlines(keepends=True)
    assert len(old_lines) == len(prior) + 1 and len(new_lines) == len(fresh) + 1
    preserved = [line for line, r in zip(old_lines[1:], prior) if r["stage"] != "50000"]
    dest = out / "anchoredRI-comparator.csv"
    dest.write_bytes(b"".join([old_lines[0], *preserved, *new_lines[1:]]))
    rows = read(dest)
    assert len(rows) == len({r["dataset_id"] for r in rows}) == 5220
    assert rows[-1440:] == fresh
    assert Counter(r["stage"] for r in rows) == {"small":2880, "grid":900, "50000":1440}
    return counts["50000"]

replaced = merge_comparator()
labels = {"anchoredRI":"Improved baseline", "armB":"Arm B", "frozenActions":"Original frozen fit, same Actions runner",
          "archivedC2":"Original archived C2", "anchorsOnly":"Anchors only", "refineR":"refineR"}
summary = read(out / "ribench-summary.csv")
summaries = {(r["method"],r["filter"],r["group"]):r for r in summary}
outcomes = read(out / "ribench-outcomes.csv")
criteria = read(out / "large-success-criteria.csv")
selections = read(out / "component-counts-and-errors.csv")
runs = read(out / "timings.csv")
assert len(outcomes) == 6 * 1440 and len(runs) == len({r["dataset_id"] for r in runs}) == 1440
assert len({r["run"] for r in runs}) == 1 and runs[0]["run"] == os.environ["GITHUB_RUN_ID"]
assert all({r["Index"] for r in outcomes if r["method"] == m} ==
           {r["dataset_id"][2:] for r in runs} for m in labels)
candidates = defaultdict(list)
for row in read(out / "candidates.csv"):
    candidates[row["dataset_id"]].append(row)
gaps = []
rank = {"identity":0, "log":1, "shifted_log":2}
for m, cap, selected in (("anchoredRI", 3, "selected_baseline"), ("armB", 5, "selected_armB")):
    for run in runs:
        key = run["dataset_id"]
        rows = [r for r in candidates[key] if int(r["k"]) <= cap and math.isfinite(number(r["bic"]))]
        rows.sort(key=lambda r:(number(r["bic"]), rank[r["transform"]], int(r["k"])))
        winners = [r for r in rows if r[selected] == "TRUE"]
        gap = dict(dataset_id=key, method=m, eligible=len(rows), first="", k="", second="", second_k="",
                   bic_gap="", lower_difference="", upper_difference="", max_difference_fraction_of_width="")
        if winners and len(rows) >= 2:
            first, second = rows[:2]
            assert len(winners) == 1 and winners[0] == first
            upper = abs(number(first["upper"]) - number(second["upper"]))
            lower = abs(number(first["lower"]) - number(second["lower"]))
            lrl = number(first["lower"])
            width = number(first["upper"]) - (lrl if math.isfinite(lrl) else 0)
            diffs = [x for x in [lower, upper] if math.isfinite(x)]
            gap.update(first=first["transform"], k=first["k"], second=second["transform"], second_k=second["k"],
                       bic_gap=number(second["bic"])-number(first["bic"]), lower_difference=lower,
                       upper_difference=upper, max_difference_fraction_of_width=max(diffs)/width if diffs and width>0 else "")
        gaps.append(gap)
write(out / "runner-up-gaps.csv", gaps)
near = []
for m in ("anchoredRI", "armB"):
    subset = [r for r in gaps if r["method"] == m and math.isfinite(number(r["bic_gap"]))]
    tied = [r for r in subset if number(r["bic_gap"]) <= 2]
    diffs = [number(r["max_difference_fraction_of_width"]) for r in tied]
    finite = [d for d in diffs if math.isfinite(d)]
    near.append(dict(method=m, attempted=1440, ranked=len(subset), within_two=len(tied),
                     median_difference_fraction_of_width=statistics.median(finite) if finite else "",
                     max_difference_fraction_of_width=max(finite) if finite else ""))
write(out / "runner-up-summary.csv", near)

def performance(group):
    rows = []
    for m, label in labels.items():
        pre, post = [summaries[m, view, group] for view in ("before", "after")]
        rows.append([label, value(pre["capped_mean"]), value(post["capped_mean"]), pre["capped_n"], post["capped_n"],
            value(pre["mean_abs_z"]), value(post["mean_abs_z"]), pre["implausible"], post["implausible"],
            f"{pre['reported']}/{pre['attempts']}", f"{pre['declined']} ({proportion(int(pre['declined']),int(pre['attempts']))})", pre["technical"]])
    return table(["Method","Native before","Native after","In cutoff before","In cutoff after","Uncapped before","Uncapped after",
                  "Implausible before","Implausible after","Reported / attempted","Declined","Technical"], rows)

components = {(r["method"],r["components"],r["filter"]):r for r in selections}
component_rows = []
for group in ("1 to 3", "4", "5", "4 or 5"):
    pre, post = (components["armB",group,f] for f in ("before","after"))
    count = int(pre["selected"])
    component_rows.append([group, f"{count}/1440 ({proportion(count,1440)})", value(pre["native_capped"]),
        value(pre["mean_abs_z"]), pre["implausible"], f"{pre['declined']}/{count} ({proportion(int(pre['declined']),count)})" if count else "NA",
        value(post["native_capped"]), pre["reported"]])

pair = read(out / "ribench-pairs.csv")
pair_rows = [[labels[r["arm"]], labels[r["comparator"]], r["filter"], r["complete"], value(r["mean_difference"]),
              f"{r['wins']} / {r['ties']} / {r['losses']}"] for r in pair if r["stratum"] == "all"]
exits = read(out / "winning-optimizer-exits.csv")
qc_rows = []
for r in near:
    m = r["method"]
    effect = number(r["median_difference_fraction_of_width"])
    qc_rows.append([labels[m], r["ranked"], f"{r['within_two']}/1440 ({proportion(r['within_two'],1440)})",
                    proportion(effect, 1) if math.isfinite(effect) else "NA",
                    sum(q["method"]==m and q["iteration_limit"]=="TRUE" for q in exits),
                    sum(q["method"]==m and q["line_search"]=="TRUE" for q in exits)])
images = Counter(r["image"] for r in runs)
gen = read(root / "generation.csv")[0]
started = min(number(r["began"]) for r in runs)
finished = max(number(r["ended"]) for r in runs)
elapsed = dict(preparation_seconds=number(gen["finished"])-number(gen["started"]),
               elapsed_fitting_window_seconds=finished-started,
               sum_fit_worker_seconds=sum(number(r["fit_seconds"]) for r in runs),
               sum_fit_cpu_seconds=sum(number(r["fit_cpu_seconds"]) for r in runs),
               sum_frozen_validation_worker_seconds=sum(number(r["frozen_seconds"]) for r in runs),
               scoring_and_export_seconds=time.time()-int((root/"scoring-start.txt").read_text()),
               elapsed_from_generation_to_report_seconds=time.time()-number(gen["started"]))
write(out / "runtime.csv", [elapsed])
run_id = runs[0]["run"]
run_url = "https://github.com/amromeo/amromeo.github.io/actions/runs/" + run_id
verdict = "; ".join(labels[r["method"]] + (" passed both criteria" if r["passed"]=="TRUE" else
           " missed " + " and ".join(label for label, field in (("the score criterion","score_pass"),("the implausible-fit criterion","count_pass")) if r[field]!="TRUE")) for r in criteria).replace("; ", ". ")
verdict = verdict.replace(". arm", ". Arm") + "."
decisions = table(["Arm","Native score before","refineR threshold","Score criterion","Implausible or failed / 1,440","Count criterion (≤14)","Overall"],
    [[labels[r["method"]], value(r["native_capped"], 8), value(r["refineR_threshold"], 8),
      "Pass" if r["score_pass"]=="TRUE" else "Miss", r["implausible_or_failed"],
      "Pass" if r["count_pass"]=="TRUE" else "Miss", "Pass" if r["passed"]=="TRUE" else "Miss"] for r in criteria])
report = [
    "# Task E at 50,000 routine results", "", verdict, "", decisions, "",
    "Each arm had 1,440 attempts. RIbench omits a dataset from its native cutoff mean when its mean absolute endpoint z deviation exceeds 5. The primary comparison uses raw fitted limits before either safeguard. Technical failures and implausible fits count once toward the limit of 14. The full-precision threshold and the prespecified additional score for missing fits are in [the criteria table](large-success-criteria.csv).", "",
    "## Accuracy and coverage", "", performance("all"), "",
    "After-safeguard scores describe the intervals each method released. RIbench uses two endpoints for nine analytes and the 95th percentile alone for CRP. The LRT declines at p < 0.01. The lower-limit check declines two-sided intervals with a lower limit below every routine observation. Heavily mixed populations can still yield inaccurate intervals that pass both checks.", "",
    "### Skewed FT4 datasets with lambda = 0.432", "", performance("lambda=0.432"), "",
    "## Component selection in Arm B", "",
    table(["Components","Share of all attempts","Native before","Uncapped before","Implausible before","Declines among selected","Native after","Reported"],component_rows), "",
    "The improved baseline allows up to three components. The pooled four-or-five row includes both individual component groups. Shares use all 1,440 attempts. Decline rates use the number selecting that component count. [Individual component counts and errors for both methods](component-counts-and-errors.csv)", "",
    "## Paired comparisons", "",
    table(["Arm","Comparator","Filter","Complete pairs","Mean change in uncapped absolute z","Wins / ties / losses"], pair_rows), "",
    "Negative differences favor the arm. After safeguarding, each pair includes only datasets for which both methods reported an interval. The paired tables also give each distribution by pathological fraction at or below 20% versus above 20%. [Paired results](ribench-pairs.csv) · [Subgroup scores and denominators](ribench-summary.csv) · [Safeguard flags by error](ribench-safeguard-by-error.csv)", "",
    "## Fit diagnostics and runtime", "",
    table(["Method","Finite runner-up","Within 2 BIC","Median endpoint movement in near ties (% of width)","Iteration-limit wins","Line-search-error wins"], qc_rows), "",
    "For one-sided CRP, the reference width for a near tie runs from zero to the fitted upper limit. The recorded exits describe the optimizer for the selected solution. [Candidate gaps](runner-up-gaps.csv) · [Winning optimizer exits](winning-optimizer-exits.csv)", "",
    f"GitHub Actions fitted all 1,440 datasets in 24 batches, with four workers per batch. Input preparation took {elapsed['preparation_seconds']/60:.1f} minutes. The fitting window was {elapsed['elapsed_fitting_window_seconds']/60:.1f} minutes and the two arms used {elapsed['sum_fit_cpu_seconds']/3600:.2f} CPU hours. Extraction, scoring and report generation took {elapsed['scoring_and_export_seconds']/60:.1f} minutes. Runner image revisions and counts were " + ", ".join(f"{k} ({v})" for k,v in sorted(images.items())) + f". [Actions run]({run_url}) · [Timing data](runtime.csv)", "",
    "RIbench supplied a single pooled routine sample and 40 independent screened anchors per dataset. Each routine vector and anchor vector matched its archived hash. The unchanged round 3 fitter considered identity, log when eligible, and shifted log. Baseline and Arm B used the same first three component candidates. The original frozen estimator and its one-group study adapter were also fitted on each Actions runner. The archived C2, anchors-only and refineR results came from [the independent-anchor C2 run](https://github.com/amromeo/amromeo.github.io/actions/runs/45101564181). The previously accepted two-setting simulation grid was evaluated separately.", "",
    f"The comparator contains 5,220 unique datasets, including all 1,440 large-N baseline fits. The Actions estimates replace {replaced} previously exported large-N rows after both input hashes match. The 2,880 small-N rows and 900 grid rows are preserved byte for byte. [Comparator CSV](anchoredRI-comparator.csv)", ""]
(out/"REPORT.md").write_text("\n".join(report))
(out/"completion.json").write_text(json.dumps(dict(run=run_id, source=os.environ["GITHUB_SHA"],
    fitter="ec3dd35197759424d32a6a341320173e5db08c89", archive_run="45101564181", datasets=len(runs),
    comparator_rows=5220, previous_large_rows=replaced, images=images, criteria=criteria), indent=2)+"\n")
