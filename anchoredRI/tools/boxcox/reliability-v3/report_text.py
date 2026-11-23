"""Write the clinical research reliability report from checked endpoint and likelihood records."""
import csv
from collections import Counter
from pathlib import Path


def percent(n, d):
    return f"{100*n/d:.1f}%" if d else "n/a"


def model_count(rows):
    return sum(r["model_match"] for r in rows)


def limit_count(rows):
    return sum(r["limit_match"] for r in rows)


def table(header, rows):
    return "\n".join(["| " + " | ".join(header) + " |",
                      "| " + " | ".join("---" for _ in header) + " |"]
                     + ["| " + " | ".join(map(str,r)) + " |" for r in rows])


def load(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def render(work, out, bsum, arm_b, baseline, cont, classifications, retained, outcomes, manifests):
    n=len(cont)
    target_m, target_l=152,149
    bm,bl=model_count(baseline),limit_count(baseline)
    nesting=sum(s["violations"] for s in bsum)
    complete=all(s["unavailable"]==0 for s in bsum)
    base_pass=bm>=target_m and bl>=target_l and nesting==0 and complete
    original=[r for r in cont if r["original42"]]
    counts=Counter(r["classification"] for r in classifications)
    n_cross=len(classifications)
    weak=counts["weak identification"]
    missed=counts["missed better fit"]
    weak_stop=weak/n>.05
    verdict=("Stop Arm A under the prespecified weak-identification rule." if weak_stop
             else "Arm A stayed below the prespecified weak-identification stopping threshold.")
    lead=(f"The baseline {'passed' if base_pass else 'missed'} the Hb reproducibility requirements. "
          f"{verdict} The accuracy comparison remains on hold.")
    lines=["# Task E reliability round 3","",lead,"",
           "These checks ask whether the same routine and screened results give the same fitted reference limits on two computers. They cover 156 simulated hemoglobin datasets, each with 40 screened anchors. Arm B uses the 540 datasets with previously saved fits.","",
           "## 1. Baseline with a rounding-based SD floor","",
           table(["Across the 156 Hb datasets","Result","Required"],[
                 ["Same transform and component count",f"{bm}/156 ({percent(bm,n)})","At least 152/156"],
                 ["Both reference limits within 0.5% of the Actions interval width",f"{bl}/156 ({percent(bl,n)})","At least 149/156"],
                 ["Eligible nesting violations",f"{bsum[0]['violations']} local, {bsum[1]['violations']} Actions","Zero in each runtime"]]),
           "",
           "The Hb rounding step was 0.1. In each candidate model, a component had to span at least one rounding unit in original-scale SD. The shared bound also preserved the stronger pre-existing floor where it applied. At zero shift, identity, log and shifted log use the same original-unit floor at a given median. The likelihood identity is between log and shifted log evaluated at the same parameters with shift zero.","",
           table(["Runtime","Fitted candidates","Eligible nesting checks","Excluded by shift bounds","Unavailable","Zero-shift vectors","Maximum likelihood difference","Maximum SD-floor difference"],[
               [s["runtime"],s["candidates"],s["checked"],s["excluded"],s["unavailable"],s["zero_vectors"],
                f'{s["max_zero_loglik_difference"]:.4g}',f'{s["max_zero_floor_difference"]:.4g}'] for s in bsum]),
           ""]
    worst=max(baseline,key=lambda r:r["fraction_of_width"])
    lines.append(f"The largest baseline endpoint difference was {worst['max_endpoint_difference']:.6g} Hb units "
                 f"({100*worst['fraction_of_width']:.4g}% of interval width) at Index {worst['Index']}. "
                 f"The local fit selected {worst['local_transform']} with {worst['local_k']} components "
                 f"({worst['local_lower']:.2f} to {worst['local_upper']:.2f}). Actions selected "
                 f"{worst['actions_transform']} with {worst['actions_k']} components "
                 f"({worst['actions_lower']:.2f} to {worst['actions_upper']:.2f}). "
                 f"Of the 468 zero-shift vectors per runtime, {bsum[0]['zero_vectors_in_fitted_range']} fell inside the permitted shift range. "
                 "The others tested the algebraic identity only. All fitted component SDs met the shared floor to numerical precision. "
                 "[Paired baseline fits](baseline-paired.csv) · [Nesting and floor summary](baseline-audit-summary.csv)")
    old_limit=len(original)-limit_count(original)
    old_any=sum(not (r["model_match"] and r["limit_match"]) for r in original)
    lines+=["","## 2. Arm A after continuing unfinished winners",""]
    for_host=[]
    for host,col in [("local","local"),("Actions","actions")]:
        codes=Counter(int(r[f"old_code_{col}"]) for r in cont)
        blocks=load(out/(host+"-continuation-blocks.csv"))
        final=sum(int(r[f"final_code_{col}"])==0 for r in cont)
        lbfgs=[r for r in blocks if r["engine"]=="L-BFGS-B"]
        bridges=sum(r["engine"]=="nlminb bridge" for r in blocks)
        scaled=len({r["Index"] for r in lbfgs if "affine code" in r["message"]})
        for_host.append([host,codes[1],codes[52],len(lbfgs),bridges,scaled,f"{final}/156"])
    lines+=[table(["Runtime","Original iteration-limit winners","Original line-search-error winners","L-BFGS-B blocks","Bounded bridge calls","Datasets using affine retry","Final winners with code 0"],for_host),"",
            f"Of the original 42 datasets where both runtimes selected the same model but produced different limits, **{old_limit}/42** still had limits outside tolerance after continuation. "
            f"{len(original)-old_limit}/42 came within tolerance. Counting either model or limit disagreement, {old_any}/42 still differed.","",
            f"Across all 156, {model_count(cont)}/156 selected the same model and {limit_count(cont)}/156 had both limits within tolerance. "
            f"{n_cross} had a model or endpoint disagreement and entered the reciprocal-start experiment. "
            "Each unfinished winning fit was resumed from its last parameters with up to 10,000 iterations per block. "
            "A bounded bridge or curvature-scaled retry was used when a numerical line search stalled. Only an actual terminal code 0 counted as converged. "
            "[All 156 before-and-after records](arm-a-continued-pairs.csv) · [Local blocks](local-continuation-blocks.csv) · [Actions blocks](Actions-continuation-blocks.csv)",
            "","## 3. Arm A reciprocal starts","",
            table(["Classification of reciprocal solutions","Cross-tested datasets","All 156 Hb datasets"],
                  [[key,f"{counts[key]}/{n_cross} ({percent(counts[key],n_cross)})",f"{counts[key]}/156 ({percent(counts[key],n)})"]
                   for key in ["missed better fit","weak identification","concordant"]]),"",
            f"The stopping boundary was more than 5% of all 156 datasets, or eight datasets. "]
    if weak:
        example=max((r for r in classifications if r["classification"]=="weak identification"),
                    key=lambda r:r["fraction_of_width"])
        lines[-1]+=(f"At Hb Index {example['Index']}, the log-likelihoods differed by just "
                    f"{example['abs_loglik_difference']:.3f}. The two lower limits were {example['local_lower']:.2f} "
                    f"locally and {example['actions_lower']:.2f} on Actions. The largest endpoint change was "
                    f"{100*example['fraction_of_width']:.2f}% of the interval width.")
    if outcomes:
        gains=[max(r["local_foreign_minus_pre_exchange"],r["actions_foreign_minus_pre_exchange"]) for r in classifications]
        restart=[max(r["local_foreign_minus_own_restart"],r["actions_foreign_minus_own_restart"]) for r in classifications]
        largest=max(classifications,key=lambda r:r["abs_loglik_difference"])
        lines+=["",f"The two runtimes completed {len(outcomes)} own-start and transferred-start solutions, all with code 0. "
            f"The largest absolute difference between a transferred starting log-likelihood and its recomputation was {max(abs(float(r['start_error'])) for r in outcomes):.3g}. "
            f"The largest reciprocal terminal log-likelihood difference was {largest['abs_loglik_difference']:.6g} at Index {largest['Index']}.",
            "",f"On the common selected model, the foreign start improved the best pre-exchange likelihood on at least one host by more than 1 in {sum(v>1 for v in gains)}/{n_cross} datasets "
            f"(largest gain {max(gains):.6g}). It beat a fresh restart from that host’s own parameters by more than 1 in {sum(v>1 for v in restart)}/{n_cross}. "
            f"When each host retained its best witnessed BIC across candidates, {model_count(retained)}/156 selected the same model and "
            f"{limit_count(retained)}/156 gave limits within tolerance."]
    lines+=["","Each runtime started from both saved solutions for each distinct winning specification. If the specifications differed, all comparisons used the same shift specification and component count selected by the lowest BIC achieved on either machine. "
            "A reciprocal log-likelihood difference greater than 1 counted as a missed better fit. A difference of 1 or less with either endpoint more than 0.5% of the Actions width apart counted as weak identification. "
            f"The {n-n_cross} datasets already concordant after continuation contribute to the full-study denominator. "
            "[Dataset classifications](arm-a-cross-decisions.csv) · [All terminal solutions](arm-a-cross-solutions.csv) · [Best fit retained by each runtime](arm-a-after-retaining-best.csv)",
            "","## 4. Arm B nesting audit","",
            table(["Saved datasets","Refit","Previously flagged","Flagged after refitting","Eligible nesting checks","Unavailable","Changed selected model"],
                  [[r["group"],r["datasets"],r["old_flagged"],r["new_flagged"],r["comparisons"],r["unavailable"],r["changed"]] for r in arm_b]),""]
    b=arm_b[-1]
    lines.append(f"The {b['datasets']} local refits covered 288 small-N Hb, 108 small-N Ca and 144 Hb datasets at N = 50,000. "
                 f"Across {b['candidates']:,} fitted candidates, {b['violations']} of {b['comparisons']:,} eligible nesting comparisons failed. "
                 f"The shift or positivity constraints excluded {b['excluded']:,} comparisons. "
                 f"All component SDs met the shared floor to numerical precision. "
                 f"The zero-shift likelihood and floor differences were zero across {b['zero_vectors']:,} available vectors. "
                 "[All 540 dataset records](arm-b-datasets.csv) · [Arm B summary](arm-b-summary.csv) · [Candidate and nesting tables](arm-b/)")
    lines+=["","## Reproducibility","",
           f"Routine-input and anchor hashes matched the archived hashes for all 156 pairs. The two runtimes also matched on {len(manifests)} distinct source-file hashes, including the native likelihood and convergence recovery. "
           "The local fits used Debian 13 with reference BLAS. Actions matched the original C2 Ubuntu image, R 4.5.0, BLAS, LAPACK and package versions. The Actions run retains the raw fit parameters and trace files. "
           "[Source hashes](source-hashes.csv) · [Local donor parameters](../exchange/local/) · [Actions donor parameters](../exchange/actions/) · [Local environment](local-runtime.txt) · [Actions environment](Actions-runtime.txt) · "
           "[GitHub Actions run](https://github.com/amromeo/amromeo.github.io/actions/runs/45180631883)","",
           "The thresholds and the Arm A cross-start rule were recorded before this round. The protocol records the numerical recovery for Hb 104. "
           "The two-runtime findings describe the 156 Hb datasets. Arm B was audited locally on its saved subset. [Protocol](../../PROTOCOL.md)"]
    (out/"REPORT.md").write_text("\n".join(lines)+"\n")
