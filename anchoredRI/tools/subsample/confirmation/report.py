#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-only
"""Reconcile and score the complete Task H synthetic confirmation after the run."""
import argparse
import csv
import gzip
import io
import json
import math
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.special import log_ndtr, ndtri, ndtri_exp
from report_records import BANKS, GENERATOR_SHA, collect, ensure, load_design, number, read, sha
from report_stats import N_BOOT, SEED, tables

HERE=Path(__file__).resolve().parent
ESTIMATORS=("package","policy","through_b")


def write_table(path, rows, fields=None):
    fields = list(fields or rows[0])
    with path.open("wb") as raw:
        stream = gzip.GzipFile(filename="",mode="wb",mtime=0,fileobj=raw) if path.suffix==".gz" else raw
        with io.TextIOWrapper(stream,encoding="utf-8",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=fields,extrasaction="ignore",lineterminator="\n")
            writer.writeheader()
            for r in rows:
                writer.writerow({k:("" if isinstance(r.get(k),float) and math.isnan(r[k]) else r.get(k,""))
                                 for k in fields})


def healthy_z(x, lam, mu, sig):
    """Independent normal/truncated-normal z for a positive original-scale value."""
    if x <= 0: return -math.inf
    if lam == 0: return (math.log(x)-mu)/sig
    a=(-1/lam-mu)/sig
    gap=math.exp(lam*math.log(x))/(lam*sig)
    t=a+gap
    log_normalizer=float(log_ndtr(-a))
    if t >= 0:
        return -float(ndtri_exp(log_ndtr(-t)-log_normalizer))
    if gap*(abs(a)+1)<1e-4:
        # Short-interval integral protects the near-lower-bound case.
        mid=a+gap/2
        lp=math.log(gap)-.5*mid**2-.5*math.log(2*math.pi)+math.log1p((mid**2-1)*gap**2/24)
    else:
        p=float(log_ndtr(t)); lp=p+math.log(-math.expm1(float(log_ndtr(a))-p))
    return float(ndtri_exp(lp-log_normalizer))


def scores(records, generator, rscript, work):
    input_records=[]
    for r in records:
        for estimator in ESTIMATORS:
            item={k:r[k] for k in ("j","bank","family","lambda","mu","sigma","truth_lower","truth_upper")}
            item.update(estimator=estimator,lower=r[estimator+"_lower"],upper=r[estimator+"_upper"])
            input_records.append(item)
    verify={"evaluations":0,"missing_endpoint_scores":0,"infinite_endpoint_scores":0,"max_absolute_finite_z_difference":0.}
    with tempfile.TemporaryDirectory(dir=work) as tmp:
        source,target=Path(tmp)/"endpoints.csv",Path(tmp)/"scores.csv"
        write_table(source,input_records)
        subprocess.run([str(rscript),"--vanilla",str(HERE/"report-score.R"),str(generator),str(source),str(target)],check=True)
        expected={(r["j"],r["bank"],est):r for r in records for est in ESTIMATORS}
        for observed in read(target):
            est=observed["estimator"]
            r=expected.pop((int(observed["j"]),observed["bank"],est))
            z=[]
            for side,p in (("lower",.025),("upper",.975)):
                x=r[est+"_"+side]; v=number(observed["z_"+side])
                if math.isfinite(x):
                    independently=float(ndtri(p))-healthy_z(x,r["lambda"],r["mu"],r["sigma"])
                    if math.isfinite(v) and math.isfinite(independently):
                        diff=abs(v-independently)
                        ensure(math.isclose(v,independently,rel_tol=1e-10,abs_tol=1e-10),
                               f"Score differs at {r['j']}/{r['bank']}/{est}/{side}")
                        verify["max_absolute_finite_z_difference"]=max(verify["max_absolute_finite_z_difference"],diff)
                    else:
                        ensure(v==independently,"Unexpected infinite endpoint score")
                        verify["infinite_endpoint_scores"]+=1
                elif math.isnan(x):
                    ensure(math.isnan(v),"Scored a missing endpoint")
                    verify["missing_endpoint_scores"]+=1
                r[est+"_z_"+side]=v
                z.append(v)
            d=number(observed["d"])
            if r[est+"_returned"] and r[est+"_invalid"]:
                d=math.inf
            elif r[est+"_returned"]:
                ensure(d==sum(map(abs,z))/2 or math.isclose(d,sum(map(abs,z))/2,rel_tol=2e-15),
                       "Dataset error differs from endpoints")
            else:
                ensure(math.isnan(d),"Nonreturn acquired a score")
            r[est+"_d"]=d
            r[est+"_extreme"]=r[est+"_returned"] and any(abs(v)>5 for v in z)
            r[est+"_all_assignment"]=min(d,5.) if r[est+"_returned"] else 5.
            verify["evaluations"]+=1
        ensure(not expected, "Stage 2 scorer omitted rows")
    return verify


def extra_outputs(rows):
    extremes=[]
    for r in rows:
        for est in ESTIMATORS:
            if est=="package" and r["bank"]=="B": continue
            for side in ("lower","upper"):
                z=r[est+"_z_"+side]
                if r[est+"_returned"] and abs(z)>5:
                    extremes.append(dict(j=r["j"],bank=r["bank"] if est!="package" else "",estimator=est,
                                         step=r["step"] if est!="package" else "",side=side,estimate=r[est+"_"+side],
                                         truth=r["truth_"+side],signed_z=z,d=r[est+"_d"],family=r["family"],
                                         fraction=r["fraction"],disease=r["disease"]))
    seeds=[]
    by={(r["j"],r["bank"]):r for r in rows}
    for j in sorted({r["j"] for r in rows}):
        a,b=(by[j,bank] for bank in BANKS)
        if a["package_returned"]: continue
        ensure(a["package_status"]==b["package_status"], "Package depends on bank")
        ra,rb=(r["step"]=="c" for r in (a,b))
        seeds.append(dict(j=j,family=a["family"],fraction=a["fraction"],disease=a["disease"],
                          seed_pattern="both" if ra and rb else ("A only" if ra else ("B only" if rb else "neither")),
                          step_A=a["step"],step_B=b["step"],majority_A=a["majority"],majority_B=b["majority"],
                          package_status=a["package_status"],package_reason=a["package_reason"],
                          fallback_reason=a["fallback_reason"]))
    return extremes,seeds


def f(x, precision=4):
    if not math.isfinite(x): return "undefined" if math.isnan(x) else ("+∞" if x>0 else "−∞")
    return f"{x:.{precision}f}"


def bounds(r,prefix):
    return f"[{f(r[prefix+'_ci_lower'])}, {f(r[prefix+'_ci_upper'])}]"


def md_table(headers, rows):
    return "\n".join("| "+" | ".join(map(str,r))+" |" for r in [headers,["---"]*len(headers),*rows])


def narrative(out,metrics,rescues,reasons,routes,gates,seeds):
    pm=[r for r in metrics if r["estimator"]=="policy" and r["axis"]=="all"]
    through=[r for r in metrics if r["estimator"]=="through_b" and r["axis"]=="all"]
    a,b=pm
    rr=[r for r in rescues if r["axis"]=="all"]
    g=gates[0]
    text=["# Task H independent confirmation at N = 50,000",
          f"The four-step policy returned {a['policy_returns']:,} of 1,440 intervals under seed set A. The package returned {a['package_returns']:,}. Both returned an interval on {a['common_returns']:,} datasets. The primary four-step judgment was **{g['complete_four_step_gate']}**. Seed set B reused the same datasets and anchors.",
          md_table(["Method","Returns","Mean error on returns","All-assignment score","Intervals with an extreme endpoint"],
              [["Package",a["package_returns"],f(a["package_mean_d_released"]),f(a["package_all_assignment"]),a["package_extremes"]]]
              +[[r["bank"]+(" through step b" if r["estimator"]=="through_b" else " through step c"),
                 r["policy_returns"],f(r["policy_mean_d_released"]),f(r["policy_all_assignment"]),r["policy_extremes"]]
                 for r in [pm[0],through[0],pm[1],through[1]]]),
          "Errors use the two true healthy endpoints for each dataset. The all-assignment score caps the mean absolute endpoint error at 5 and assigns 5 to every nonreturn. An extreme interval has at least one endpoint with absolute healthy z error above 5. [All methods, breakdowns and bootstrap denominators](metrics.csv)",
          "## Paired accuracy and decisions",
          md_table(["Seed set","Common returns","Policy mean d","Package mean d","Uncapped difference [95% interval]","Capped difference [95% interval]","Better / tied / worse"],
                   [[r["bank"],r["common_returns"],f(r["policy_common_mean_d"]),f(r["package_common_mean_d"]),
                     f(r["common_uncapped_difference"])+" "+bounds(r,"uncapped"),f(r["common_capped_difference"])+" "+bounds(r,"capped"),
                     f"{r['wins']} / {r['ties']} / {r['losses']}"] for r in pm]),
          md_table(["Seed set","Lower observed mean error","Improved accuracy","No worse","Accuracy gate","Coverage and extremes","All-assignment gate","Step c","Four-step judgment"],
                   [[r["bank"],*["yes" if r[k] else "no" for k in ("lower_observed_mean_error","improved_accuracy","no_worse")],
                     *[r[k] for k in ("accuracy_gate","coverage_and_extremes_gate","all_assignment_gate","step_c_evidence","complete_four_step_gate")]] for r in gates]),
          "The paired differences subtract the package error. Improved accuracy requires the entire two-sided cell-resampled 95% interval for the uncapped difference below zero. No worse requires the upper 95% bound for the capped difference below +0.02. The accuracy gate also checks the capped point estimates in both pathology strata. The bootstrap resampled factorial cells 2,000 times and kept the eight replicates in each sampled cell together. [Decision table](gates.csv) · [Protocol](https://github.com/amromeo/amromeo.github.io/blob/140bd43/anchoredRI/tools/subsample/confirmation/PROTOCOL.md)",
          md_table(["Set A stratum","Common n","Capped change","Policy extremes","Package extremes","Lost package returns"],
                   [[r["group"],r["common_returns"],f(r["common_capped_difference"]),r["policy_extremes"],r["package_extremes"],r["lost_package_returns"]]
                    for r in metrics if (r["bank"],r["estimator"],r["axis"])==("A","policy","pathology")]),
          f"Across all assignments the score difference was {f(a['all_assignment_difference'])} {bounds(a,'all_assignment')} for A, and {f(b['all_assignment_difference'])} {bounds(b,'all_assignment')} for B. All bounds and gate decisions use full precision. Signed infinite differences enter the bootstrap, and the tables count undefined replicates. [Common-return, finite-pair and all-assignment statistics](metrics.csv)",
          "## How the policy routed datasets",
          md_table(["Step","Output","A assignments","B assignments","A returns","B returns"],
                [[step,label,*[next(r[field] for r in routes if r["bank"]==bank and r["axis"]=="all" and r["step"]==step)
                              for field in ("n","released") for bank in BANKS]]
                    for step,label in zip("abcd",["Median after full-data pass","Package fallback","Median after package nonreturn","No interval"])]),
          md_table(["Set","Package nonreturns","Step c returns","Mean d","Median d","Maximum d","Extreme intervals","Extreme rate, exact 95% interval","Judgment"],
                   [[r["bank"],r["package_nonreturns"],r["rescues"],f(r["rescue_mean_d"]),f(r["rescue_median_d"]),f(r["rescue_max_d"]),
                     r["rescue_endpoint_extremes"],f"[{f(r['exact95_low'])}, {f(r['exact95_high'])}]",gg["step_c_evidence"]]
                     for r,gg in zip(rr,gates)]),
          "With fewer than 10 recoveries the step-c evidence is insufficient. At 10 or more, success requires mean uncapped error below 0.5 and zero extreme endpoints. The exact binomial interval describes uncertainty among the observed recoveries in this factorial design. Each package nonreturn’s flags, votes, contributing draws and reasons are included. [Each package nonreturn](package-nonreturns.csv.gz) · [Reasons and majority outcomes](nonreturn-reasons.csv) · [Extreme endpoints in all methods](extreme-endpoints.csv) · [Step-c summaries](recoveries.csv) · [Routes by factor](routes.csv)",
          "The two seed sets recovered the same "+str(sum(r["seed_pattern"]=="both" for r in seeds))+" package nonreturns. A alone recovered "+str(sum(r["seed_pattern"]=="A only" for r in seeds))+", B alone recovered "+str(sum(r["seed_pattern"]=="B only" for r in seeds))+" and neither recovered "+str(sum(r["seed_pattern"]=="neither" for r in seeds))+". [Recovery status by dataset](seed-comparison.csv)",
          "## Healthy family",
          md_table(["Set A family","Assigned","Common returns","Uncapped change [95% interval]","All-assignment change","Additional returns"],
                   [[r["group"],r["assigned"],r["common_returns"],f(r["common_uncapped_difference"])+" "+bounds(r,"uncapped"),
                     f(r["all_assignment_difference"]),r["added_returns"]] for r in metrics
                     if (r["bank"],r["estimator"],r["axis"])==("A","policy","family")]),
          "The power-0.5 generating family has no exact match among anchoredRI’s fitting transforms. These family comparisons are descriptive. [Exact fraction, disease pattern, rounding and family tables](metrics.csv) · [Per-case decisions, endpoints and scores](cases.csv.gz) · [Readback and independent score verification](verification.json)"]
    (out/"REPORT.md").write_text("\n\n".join(text)+"\n")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ("records","design","generator","reference-hashes","output"):
        parser.add_argument("--"+key,required=True,type=Path)
    parser.add_argument("--rscript",default="Rscript")
    args=parser.parse_args()
    out=args.output.resolve()
    out.mkdir(parents=True,exist_ok=True)
    ensure(sha(args.generator)==GENERATOR_SHA,"Unexpected Stage 2 generator source")
    ds=load_design(args.design,args.reference_hashes)
    rows=collect(args.records,ds)
    cells={r["cell_id"] for r in rows}
    ensure(len(cells)==180 and all(sorted(r["replicate"] for r in rows if r["cell_id"]==c and r["bank"]=="A")==list(range(1,9)) for c in cells),
           "Incomplete Stage 2 factorial cell structure")
    verify=scores(rows,args.generator.resolve(),args.rscript,out)
    metrics,rescues,reasons,routes,gates=tables(rows)
    extremes,seeds=extra_outputs(rows)
    write_table(out/"cases.csv.gz",rows)
    write_table(out/"package-nonreturns.csv.gz",[r for r in rows if not r["package_returned"]],rows[0].keys())
    for fn,data in [("metrics.csv",metrics),("gates.csv",gates),("routes.csv",routes),("recoveries.csv",rescues)]:
        write_table(out/fn,data)
    write_table(out/"nonreturn-reasons.csv",reasons,["bank","first_reason","package_status","package_reason","fallback_status","fallback_reason","step","majority","n"])
    write_table(out/"extreme-endpoints.csv",extremes,["j","bank","estimator","step","side","estimate","truth","signed_z","d","family","fraction","disease"])
    write_table(out/"seed-comparison.csv",seeds,["j","family","fraction","disease","seed_pattern","step_A","step_B","majority_A","majority_B","package_status","package_reason","fallback_reason"])
    verify.update(datasets=len(ds),cells=len(cells),banks=2,draws=len(rows)*50,full_fit_records=len(ds),
                  reconstructed_routes=len(rows),input_hash_pairs_matched=len(ds),
                  source_design_sha256=sha(args.design),source_generator_sha256=sha(args.generator),
                  reference_hashes_sha256=sha(args.reference_hashes),bootstrap_draws=N_BOOT,bootstrap_seed=SEED,
                  numpy_version=np.__version__)
    (out/"verification.json").write_text(json.dumps(verify,indent=2)+"\n")
    narrative(out,metrics,rescues,reasons,routes,gates,seeds)
    print(json.dumps({"datasets":len(ds),"primary":gates[0],"sensitivity":gates[1]},indent=2))


if __name__=="__main__":
    main()
