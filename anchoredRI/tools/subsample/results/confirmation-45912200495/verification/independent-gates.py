#!/usr/bin/env python3
"""Recompute confirmation gates from the 1,440 released worker tables. No fitting or study-code imports."""
import argparse
import csv
from collections import Counter
from dataclasses import dataclass
import json
import math
from pathlib import Path
from statistics import median

import numpy as np
from scipy.special import betaincinv, ndtri, ndtri_exp, log_ndtr
from scipy.stats import truncnorm


def load(path):
    with path.open(newline="") as src:
        return list(csv.DictReader(src))


def assert_(value, detail):
    if not value:
        raise ValueError(detail)


def parse(x):
    return float(x) if x and x != "NA" else math.nan


def quantile(v,q):
    vals=sorted(float(x) for x in v if not math.isnan(x))
    if not vals: return math.nan
    pos=(len(vals)-1)*q; lo=math.floor(pos); hi=math.ceil(pos)
    left,right=vals[lo],vals[hi]
    if lo==hi or left==right: return left
    if left==-math.inf: return math.nan if right==math.inf else left
    if right==math.inf:return right
    return left+(pos-lo)*(right-left)


def match(value, field, label):
    if isinstance(value,(bool,str)):
        ok=(str(value) == field)
    else:
        reported=parse(field)
        ok=value==reported or (math.isnan(value) and math.isnan(reported)) or math.isclose(value,reported,abs_tol=1e-10,rel_tol=1e-10)
    if not ok:
        raise ValueError(f"{label} independently {value!r} versus published {field!r}")


def healthy(x, design):
    if x<=0: return -math.inf
    lam,mu,sd=(float(design[k]) for k in ("lambda","mu","sigma"))
    if lam==0: return (math.log(x)-mu)/sd
    # Condition a standard Gaussian above the transform boundary.
    lo=(-1./lam-mu)/sd
    y=(math.expm1(lam*math.log(x))/lam-mu)/sd
    if y==lo:
        # As x->0 analytic small-width log integral avoids floating-point cancellation.
        log_width=lam*math.log(x)-math.log(lam*sd)
        return float(ndtri_exp(log_width-lo**2/2-.5*math.log(2*math.pi)-log_ndtr(-lo)))
    if y<0:
        return float(ndtri_exp(truncnorm.logcdf(y,lo,np.inf)))
    return -float(ndtri_exp(truncnorm.logsf(y,lo,np.inf)))


@dataclass
class Result:
    d:float
    maximum:float
    valid:bool


def score(ends, design):
    if ends is None:
        return None
    lo,hi=ends
    if not math.isfinite(lo) or not math.isfinite(hi) or lo>=hi:
        return Result(math.inf,math.inf,False)
    z=[abs(t-healthy(x,design)) for x,t in zip(ends,(-1.959963984540054,1.959963984540054))]
    return Result(sum(z)/2,max(z),True)


def aggregate(j, rows, design):
    assert_(len(rows)==101 and all(int(x["j"])==j for x in rows),f"Wrong case fit count {j}")
    full,=[x for x in rows if x["fit_type"]=="full"]
    pkg=score((parse(full["lower"]),parse(full["upper"])),design) if full["status"]=="reported" else None
    results={}
    for bank in ("A","B"):
        subset=[x for x in rows if x["bank"]==bank]
        assert_(sorted(int(x["draw"]) for x in subset)==list(range(1,51)), f"Incomplete {j} {bank}")
        pairs=[(parse(x["raw_lower"]),parse(x["raw_upper"])) for x in subset]
        pairs=[(a,b) for a,b in pairs if math.isfinite(a) and math.isfinite(b) and a<b]
        median_ends=tuple(median(side) for side in zip(*pairs)) if pairs else None
        # At most 25 possible declines is equivalent to having at least 25 definite passes.
        approve=sum(x["first_status"]=="reported" for x in subset)>=25
        if full["first_status"]=="reported":
            s,policy="a",score(median_ends,design)
        elif pkg is not None:
            s,policy="b",pkg if pkg.valid else None
        elif approve and median_ends:
            s,policy="c",score(median_ends,design)
        else:
            s,policy="d",None
        results[bank]=dict(j=j,step=s,pkg=pkg,policy=policy,frac=float(design["fraction"]))
    return results


def error(result):
    return math.nan if result is None else result.d


def metric(rows, ci=False):
    package=np.array([error(r["pkg"]) for r in rows]);candidate=np.array([error(r["policy"]) for r in rows])
    e1=np.array([bool(r["pkg"] and r["pkg"].maximum>5) for r in rows])
    e2=np.array([bool(r["policy"] and r["policy"].maximum>5) for r in rows])
    p,c=~np.isnan(package),~np.isnan(candidate);common=p&c
    all_p=np.where(p,np.minimum(package,5),5);all_c=np.where(c,np.minimum(candidate,5),5)
    def paired(ix):
        selected=ix[common[ix]]
        if len(selected):
            with np.errstate(invalid="ignore"):
                dx=candidate[selected].mean()-package[selected].mean()
                dxcap=(np.clip(candidate[selected],0,5)-np.clip(package[selected],0,5)).mean()
        else:
            dx=dxcap=math.nan
        return dx,dxcap,(all_c[ix]-all_p[ix]).mean()
    ids=np.arange(len(rows));d,dc,da=paired(ids)
    out=dict(assigned=len(rows),common_returns=int(common.sum()),policy_returns=int(c.sum()),package_returns=int(p.sum()),
             lost_package_returns=int((p & ~c).sum()),added_returns=int((c & ~p).sum()),
             policy_invalid_returns=sum(r["policy"] is not None and not r["policy"].valid for r in rows),
             common_uncapped_difference=d,common_capped_difference=dc,all_assignment_difference=da,
             policy_extremes=int(e2.sum()),package_extremes=int(e1.sum()),
             policy_all_assignment=float(all_c.mean()),package_all_assignment=float(all_p.mean()))
    if ci:
        # Explicitly copy 8 original *rows* for each sampled cell. No aggregation of cell numerator/denominator.
        gids=sorted({(r["j"]-1)//8 for r in rows})
        grouped=np.array([[i for i,r in enumerate(rows) if (r["j"]-1)//8 == gid] for gid in gids],int)
        assert_(grouped.shape==(len(gids),8),"Expected eight independent replicates per cell")
        rng=np.random.Generator(np.random.PCG64(37000000)); distributions=[[] for _ in range(3)]
        for _ in range(2000):
            drawn=rng.choice(len(gids),size=len(gids),replace=True)
            for bucket,x in zip(distributions,paired(grouped[drawn].ravel())): bucket.append(x)
        for key,dist in zip(("uncapped","capped","all_assignment"),distributions):
            out[key+"_ci_lower"]=quantile(dist,.025)
            out[key+"_ci_upper"]=quantile(dist,.975)
            out[key+"_bootstrap_undefined"]=sum(math.isnan(x) for x in dist)
    return out


def json_safe(data):
    if isinstance(data,dict): return {str(k):json_safe(v) for k,v in data.items()}
    if isinstance(data,(float,np.floating)) and not math.isfinite(data): return str(data)
    if isinstance(data,np.generic):return data.item()
    return data


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("records","design","report","output"): parser.add_argument("--"+name,required=True,type=Path)
    a=parser.parse_args()
    ds_rows=[r for r in load(a.design) if int(r["N"])==50000]
    by={int(r["j"]):r for r in ds_rows}
    fits_files=[x for x in a.records.rglob("fits.csv") if x.parent.name.isdecimal()]
    files={int(x.parent.name):x for x in fits_files}
    assert_(len(ds_rows)==len(by)==len(fits_files)==len(files)==1440 and set(by)==set(files),"Requires exactly the complete 1,440 cases")
    archive=[aggregate(j,load(files[j]),by[j]) for j in sorted(by)]
    metrics=load(a.report/"metrics.csv")
    rescue=load(a.report/"recoveries.csv");gates=load(a.report/"gates.csv")
    checks=0; result={}
    for bank in "AB":
        raw=[r[bank] for r in archive]
        subs={"all":raw,"<=20%":[r for r in raw if r["frac"]<=.2],">20%":[r for r in raw if r["frac"]>.2]}
        computed={}
        for group,rows in subs.items():
            axis="all" if group=="all" else "pathology"
            target,=[r for r in metrics if (r["estimator"],r["bank"],r["group"],r["axis"])==("policy",bank,group,axis)]
            computed[group]=metric(rows,ci=True)
            for k,v in computed[group].items():match(v,target[k],f"{bank}/{group}/{k}");checks+=1
        rescued=[r["policy"] for r in raw if r["step"]=="c"]
        n=len(rescued); e=sum(r.maximum>5 for r in rescued); mean_d=np.mean([r.d for r in rescued]) if n else math.nan
        exactlo=0. if n and e==0 else (betaincinv(e,n-e+1,.025) if n else math.nan)
        exacthi=1. if n and e==n else (betaincinv(e+1,n-e,.975) if n else math.nan)
        resc={"rescues":n,"rescue_endpoint_extremes":e,"rescue_mean_d":mean_d,"exact95_low":exactlo,"exact95_high":exacthi}
        published,=[r for r in rescue if r["axis"]=="all" and r["bank"]==bank]
        for k,v in resc.items():match(v,published[k],f"{bank}/step c/{k}");checks+=1
        stepgate="insufficient evidence" if n<10 else ("pass" if mean_d<.5 and e==0 else "fail")
        pop=computed["all"]; rawbetter=pop["common_uncapped_difference"]<0
        def status(y): return "pass" if y else "fail"
        accuracy=status(rawbetter and all(m["common_returns"]>0 and m["common_capped_difference"]<=.02 for key,m in computed.items() if key!="all"))
        coverage=status(pop["policy_invalid_returns"]==0 and pop["lost_package_returns"]==0 and
                        all(m["policy_extremes"]<=m["package_extremes"] for m in computed.values()))
        allscore=status(pop["all_assignment_difference"]<0)
        valid_u=not any(math.isnan(pop[k]) for k in ("uncapped_ci_lower","uncapped_ci_upper","common_uncapped_difference"))
        valid_c=not any(math.isnan(pop[k]) for k in ("capped_ci_lower","capped_ci_upper","common_capped_difference"))
        labels=dict(lower_observed_mean_error=bool(rawbetter), improved_accuracy=bool(valid_u and pop["uncapped_ci_upper"]<0),
                    no_worse=bool(valid_c and pop["capped_ci_upper"]<.02),accuracy_gate=accuracy,coverage_and_extremes_gate=coverage,
                    all_assignment_gate=allscore,step_c_evidence=stepgate,
                    complete_four_step_gate="fail" if "fail" in (accuracy,coverage,allscore,stepgate) else stepgate)
        gate,=[r for r in gates if r["bank"]==bank]
        for k,v in labels.items():match(v,gate[k],f"{bank}/gate/{k}");checks+=1
        result[bank]=dict(routes=dict(Counter(r["step"] for r in raw)),metrics=computed,rescue=resc,gates=labels)
    output=dict(matched=True,cases=len(by),draws_and_full_fits=1440*101,checks=checks,bootstrap_draws_per_population=2000,
                bootstrap_seed=37000000,manual_row_level_bootstrap=True,results=result)
    a.output.write_text(json.dumps(json_safe(output),indent=2,allow_nan=False)+"\n")
    print(f"Matched {checks} independently reconstructed metrics, confidence bounds and gates on 1,440 datasets.")

if __name__=="__main__":
    main()
