"""Audit fitted likelihoods. No optimizer or reference limits are used."""
import argparse
import csv
import gzip
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

ORDER = ("identity", "log", "shifted_log", "boxcox_unshifted", "boxcox_shifted")
FITTED_SHIFT = {"shifted_log", "boxcox_shifted"}
FREE_SHAPE = {"boxcox_unshifted", "boxcox_shifted"}


def number(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return math.nan


def records(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def edges(family, k):
    if k > 1:
        yield family, k - 1, "adjacent K", False
    if family == "shifted_log":
        yield "log", k, "shifted vs unshifted", True
    if family in FREE_SHAPE:
        yield "identity", k, "free shape vs identity", False
        yield "log", k, "free shape vs log", family == "boxcox_shifted"
    if family == "boxcox_shifted":
        yield "shifted_log", k, "free shape vs shifted log", False
        yield "boxcox_unshifted", k, "shifted vs unshifted", True


def parameters(node, groups):
    f, k = node
    return 2*k + groups*(k-1) + (f in FREE_SHAPE) + (f in FITTED_SHIFT)


def write(path, rows):
    if not rows:
        return
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", newline="") as out:
        w = csv.DictWriter(out, fieldnames=list(rows[0]),lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def audit_fit(meta, method, tables):
    n, groups = int(meta["n"]), int(meta["groups"])
    positive = meta["positive"].upper() == "TRUE"
    zero = meta["zero_shift"].upper() == "TRUE"
    own = tables[method]
    candidates = dict(tables.get("baseline", {})) if method == "armA" else {}
    candidates.update(own)
    original = {}; lower = {}; comparisons = []
    for node, item in candidates.items():
        bic = number(item["bic"])
        if math.isfinite(bic) and item["eligible"].upper() == "TRUE":
            original[node] = (parameters(node, groups)*math.log(n) - bic)/2
    for node in sorted(candidates, key=lambda t:(t[1], ORDER.index(t[0]))):
        f, k = node
        item = candidates[node]
        if item["eligible"].upper() != "TRUE":
            continue
        best = original.get(node, math.nan)
        for sub, j, kind, needs_zero in edges(f, k):
            other = (sub, j)
            raw = original.get(node, math.nan); nested = original.get(other, math.nan)
            status = "checked"
            if sub in {"log","boxcox_unshifted"} and not positive:
                status = "transform ineligible"
            elif needs_zero and not zero:
                status = "zero shift outside bounds"
            elif not (math.isfinite(raw) and math.isfinite(nested)):
                status = "saved likelihood unavailable"
            loss = nested-raw if status == "checked" else math.nan
            tol = 128*sys.float_info.epsilon*max(1,n,abs(raw),abs(nested)) if status == "checked" else math.nan
            violation = status == "checked" and loss > tol
            if node in own:
                comparisons.append(dict(transform=f,k=k,nested=sub,nested_k=j,comparison=kind,
                    status=status,loglik=raw,nested_loglik=nested,loss=loss,roundoff=tol,violation=int(violation)))
            if status in ("checked","saved likelihood unavailable") and other in lower:
                best = max(best,lower[other]) if math.isfinite(best) else lower[other]
        if math.isfinite(best):
            lower[node] = best
    def choice(ll):
        eligible=sorted((node for node in own if node in ll),key=lambda t:(ORDER.index(t[0]),t[1]))
        return min(eligible,key=lambda z:parameters(z,groups)*math.log(n)-2*ll[z]) if eligible else None
    old, corrected = choice(original), choice(lower)
    label = lambda z: f"{z[0]}:K{z[1]}" if z else ""
    good = [x for x in comparisons if x["status"] == "checked"]
    return comparisons, dict(pairs=len(good),violations=sum(x["violation"] for x in good),
        max_loglik_loss=max((max(0,r["loss"]) for r in good),default=math.nan),
        unavailable=sum(x["status"]=="saved likelihood unavailable" for x in comparisons),
        excluded=sum(x["status"] in ("transform ineligible","zero shift outside bounds") for x in comparisons),
        saved_winner=label(old),winner_using_nested_bounds=label(corrected),
        changes_under_bounds=int(old!=corrected),projected_winner_previously_missing=int(corrected is not None and corrected not in original),
        incomplete_seed_inventory=int(any(r.get("source")=="seed not retained" for r in own.values())))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("data",type=Path);parser.add_argument("out",type=Path)
    a=parser.parse_args(); a.out.mkdir(parents=True,exist_ok=True)
    metadata={r["key"]:r for r in records(a.data/"datasets.csv")}
    saved=defaultdict(lambda:defaultdict(dict))
    for r in records(a.data/"saved-bics.csv"):
        bucket=saved[(r["runtime"],r["key"])][r["method"]]
        node=(r["transform"],int(r["k"]))
        if node in bucket: raise ValueError(f"Duplicate {r}")
        bucket[node]=r
    fits=[]; all_edges=[]; bytype=defaultdict(Counter)
    for (runtime,key),tables in saved.items():
        meta=metadata[key]; idx=number(meta["Index"])
        for method in tables:
            checks, result=audit_fit(meta,method,tables)
            group=meta["group"]; cohort="replay 1-156" if 1<=idx<=156 else "other"
            label=dict(runtime=runtime,method=method,group=group,cohort=cohort,key=key)
            fits.append({**label,**result})
            for row in checks:
                all_edges.append({**label,**row})
                cell=bytype[(runtime,method,group,cohort,row["comparison"])]
                cell[row["status"]]+=1;cell["violations"]+=row["violation"]
    groups=defaultdict(list)
    for row in fits:
        groups[tuple(row[k] for k in ("runtime","method","group","cohort"))].append(row)
    summary=[]
    for (runtime,method,group,cohort),rows in sorted(groups.items()):
        flagged=sum(r["violations"]>0 for r in rows)
        summary.append(dict(runtime=runtime,method=method,group=group,cohort=cohort,
            fits=len(rows),flagged=flagged,rate=flagged/len(rows),
            eligible_checked=sum(r["pairs"] for r in rows),violations=sum(r["violations"] for r in rows),
            unavailable=sum(r["unavailable"] for r in rows),excluded=sum(r["excluded"] for r in rows),
            changes_under_bounds=sum(r["changes_under_bounds"] for r in rows)))
    strata=[dict(runtime=rt,method=m,group=g,cohort=p,comparison=c,checked=z["checked"],
        violations=z["violations"],unavailable=z["saved likelihood unavailable"],
        excluded=z["transform ineligible"]+z["zero shift outside bounds"])
        for (rt,m,g,p,c),z in sorted(bytype.items())]
    write(a.out/"datasets-audited.csv.gz",fits)
    detail=a.out/"pairs"; detail.mkdir(exist_ok=True)
    for i in range(0,len(all_edges),7000):
        write(detail/f"part-{i//7000+1:02d}.csv.gz",all_edges[i:i+7000])
    write(a.out/"violation-summary.csv",summary)
    write(a.out/"by-comparison.csv",strata)


if __name__=="__main__":
    main()
