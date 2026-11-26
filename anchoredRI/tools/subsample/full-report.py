"""Run the frozen Task E scoring procedure on saved Task H estimates."""
import argparse
import csv
import gzip
import importlib.util
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
import subprocess
import tempfile

from full_report_tables import analyses, write_rows, num, Score, Check, finite, logical

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("task_h_pilot", HERE / "pilot-report.py")
pilot = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pilot)


def reader(path):
    with (gzip.open if str(path).endswith(".gz") else open)(path, "rt", newline="") as f:
        yield from csv.DictReader(f)


def keyed(path, fields=("Index",)):
    return pilot.by_key(reader(path), *fields)


def checked(status, scope, row=None, route="", votes=None):
    if status not in ("reported", "declined", "withheld"):
        raise ValueError(f"Unexpected release status {status}")
    row = row or {}
    h = lambda k: logical(row.get(k))
    return Check(status, scope, h("lrt_flag"), h("lower_flag"), str(row.get("p", "")), route,
                 *(votes or ("", "", "", "")))


def cap_score(cap):
    return "standard_raw" if cap == "5" else "full3_raw"


def archive(inputs, design, ids, add_ep, add_dec):
    all_ids = set(design)
    e = keyed(inputs / "design-large.csv")
    g = keyed(inputs / "task-g-index.csv")
    if set(g) != {(idx,) for idx in all_ids} or g.keys() != e.keys():
        raise ValueError("Task G and Task E source rosters disagree")
    def verify(idx, row):
        dr = design[idx]
        truth = e[idx,]
        sd = g[idx,]
        if not all(row[k] == truth[k] == sd[k] for k in ("routine_hash", "anchor_hash")):
            raise ValueError(f"Comparator vector mismatch at {idx}")
        if row["Analyte"] != dr["Analyte"] or int(dr["N"]) != 50000:
            raise ValueError(f"Unexpected comparator design at {idx}")
    for file, arm in (("armB-fits.csv", "TaskE_armB"), ("baseline-fits.csv", "TaskE_baseline")):
        fits = keyed(inputs / file)
        if fits.keys() != g.keys():
            raise ValueError(f"Incomplete {arm} archive")
        for idx in ids:
            r = fits[idx,]
            verify(idx, r)
            add_ep(idx, arm, r["lower"], r["upper"])
            flags = dict(lrt_flag=str(num(r["lrt_p"]) < .01), lower_flag=r["lower_flag"], p=r["lrt_p"])
            add_dec(idx, arm, arm, checked(r["report_status"], "archived anchoredRI first stage", flags))
    fits = keyed(inputs / "anchorRefine-fits.csv.gz", ("Index", "arm"))
    selected = {(r["Index"],r["arm"]):r for r in reader(inputs / "anchorRefine-checks.csv.gz")
                if r["panel"] == "guidance" and r["seed_kind"] == "primary"}
    for src, arm in (("baseline", "refineR_native"), ("w1", "anchorRefine_w1"),
                     ("t0.25", "anchorRefine_t0.25"), ("anchorsOnly", "anchors_only")):
        if {idx for idx, a in fits if a == src} != all_ids or {idx for idx, a in selected if a == src} != all_ids:
            raise ValueError(f"Incomplete Task G {src} archive")
        for idx in ids:
            r = fits[idx,src]
            verify(idx, r)
            decision = selected[idx,src]
            if int(decision["seed"]) != 38000000 + int(idx):
                raise ValueError(f"Incorrect primary PIT seed at {idx}")
            valid = logical(r["native_nonfinite"]) is False and r["fit_status"] == "ResultFound"
            stat = "reported" if valid else "withheld"
            if src in ("w1", "t0.25") and valid:
                stat = {"keep":"reported", "decline":"declined", "technical":"withheld"}[decision["status"]]
                for side in ("lower", "upper"):
                    if side == "lower" and r["Analyte"] == "CRP":
                        continue
                    released = num(decision["released_native_"+side])
                    if ((stat == "reported" and not math.isclose(released, num(r["native_"+side]), rel_tol=1e-12, abs_tol=1e-12)) or
                        (stat != "reported" and finite(released))):
                        raise ValueError(f"Stored release differs at {idx}/{src}/{side}")
            flags = dict(p=decision["p"], lrt_flag=str(stat=="declined")) if src in ("w1","t0.25") else {}
            add_ep(idx, arm, r["native_lower"], r["native_upper"])
            add_dec(idx, arm, arm, checked(stat, "archived primary anchor check" if flags else "valid raw estimate", flags))
    return fits


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("records","inputs","output"):
        parser.add_argument(name, type=Path)
    parser.add_argument("--pilot", action="store_true", help="Restrict validation to the frozen 20 pilot records")
    args=parser.parse_args()
    output=args.output
    output.mkdir(parents=True, exist_ok=True)
    design = {r["Index"]:r for r in reader(args.inputs / "testsets.csv") if r["N"]=="50000"}
    canonical = keyed(args.inputs / "task-g-index.csv")
    ids = sorted((r["Index"] for r in reader(HERE/"pilot-roster.csv")) if args.pilot else design, key=int)
    if len(design)!=1440 or len(ids)!=(20 if args.pilot else 1440) or set(ids)!={d.name for d in args.records.iterdir() if d.is_dir() and d.name.isdecimal()}:
        raise ValueError("Record directories do not match the required cohort")
    endpoints, index_states, release_base = defaultdict(dict), defaultdict(dict), {}
    counts, study, flags = Counter(), {}, []
    def add_ep(idx, arm, lower, upper):
        if idx in endpoints[arm]:
            raise ValueError(f"Repeated raw estimate {idx}/{arm}")
        endpoints[arm][idx] = (math.nan if design[idx]["Analyte"]=="CRP" else num(lower), num(upper))
    def add_dec(idx, arm, base, check):
        if idx in index_states[arm] or release_base.setdefault(arm,base) != base:
            raise ValueError(f"Inconsistent checked arm {arm}/{idx}")
        index_states[arm][idx] = check
    taskg_fits = archive(args.inputs,design,ids,add_ep,add_dec)
    for idx in ids:
        base = args.records/idx
        if not (base/"COMPLETE").is_file():
            raise ValueError(f"Incomplete Index {idx}")
        fits = list(reader(base/"fits.csv"))
        m = pilot.by_key(reader(base/"medians.csv"),"cap","seed_set","definition")
        ds = pilot.by_key(reader(base/"decisions.csv"),"cap","seed_set","definition","estimator","rule")
        draw = pilot.by_key(fits,"fit_type","cap","seed_set","draw")
        top = {(r["fit_type"],r["cap"]):r for r in fits if r["fit_type"]!="subsample"}
        full = top["full","5"]
        if len(draw)!=204 or len(top)!=4 or len(m)!=8 or len(ds)!=48:
            raise ValueError(f"Incomplete fit or release set for {idx}")
        meta = canonical[idx,]
        for row in fits:
            if row["Index"]!=idx or any(row[k]!=meta[k] for k in ("routine_hash","anchor_hash")):
                raise ValueError(f"Routine or anchor identity mismatch for {idx}")
        study[idx] = (meta,full,top,m,ds)
        add_ep(idx,"package",full["lower"],full["upper"])
        add_dec(idx,"package","package",checked(full["status"],"standard package",full,full["method"]))
        for cap in ("3","5"):
            f,w=top["full",cap],top["W",cap]
            for stem,r in ((cap_score(cap),f),("W"+cap,w)):
                add_ep(idx,stem,r["raw_lower"],r["raw_upper"])
                add_dec(idx,stem,stem,checked(r["first_status"],"own first stage",r))
            add_dec(idx,"W"+cap+"_full","W"+cap,checked(ds[cap,"A","all","W","full"]["status"],"full fit",f))
            for bank in ("A","B"):
                subs=[draw["subsample",cap,bank,str(i)] for i in range(1,51)]
                votes = (sum(logical(r["lrt_flag"]) is True for r in subs),sum(logical(r["lower_flag"]) is True for r in subs),
                         sum(logical(r["lrt_flag"]) is True and logical(r["lower_flag"]) is True for r in subs),
                         sum(r["first_status"]=="withheld" for r in subs))
                s = ds[cap,bank,"all","W","majority"]
                add_dec(idx,"W"+cap+"_majority_"+bank,"W"+cap,checked(s["status"],"50-draw vote",votes=votes))
                for i,r in enumerate(subs,1):
                    add_ep(idx,f"draw_{cap}_{bank}_{i}",r["raw_lower"],r["raw_upper"])
                for definition in ("all","passing"):
                    key = cap,bank,definition
                    suffix = "_".join(key)
                    z = m[key]
                    for arm,lo,hi in (("median_"+suffix,z["lower"],z["upper"]),
                                      ("matched_"+suffix,z["matched_lower"],z["matched_upper"])):
                        add_ep(idx,arm,lo,hi)
                    med = "median_"+suffix
                    for rule in ("full","majority"):
                        det = ds[key+("median",rule)]
                        state = checked(det["status"],"full fit",f) if rule=="full" else checked(det["status"],"50-draw vote",votes=votes)
                        add_dec(idx,med+"_"+rule,med,state)
                    comb = ds[key+("median","majority_or_package")]
                    method = "policy_"+suffix
                    add_ep(idx,method,comb["lower"],comb["upper"])
                    add_dec(idx,method,method,checked(comb["status"],"majority or full package",route=comb["route"],votes=votes))
                    flags.append(dict(Index=idx,cap=cap,seed_set=bank,definition=definition,
                                      contributors=z["contributors"],share_contributors=z["share_n"],
                                      matched_share_contributors=z["matched_share_n"],votes_decline=z["votes_decline"],
                                      votes_pass=z["votes_pass"],votes_technical=z["votes_technical"],
                                      lrt_votes=votes[0],lower_votes=votes[1],intersection_votes=votes[2],
                                      majority_status=ds[key+("median","majority")]["status"],
                                      full_status=ds[key+("median","full")]["status"],combined_route=comb["route"]))
        for fit_type in ("full","W","subsample"):
            both_caps={(bank,b):(r,draw[fit_type,"5",bank,b]) for (kind,cap,bank,b),r in draw.items()
                       if kind==fit_type and cap=="3"}
            for (bank,b),(c3,c5) in both_caps.items():
                if c3["position_hash"]!=c5["position_hash"] or c3["sampling_seed"]!=c5["sampling_seed"]:
                    raise ValueError(f"Caps are not paired on the same input {idx}/{bank}/{b}")
                counts[fit_type,bank,c3["components"],c5["components"]] += 1

    historical=defaultdict(dict)
    for r in reader(args.inputs / "ribench-outcomes.csv"):
        idx=r["Index"]
        if idx not in study: continue
        arm="historical_"+r["method"]
        if idx in historical[arm]:raise ValueError(f"Duplicate Task E archive Index {idx}")
        historical[arm][idx]=r
        add_ep(idx,arm,r["LRL"],r["URL"])
        stat="reported" if logical(r["reported"]) else ("declined" if logical(r["declined"]) else "withheld")
        add_dec(idx,arm,arm,checked(stat,"Task E archived decision",
                                   dict(p=r["lrt_p"],lrt_flag=str(num(r["lrt_p"])<.01),lower_flag=r["lower_flag"])))
        if r["method"]=="refineR":
            native = taskg_fits[idx,"baseline"]
            for side, field in (("lower","LRL"),("upper","URL")):
                if side=="lower" and r["Analyte"]=="CRP":continue
                if not math.isclose(num(native["historical_"+side]),num(r[field]),rel_tol=1e-10,abs_tol=1e-10):
                    raise ValueError(f"Task E and Task G historical refineR disagree for {idx}/{side}")
    del taskg_fits
    if len(historical)!=6 or any(len(v)!=len(ids) for v in [*endpoints.values(),*index_states.values()]):
        raise ValueError("Incomplete raw or checked method cohort")
    with tempfile.TemporaryDirectory(dir=output) as tmp:
        ear=Path(tmp)/"endpoints.csv"
        write_rows(ear, ["Index","method","lower","upper"],
                   (dict(Index=idx,method=method,lower=v[0],upper=v[1]) for method,ds in endpoints.items() for idx,v in ds.items()))
        scores_csv = Path(tmp)/"scores.csv"
        # This entry point calls native RIbench 1.0.2 computePerfMeasAll, as in Task E.
        subprocess.run(["Rscript",str(HERE/"pilot-score.R"),str(ear),str(args.inputs/"testsets.csv"),str(scores_csv)],check=True)
        scores=defaultdict(dict)
        for r in reader(scores_csv):
            idx,arm=r["Index"],r["method"]
            if idx in scores[arm]:raise ValueError("Repeated score")
            scores[arm][idx]=Score(*endpoints[arm][idx],*(num(r[k]) for k in ("z_lower","z_upper","d","native_cutoff")))
    if any(scores[arm].keys()!=cohort.keys() for arm,cohort in endpoints.items()):
        raise ValueError("RIbench scorer lost endpoint rows")
    del endpoints
    write_rows(output/"endpoint-scores.csv.gz",["Index","method",*Score._fields],
               (dict(Index=i,method=a,**scores[a][i]._asdict()) for a in sorted(scores) for i in ids))
    missing_historical = Counter()
    for arm,hist in historical.items():
        for idx,r in hist.items():
            for k,orig in (("z_lower","zzDev_LRL"),("z_upper","zzDev_URL"),("d","zzDevAbs_Ov"),("native_cutoff","zzDevAbsCutoff_Ov")):
                if k=="z_lower" and r["Analyte"]=="CRP":continue
                v,w=getattr(scores[arm][idx],k),num(r[orig])
                if (math.isnan(v) and math.isnan(w)) or math.isclose(v,w,rel_tol=1e-9,abs_tol=1e-9):continue
                missing_historical[arm]+=1
                if not args.pilot:raise ValueError(f"Complete-cohort Task E score readback failed: {arm}/{idx}/{k}")
    write_rows(output/"gate-records.csv.gz",flags[0].keys(),flags)
    write_rows(output/"cap-component-counts.csv",["fit_type","seed_set","cap3_components","standard_components","fits"],
               (dict(fit_type=k[0],seed_set=k[1],cap3_components=k[2],standard_components=k[3],fits=v) for k,v in sorted(counts.items())))
    # Keep the established split table implementation identical to the pilot.
    paired_arms = ["package","standard_raw","W3","W5"]+[m for m in scores if m.startswith("median_")]
    lookup={(i,m):scores[m][i]._asdict() for m in paired_arms for i in ids}
    split=pilot.first_table(study,lookup,lookup,output)
    analyses(output,ids,design,study,scores,index_states,release_base,split,args.records,args.pilot,counts)
    (output/"report-cohort.json").write_text(json.dumps(dict(
        datasets=len(ids), raw_methods=len(scores), release_methods=len(index_states),
        individual_fit_scores=len(ids)*204, cohort="pilot rehearsal" if args.pilot else "full",
        historical_field_differences_on_authoritative_cohort=dict(missing_historical)),indent=2)+"\n")


if __name__=="__main__":
    main()
