"""Publish the synthetic benchmark in text tables with explicit, checked columns."""
import argparse
import csv
import gzip
import itertools
import math
import re
from collections import defaultdict
from pathlib import Path


def rows(path):
    op=gzip.open if path.suffix==".gz" else open
    with op(path,"rt") as f: yield from csv.DictReader(f)


def table(path, fields, records):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,lineterminator="\n",extrasaction="ignore")
        w.writeheader();w.writerows(records)


def saved(src, dataset_manifest=None):
    allfits=list(rows(src/"datasets-audited.csv.gz"))
    parts=defaultdict(list)
    for r in allfits:
        assert re.fullmatch(r"(r-[0-9]+|s-[0-9]{2}-[0-9]{2})",r["key"])
        assert 0<=int(r["violations"])<=int(r["pairs"])
        parts[(r["runtime"],r["method"],r["group"])].append(r)
    if dataset_manifest:
        metadata={r["key"]:r for r in rows(dataset_manifest)}
        strata=defaultdict(list)
        for r in allfits:
            strata[r["runtime"],r["method"],r["group"],metadata[r["key"]]["subtype"],r["cohort"]].append(r)
        detail=[]
        for (rt,m,g,sub,cohort),items in sorted(strata.items()):
            flagged=sum(int(x["violations"])>0 for x in items)
            detail.append(dict(runtime=rt,method=m,group=g,analyte_or_simulation_family=sub,cohort=cohort,
                datasets=len(items),with_violation=flagged,rate=flagged/len(items),
                changes_under_bounds=sum(x["saved_winner"]!=x["winner_using_nested_bounds"] for x in items)))
        table(src/"by-analyte-or-simulation-family.csv",list(detail[0]),detail)
    cols=("key","cohort","pairs","violations","unavailable","excluded","max_loglik_loss","saved_winner","winner_using_nested_bounds","incomplete_seed_inventory")
    for group,items in parts.items():
        name="-".join(group).replace(" ","-")+".csv"
        table(src/"fit-decisions"/name,cols,items)
    pairs=[r for path in sorted((src/"pairs").glob("*.gz")) for r in rows(path) if r["violation"]=="1"]
    assert len(pairs)==sum(int(r["violations"]) for r in allfits)
    cols=("runtime","method","group","key","transform","k","nested","nested_k","comparison","loglik","nested_loglik","loss","roundoff")
    table(src/"violations.csv",cols,pairs)


def improved(src, local, actions):
    allfits=list(rows(src/"candidate-fits.csv.gz"))
    assert len(allfits)==156*15*2
    fits=defaultdict(list)
    for r in allfits: fits[r["Index"],r["runtime"],r["method"]].append(r)
    chosen={}
    for key,items in fits.items():
        ordered=sorted(items,key=lambda r:float(r["bic"]))
        assert ordered[0]["selected"]=="TRUE"
        winner=float(ordered[0]["bic"])
        gap=float(ordered[1]["bic"])-winner
        for rank,r in enumerate(ordered,1):
            r.update(rank=rank,delta_bic=float(r["bic"])-winner,runner_up_gap=gap)
        chosen[key]=ordered[0]
    paired=[]
    for r in rows(src/"paired-replay.csv"):
        for rt,suffix in (("local","local"),("Actions","actions")):
            pick=chosen[r["Index"],rt,r["method"]]
            for field in ("lambda","shift","runner_up_gap"):r[field+"_"+suffix]=pick[field]
        paired.append(r)
    table(src/"paired-replay.csv",list(paired[0]),paired)
    cols=("Index","transform","k","selected","rank","loglik","bic","delta_bic","runner_up_gap","lambda","shift","lower","upper","origin","code","initial_retained","at_bound","starts","nonzero_codes","perturbation_seed")
    for runtime, method in itertools.product(("local","Actions"),("baseline","armA")):
        sel=[r for r in allfits if r["runtime"]==runtime and r["method"]==method]
        assert len(sel)==156*(9 if method=="baseline" else 6)
        table(src/"candidates"/f"{runtime}-{method}.csv",cols,sel)
    details=list(rows(src/"selection-audit.csv.gz"))
    table(src/"original-and-improved-selections.csv",("Index","runtime","method","old_transform","old_k",
            "selected_transform","selected_k","old_model_changed","violations","max_loglik_loss"),details)
    for runtime, folder in (("local",local),("Actions",actions)):
        traces=defaultdict(list)
        for path in sorted(folder.glob("*.gz")):
            for r in rows(path):
                assert r["runtime"]==runtime and 1<=int(r["Index"])<=156
                traces[int(r["Index"]),r["transform"],int(r["k"])].append(r)
        compact=[]
        for (i,tr,k),group in sorted(traces.items()):
            assert [int(r["start"]) for r in group]==list(range(1,len(group)+1))
            assert sum(r["source"].startswith("seeded perturbation") for r in group)==20
            chosen=[r for r in group if r["chosen"]=="TRUE"]
            assert len(chosen)==1
            nonzero=[r for r in group if int(r["convergence"])!=0]
            compact.append(dict(Index=i,transform=tr,k=k,seed=group[0]["seed"],n_starts=len(group),
                chosen_start=chosen[0]["start"],chosen_code=chosen[0]["convergence"],chosen_source=chosen[0]["source"],
                n_nonzero=len(nonzero),nonzero_start_and_code=" ".join(r["start"]+":"+r["convergence"] for r in nonzero)))
        assert len(compact)==156*15
        table(src/"convergence-codes"/f"{runtime}.csv",list(compact[0]),compact)


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("mode",choices=("saved","improved"));p.add_argument("path",type=Path,nargs="+")
    a=p.parse_args()
    (saved if a.mode=="saved" else improved)(*a.path)
