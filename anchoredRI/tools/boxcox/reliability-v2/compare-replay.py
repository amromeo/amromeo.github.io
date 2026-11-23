"""Apply only the preregistered reliability gates to Index 1 through 156."""
import argparse
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path
from nesting import audit_fit, number, records, write, ORDER


def load(path, runtime):
    if not (path/"DONE.txt").is_file(): raise ValueError(f"{runtime} replay is incomplete")
    source={r["file"]:r["sha256"] for r in records(path/"sources.csv")}
    meta={int(r["Index"]):r for r in records(path/"datasets.csv")}
    assert set(meta)==set(range(1,157)) and len(source)==8
    result=defaultdict(lambda:defaultdict(dict))
    for r in records(path/"candidates.csv"):
        i=int(r["Index"]); assert r["runtime"]==runtime and r["key"]==f"r-{i}"
        node=(r["transform"],int(r["k"])); group=result[i][r["method"]]
        assert node not in group; group[node]=r
    assert set(result)==set(meta)
    for i, methods in result.items():
        assert len(methods["baseline"])==9 and len(methods["armA"])==6
        for rows in methods.values():
            assert sum(row["selected"]=="TRUE" for row in rows.values())==1
            for row in rows.values():
                assert int(row["starts"])>=21 and len(row["routine_hash"])==64 and len(row["anchor_hash"])==64
    return source,meta,result


def original_choices(saved, archive):
    inventory=defaultdict(list)
    for r in records(saved):
        i=r["key"].split("-")
        if i[0]!="r" or int(i[1])>156 or r["method"] not in ("baseline","armA") or r["runtime"] not in ("local","Actions"):continue
        if r["runtime"]=="Actions" and r["method"]=="baseline":continue
        if math.isfinite(number(r["bic"])):
            inventory[int(i[1]),r["runtime"],r["method"]].append(r)
    choices={}
    for key, rows in inventory.items():
        r=min(rows,key=lambda r:(number(r["bic"]),ORDER.index(r["transform"]),int(r["k"])))
        choices[key]=(r["transform"],int(r["k"]))
    for r in records(archive):
        if int(r["Index"])<=156:
            choices[int(r["Index"]),"Actions","baseline"]=(r["archived_transform"],int(r["archived_components"]))
    return choices


def main():
    p=argparse.ArgumentParser()
    for arg in ("local","actions","saved_bics","archive_checks","out"):p.add_argument(arg,type=Path)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    local=load(a.local,"local"); actions=load(a.actions,"Actions")
    assert local[0]==actions[0],"fitting source mismatch"
    structural=("key","Index","group","subtype","routine_n","anchor_n","n","groups","positive","zero_shift")
    for i in range(1,157):
        assert all(local[1][i][k]==actions[1][i][k] for k in structural),"input metadata mismatch"
    original=original_choices(a.saved_bics,a.archive_checks)
    assert len(original)==624
    out=[]; details=[]; candidates=[]; nesting=[]; totals=defaultdict(Counter)
    for i in range(1,157):
        selections={}
        for name,(_,meta,allfits) in (("local",local),("Actions",actions)):
            fits=allfits[i]
            ids={(r["routine_hash"],r["anchor_hash"]) for rows in fits.values() for r in rows.values()}
            assert len(ids)==1
            if name=="local": local_ids=ids
            else: assert local_ids==ids
            for m, rows in fits.items():
                node,r=next((z,row) for z,row in rows.items() if row["selected"]=="TRUE")
                selections[name,m]=(node,r)
                ch,result=audit_fit(meta[i],m,fits)
                for row in ch:
                    nesting.append(dict(Index=i,runtime=name,method=m,**row))
                assert node==(result["saved_winner"].split(":")[0],int(result["saved_winner"].split(":K")[1]))
                assert len({row["perturbation_seed"] for row in rows.values()})==len(rows)
                cnt=totals[m,name];cnt["nesting_violations"]+=result["violations"];cnt["checked"]+=result["pairs"]
                cnt["unavailable"]+=result["unavailable"]
                cnt["changed_from_old"]+=node!=original[i,name,m]
                cnt["selected_nonzero"]+=r["code"]!="0"
                cnt["initial_winner"]+=r["initial_retained"]=="TRUE"
                cnt["all_starts"]+=sum(int(x["starts"]) for x in rows.values())
                cnt["nonzero_starts"]+=sum(int(x["nonzero_codes"]) for x in rows.values())
                for x in rows.values():candidates.append(x)
                details.append(dict(Index=i,runtime=name,method=m,old_transform=original[i,name,m][0],
                    old_k=original[i,name,m][1],selected_transform=node[0],selected_k=node[1],
                    old_model_changed=int(node!=original[i,name,m]),**result))
        for m in ("baseline","armA"):
            ln,lr=selections["local",m];rn,rr=selections["Actions",m]
            width=number(rr["upper"])-number(rr["lower"])
            dl=abs(number(lr["lower"])-number(rr["lower"]));du=abs(number(lr["upper"])-number(rr["upper"]))
            finite=all(map(math.isfinite,(width,dl,du))) and width>0
            worst=max(dl,du) if finite else math.nan
            out.append(dict(Index=i,method=m,local_transform=ln[0],local_k=ln[1],actions_transform=rn[0],actions_k=rn[1],
                model_match=int(ln==rn),lower_local=lr["lower"],upper_local=lr["upper"],lower_actions=rr["lower"],
                upper_actions=rr["upper"],width_actions=width,lower_difference=dl,upper_difference=du,max_difference=worst,
                fraction_of_width=worst/width if finite else math.nan,finite_pair=int(finite),
                within_0_5pct=int(finite and dl<=.005*width and du<=.005*width),
                routine_hash=lr["routine_hash"],anchor_hash=lr["anchor_hash"]))
    summary=[]
    for m in ("baseline","armA"):
        rows=[r for r in out if r["method"]==m];assert len(rows)==156
        model=sum(r["model_match"] for r in rows);lim=sum(r["within_0_5pct"] for r in rows)
        loc=totals[m,"local"];act=totals[m,"Actions"]
        ok=model>=152 and lim>=149 and all(c["nesting_violations"]==c["unavailable"]==0 for c in (loc,act))
        usable=[r for r in rows if r["finite_pair"]]
        summary.append(dict(method=m,attempts=156,model_matches=model,model_required=152,
             finite_limit_pairs=len(usable),within_tolerance=lim,limits_required=149,
             local_nesting_violations=loc["nesting_violations"],actions_nesting_violations=act["nesting_violations"],
             local_comparisons=loc["checked"],actions_comparisons=act["checked"],
             unavailable_local=loc["unavailable"],unavailable_actions=act["unavailable"],
             max_absolute_difference=max((r["max_difference"] for r in usable),default=math.nan),
             max_fraction_of_width=max((r["fraction_of_width"] for r in usable),default=math.nan),passed=int(ok)))
    write(a.out/"gates.csv",summary);write(a.out/"paired-replay.csv",out);write(a.out/"selection-audit.csv.gz",details)
    write(a.out/"candidate-fits.csv.gz",candidates);write(a.out/"nesting.csv.gz",nesting)
    write(a.out/"convergence.csv",[dict(method=m,runtime=t,**v) for (m,t),v in sorted(totals.items())])
    write(a.out/"source-hashes.csv",[dict(file=f,sha256=h) for f,h in local[0].items()])


if __name__=="__main__":main()
