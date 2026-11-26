"""Task H tables calculated from paired RIbench scores and saved check decisions."""
import csv
import gzip
import math
import statistics
from collections import Counter, defaultdict, namedtuple

Score = namedtuple("Score", "lower upper z_lower z_upper d native_cutoff")
Check = namedtuple("Check", "status basis lrt_flag lower_flag p route lrt_votes lower_votes intersection_votes unresolved_votes")
finite = math.isfinite


def num(x):
    return math.nan if x in ("", None, "NA", "NaN") else float(x)


def logical(x):
    return {"TRUE":True,"FALSE":False,"True":True,"False":False,True:True,False:False}.get(x)


def write_rows(path, fields, rows):
    with (gzip.open if str(path).endswith(".gz") else open)(path,"wt",newline="") as stream:
        out=csv.DictWriter(stream,list(fields),lineterminator="\n")
        out.writeheader()
        for r in rows:
            out.writerow({k:"" if isinstance(v,float) and math.isnan(v) else v for k,v in r.items()})


def save(path,rows):
    write_rows(path,rows[0].keys(),rows)


def avg(values):
    v=list(values)
    return statistics.fmean(v) if v else math.nan


def group_ids(design, ids):
    groups=defaultdict(list)
    for i in ids:
        r=design[i]
        stratum="<=20%" if float(r["fractionPathol"])<=.2 else ">20%"
        dist=r["Distribution"]
        for group in (("all","all"),("pathology",stratum),("distribution",dist),
                      ("distribution by pathology",dist+" / "+stratum)):
            groups[group].append(i)
        if abs(float(r["nonp_lambda"])-.432)<1e-12:
            groups["lambda","0.432"].append(i)
    return groups


def at5(d):
    return min(d,5) if finite(d) else 5


def summary(scores,states,base,groups):
    metrics=[];flagrows=[]
    for arm,check in states.items():
        ds=scores[base[arm]]
        for (axis,group),ids in groups.items():
            tally=Counter(check[i].status for i in ids)
            for phase in ("before","after"):
                sample=[i for i in ids if phase=="before" or check[i].status=="reported"]
                v=[ds[i].d for i in sample if finite(ds[i].d)]
                cutoff=[ds[i].native_cutoff for i in sample if finite(ds[i].native_cutoff)]
                scored=set(sample)
                metrics.append(dict(method=arm,raw_estimator=base[arm],axis=axis,group=group,phase=phase,
                     attempted=len(ids),returned=tally["reported"],statistical_declines=tally["declined"],
                     technical=tally["withheld"],return_rate=tally["reported"]/len(ids),
                     decline_rate=tally["declined"]/len(ids),technical_rate=tally["withheld"]/len(ids),
                     evaluated=len(sample),scorable=len(v),missing_among_evaluated=len(sample)-len(v),mean_d=avg(v),
                     native_cutoff_n=len(cutoff),native_cutoff_mean=avg(cutoff),
                     all_attempt_capped_mean=avg(at5(ds[i].d) if i in scored else 5 for i in ids),
                     any_extreme=sum(abs(ds[i].z_lower)>5 or abs(ds[i].z_upper)>5 for i in sample),
                     extreme_lower=sum(abs(ds[i].z_lower)>5 for i in sample),extreme_upper=sum(abs(ds[i].z_upper)>5 for i in sample),
                     d_above_5=sum(ds[i].d>5 for i in sample)))
        for label,test in (("all",lambda d:True),("d<0.25",lambda d:finite(d) and d<.25),
                           ("d>5",lambda d:d>5),("unscored",lambda d:not finite(d))):
            selected=[check[i] for i in groups["all","all"] if test(ds[i].d)]
            flagrows.append(dict(method=arm,error_group=label,attempted=len(selected),
              dataset_lrt_available=sum(c.lrt_flag is not None for c in selected),dataset_lrt_flag=sum(c.lrt_flag is True for c in selected),
              dataset_lower_available=sum(c.lower_flag is not None for c in selected),dataset_lower_flag=sum(c.lower_flag is True for c in selected),
              both_dataset_flags=sum(c.lower_flag is True and c.lrt_flag is True for c in selected),
              individual_lrt_votes=sum(c.lrt_votes for c in selected if c.lrt_votes!=""),
              individual_lower_votes=sum(c.lower_votes for c in selected if c.lower_votes!=""),
              reported=sum(c.status=="reported" for c in selected),declined=sum(c.status=="declined" for c in selected),
              withheld=sum(c.status=="withheld" for c in selected)))
    return metrics,flagrows


def comparisons(states):
    references=("standard_raw","package","TaskE_armB","TaskE_baseline","refineR_native",
                "anchorRefine_w1","anchorRefine_t0.25","anchors_only")
    results=set()
    for a in states:
        if a.startswith("historical_"):continue
        main = a not in references
        for r in references:
            if a!=r and (main or r in ("standard_raw","package")):
                results.add((a,r,"common comparator"))
        if a.startswith(("median_", "policy_")):
            parts=a.split("_")
            if parts[3]=="passing":
                results.add((a,a.replace("_passing","_all"),"passing minus all"))
            if parts[1]=="3":
                results.add((a,a.replace("_3_","_5_",1),"cap 3 minus standard"))
            if parts[2]=="B":
                results.add((a,a.replace("_B_","_A_",1),"bank B minus A"))
            if parts[0]=="median":
                if parts[-1]=="majority":
                    results.add((a,a[:-len("majority")]+"full","majority minus full gate"))
                    results.add((a,"W"+parts[1],"median minus W"))
                else:
                    results.add((a,"W"+parts[1]+"_full","median minus W"))
        if main and (a.startswith("W") or a.startswith("median_") or a.startswith("policy_")):
            results.add((a,"historical_refineR","historical refineR extraction"))
        if a.startswith("W3"):
            results.add((a,a.replace("W3","W5",1),"cap 3 minus standard"))
    results.add(("full3_raw","standard_raw","cap 3 minus standard"))
    for a in ("historical_refineR","historical_anchorsOnly"):
        for r in ("historical_armB","historical_anchoredRI"):
            results.add((a,r,"Task E historical"))
    return sorted(results)


def paired(scores,states,base,groups):
    rows=[];changes=[]
    for a,b,comparison in comparisons(states):
        ac,bc=states[a],states[b];xs,ys=scores[base[a]],scores[base[b]]
        for (axis,group),ids in groups.items():
            ar={i for i in ids if ac[i].status=="reported"}
            br={i for i in ids if bc[i].status=="reported"}
            for phase in ("before","after"):
                both=[i for i in ids if (phase=="before" or i in ar and i in br) and finite(xs[i].d) and finite(ys[i].d)]
                delta=[xs[i].d-ys[i].d for i in both]
                native=[xs[i].d-ys[i].d for i in both if finite(xs[i].native_cutoff) and finite(ys[i].native_cutoff)]
                rows.append(dict(candidate=a,reference=b,contrast=comparison,axis=axis,group=group,phase=phase,
                    assigned=len(ids),candidate_releases=len(ar),reference_releases=len(br),
                    release_rate_difference=(len(ar)-len(br))/len(ids),joint_releases=len(ar & br),pairs=len(both),
                    mean_candidate_d=avg(xs[i].d for i in both),mean_reference_d=avg(ys[i].d for i in both),
                    paired_difference=avg(delta),wins=sum(d<0 for d in delta),ties=sum(d==0 for d in delta),
                    losses=sum(d>0 for d in delta),both_native_cutoff_n=len(native),native_paired_difference=avg(native),
                    paired_all_assignments_capped=avg((at5(xs[i].d) if phase=="before" or i in ar else 5)-
                                                       (at5(ys[i].d) if phase=="before" or i in br else 5) for i in ids)))
        if comparison in ("cap 3 minus standard","bank B minus A","majority minus full gate","passing minus all"):
            outcomes=Counter((ac[i].status,bc[i].status,ac[i].route,bc[i].route) for i in groups["all","all"])
            for (ast,bst,route_a,route_b),n in sorted(outcomes.items()):
                changes.append(dict(candidate=a,reference=b,contrast=comparison,candidate_status=ast,reference_status=bst,
                                    candidate_route=route_a,reference_route=route_b,datasets=n))
    return rows,changes


def extreme_exports(out,ids,design,states,base,scores,records):
    fields=["Index","Analyte","method","raw_estimator","cap","seed_set","draw","status","basis","route","lrt_flag","lower_flag",
            "p","lrt_votes","lower_votes","intersection_votes","unresolved_votes","d","endpoint","estimate","truth","signed_z"]
    highfields=fields[:-4]+["lower","upper","z_lower","z_upper"]
    eh=[];dh=[]
    def add(i,arm,raw,s,c,cap="",bank="",draw=""):
        any_bad=abs(s.z_lower)>5 or abs(s.z_upper)>5
        if not any_bad and not s.d>5:return
        d=design[i]
        common=dict(Index=i,Analyte=d["Analyte"],method=arm,raw_estimator=raw,cap=cap,seed_set=bank,draw=draw,**c._asdict(),d=s.d)
        for side,k in (("lower","GT_LRL"),("upper","GT_URL")):
            z=getattr(s,"z_"+side)
            if abs(z)>5:
                eh.append(dict(**common,endpoint=side,estimate=getattr(s,side),truth=d[k],signed_z=z))
        if s.d>5:
            dh.append(dict(**common,lower=s.lower,upper=s.upper,z_lower=s.z_lower,z_upper=s.z_upper))
    for arm,checks in states.items():
        for i,c in checks.items():
            add(i,arm,base[arm],scores[base[arm]][i],c)
    fit_fields=["Index","Analyte","fit_type","cap","seed_set","draw","components","transform","optimizer_status","optimizer_code",
                "first_status","lrt_flag","lower_flag","p","raw_lower","raw_upper","z_lower","z_upper","d","native_cutoff","any_extreme"]
    def fit_records():
        for i in ids:
            with (records/i/"fits.csv").open(newline="") as f:
                for r in csv.DictReader(f):
                    kind,cap,bank,b=r["fit_type"],r["cap"],r["seed_set"],r["draw"]
                    arm=f"draw_{cap}_{bank}_{b}" if kind=="subsample" else (
                        "W"+cap if kind=="W" else ("standard_raw" if cap=="5" else "full3_raw"))
                    s=scores[arm][i]
                    if kind=="subsample":
                        c=Check(r["first_status"],"own first stage",logical(r["lrt_flag"]),logical(r["lower_flag"]),r["p"],"","","","","")
                        add(i,arm,arm,s,c,cap,bank,b)
                    yield {**{k:r[k] for k in fit_fields if k in r},"Analyte":design[i]["Analyte"],"z_lower":s.z_lower,
                           "z_upper":s.z_upper,"d":s.d,"native_cutoff":s.native_cutoff,"any_extreme":abs(s.z_lower)>5 or abs(s.z_upper)>5}
    write_rows(out/"fit-scores.csv.gz",fit_fields,fit_records())
    write_rows(out/"extreme-endpoints.csv.gz",fields,eh)
    write_rows(out/"datasets-above-five.csv.gz",highfields,dh)


def quantile(sorted_values,p):
    if not sorted_values:return math.nan
    j=(len(sorted_values)-1)*p;ix=int(j)
    return sorted_values[ix] + (j-ix)*(sorted_values[min(ix+1,len(sorted_values)-1)]-sorted_values[ix])


def seed_stability(out,ids,design,study,scores):
    changes=[];data=[];aggregates=defaultdict(list)
    for cap in ("3","5"):
        for definition in ("all","passing"):
            for rule in ("full","majority"):
                counts=Counter()
                for idx in ids:
                    m,dec = study[idx][3:]
                    ar,br=(dec[cap,b,definition,"median",rule]["status"] for b in ("A","B"))
                    am,bm=(m[cap,b,definition] for b in ("A","B"))
                    counts[ar,br]+=1
                    a,b=(scores[f"median_{cap}_{bank}_{definition}"][idx] for bank in ("A","B"))
                    for side in ("lower","upper"):
                        if side=="lower" and design[idx]["Analyte"]=="CRP":continue
                        av,bv=getattr(a,side),getattr(b,side)
                        az,bz=getattr(a,"z_"+side),getattr(b,"z_"+side)
                        movement=bv-av; zdelta=bz-az
                        row=dict(Index=idx,Analyte=design[idx]["Analyte"],unit=design[idx]["Unit"],cap=cap,definition=definition,
                                 rule=rule,endpoint=side,A_status=ar,B_status=br,A_votes_decline=am["votes_decline"],
                                 B_votes_decline=bm["votes_decline"],A_contributors=am["contributors"],B_contributors=bm["contributors"],
                                 A=av,B=bv,B_minus_A=movement,abs_movement=abs(movement),abs_z_movement=abs(zdelta))
                        data.append(row)
                        views=(["raw"] if rule=="majority" else [])+(["released_both_"+rule] if ar==br=="reported" else [])
                        for view in views:
                            for group,measure,val in ((design[idx]["Analyte"],"signed units",movement),
                                                      (design[idx]["Analyte"],"absolute units",abs(movement)),
                                                      ("all analytes","absolute z",abs(zdelta))):
                                aggregates[cap,definition,view,side,group,measure].append((val,idx))
                for (a,b),n in sorted(counts.items()):
                    changes.append(dict(cap=cap,definition=definition,rule=rule,A_status=a,B_status=b,datasets=n))
    summary=[]
    for key, values in sorted(aggregates.items()):
        vals=sorted(x for x,i in values if finite(x));top=max((abs(x) for x,i in values if finite(x)),default=math.nan)
        row=dict(zip(("cap","definition","view","endpoint","group","measure"),key))
        summary.append(dict(**row,eligible=len(values),pairs=len(vals),median=quantile(vals,.5),q90=quantile(vals,.9),
                            q95=quantile(vals,.95),maximum=max(vals) if vals else math.nan,
                            largest_absolute=top,largest_absolute_indices=" ".join(i for x,i in values if abs(x)==top)))
    save(out/"seed-movement.csv.gz",data);save(out/"seed-movement-summary.csv",summary);save(out/"seed-decisions.csv",changes)


def md(x):
    return f"{x:.3f}" if finite(x) else "—"


def analyses(out,ids,design,study,scores,states,base,split,records,is_pilot,component_counts):
    groups=group_ids(design,ids)
    mets,flags=summary(scores,states,base,groups)
    rows,cross=paired(scores,states,base,groups)
    current=[r for r in mets if not r["method"].startswith("historical_")]
    historical=[r for r in mets if r["method"].startswith("historical_")]
    save(out/"task-e-metrics.csv.gz",current)
    save(out/"task-e-historical-reproduction.csv",historical)
    save(out/"paired-contrasts.csv.gz", [r for r in rows if "historical" not in r["contrast"].lower()])
    save(out/"historical-paired-contrasts.csv.gz", [r for r in rows if "historical" in r["contrast"].lower()])
    save(out/"release-changes.csv",cross)
    save(out/"flags-by-error.csv",flags)
    seed_stability(out,ids,design,study,scores);extreme_exports(out,ids,design,states,base,scores,records)
    releases=(dict(Index=idx,Analyte=design[idx]["Analyte"],method=a,raw_estimator=base[a],**states[a][idx]._asdict())
              for a in sorted(states) for idx in ids)
    write_rows(out/"method-decisions.csv.gz",["Index","Analyte","method","raw_estimator",*Check._fields],releases)
    incom=[];routes=Counter();gates=Counter()
    for idx,(meta,f,fc,med,dec) in study.items():
        for key,m in med.items():
            cap,b,df=key
            returned=dec[key+("median","majority_or_package")]
            maj=dec[key+("median","majority")];full=dec[key+("median","full")]
            routes[cap,b,df,returned["route"],returned["status"]]+=1
            gates[cap,b,df,full["status"],maj["status"]]+=1
            if int(m["contributors"])<50:
                incom.append(dict(Index=idx,Analyte=meta["Analyte"],cap=cap,seed_set=b,definition=df,
                                   contributors=m["contributors"],votes_decline=m["votes_decline"],technical_votes=m["votes_technical"],
                                   majority_status=maj["status"],full_status=full["status"]))
    write_rows(out/"median-contributors-below-50.csv.gz",["Index","Analyte","cap","seed_set","definition","contributors",
                                                        "votes_decline","technical_votes","majority_status","full_status"],incom)
    write_rows(out/"policy-routes.csv",["cap","seed_set","definition","route","status","datasets"],
               (dict(cap=k[0],seed_set=k[1],definition=k[2],route=k[3],status=k[4],datasets=v) for k,v in sorted(routes.items())))
    write_rows(out/"full-vs-majority.csv",["cap","seed_set","definition","full_status","majority_status","datasets"],
               (dict(cap=k[0],seed_set=k[1],definition=k[2],full_status=k[3],majority_status=k[4],datasets=v) for k,v in sorted(gates.items())))
    component_changes = []
    for kind,bank in sorted({(k[0],k[1]) for k in component_counts}):
        path = [(c3,c5,n) for (ft,b,c3,c5),n in component_counts.items() if (ft,b)==(kind,bank)]
        component_changes.append(dict(fit_type=kind,seed_set=bank,attempts=sum(n for _,_,n in path),
            both_selected=sum(n for c3,c5,n in path if c3 and c5),
            selected_count_changed=sum(n for c3,c5,n in path if c3 and c5 and c3!=c5),
            selected_count_equal=sum(n for c3,c5,n in path if c3 and c5 and c3==c5),
            missing_either_count=sum(n for c3,c5,n in path if not c3 or not c5)))
    save(out/"cap-selection-changes.csv",component_changes)
    stage_n=Counter(study[i][1]["first_status"] for i in ids)
    report=["# Task H subsampling results" if not is_pilot else "# Task H pilot rehearsal of the full-cohort report","",
            f"The cohort contains {len(ids):,} RIbench datasets at nominal N = 50,000. "
            "Each dataset contributes one equally weighted comparison. Each subsample draws 5,000 routine observations and retains its 40 original anchors.",
            "", "| Standard full first stage | Cap | Eligible | Majority releases | W releases | Common triples | Mean error change, median | Mean error change, W | Extreme intervals, reference / median / W |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---|"]
    for r in split:
        if r["bank"]!="A" or r["definition"]!="all":continue
        report.append(f"| {'Release' if r['panel']=='A' else 'Decline'} | {r['cap']} | {r['eligible']} | "
                      f"{r['majority_release']} | {r['W_release']} | {r['common']} | {md(r['median_paired_delta_d'])} | "
                      f"{md(r['W_paired_delta_d'])} | {r['ref_extreme_any']} / {r['median_extreme_any']} / {r['W_extreme_any']} |")
    report += ["", "Seed set A uses the coordinatewise medians of all usable raw draws. Panel A compares with the five-component-cap full "
               "first-stage estimate. In panel B the median must pass its majority check and the full package must return a fallback. "
               f"{stage_n['withheld']:,} full first stages had technical withholdings. The endpoint columns count intervals with any requested |z| > 5. "
               "W contributes its raw limits to the common triples. Its own two checks determine its release count. "
               "[All split panels and W’s standalone releases](paired-summary.csv) · [Each pair](paired-index.csv) · [Limit changes by analyte](paired-limits-by-analyte.csv)"]
    if is_pilot: report+=["","This rehearsal uses only the 20 preselected pilot datasets. Analyte-wide scoring decisions use just those 20."]
    overall={(r["method"],r["phase"]):r for r in mets if r["axis"]=="all"}
    table=["", "## Accuracy and release", "",
           "| Estimator and gate | Returned / assigned | Declined | Technical | Pre-check mean error | Released mean error | Released cutoff mean (n) | All-assignment score after checking | Released extreme intervals |",
           "|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    labels=[("standard_raw","Full, first stage"),("package","Full, package"),("W5","W, cap 5, own checks"),
            ("median_5_A_all_full","Median, cap 5, full gate"),("median_5_A_all_majority","Median, cap 5, majority"),
            ("policy_5_A_all","Majority median or package, cap 5"),("full3_raw","Full, cap 3, first stage"),
            ("W3","W, cap 3, own checks"),("median_3_A_all_full","Median, cap 3, full gate"),
            ("median_3_A_all_majority","Median, cap 3, majority"),("policy_3_A_all","Majority median or package, cap 3"),
            ("TaskE_armB","Archived Task E Arm B"),("TaskE_baseline","Archived Task E baseline"),
            ("refineR_native","refineR, native limits"),("anchorRefine_w1","anchorRefine, w = 1"),
            ("anchorRefine_t0.25","anchorRefine, t = 0.25"),("anchors_only","Anchors only")]
    for method,label in labels:
        raw=overall[method,"before"];q=overall[method,"after"]
        table.append(f"| {label} | {q['returned']:,} / {len(ids):,} | {q['statistical_declines']} | {q['technical']} | "
                     f"{md(raw['mean_d'])} | {md(q['mean_d'])} | {md(q['native_cutoff_mean'])} ({q['native_cutoff_n']:,}) | "
                     f"{md(q['all_attempt_capped_mean'])} | {q['any_extreme']} |")
    report+=table+["","The dataset error is the mean absolute deviation on RIbench’s z scale over requested endpoints, with the upper limit alone for CRP. "
                   "The native cutoff excludes errors above 5. The all-assignment score caps errors at 5 and assigns 5 to a nonreturn. "
                   "For the package and combined policy, the score uses the actual returned endpoints over the entire cohort. Native refineR "
                   "limits are primary. Archived anchorRefine rows use the primary Stage 1 check. The historical table reproduces Task E’s original extraction. "
                   "[All strata and filters](task-e-metrics.csv.gz) · [Historical Task E results](task-e-historical-reproduction.csv) · [Historical paired differences](historical-paired-contrasts.csv.gz) · [Paired differences, wins, ties, losses and coverage](paired-contrasts.csv.gz)",
                   "", "## Safeguards, component caps and seed stability", ""]
    for cap in ("5","3"):
        q=[r for r in split if r["bank"]=="A" and r["definition"]=="all" and r["cap"]==cap and r["panel"]=="B"][0]
        all_routes=Counter()
        for (c,s,d,route,status),n in routes.items():
            if (c,s,d)==(cap,"A","all"):all_routes[route]+=n
        ids_changed=[i for i in ids if states[f"median_{cap}_A_all_majority"][i].status != states[f"median_{cap}_B_all_majority"][i].status]
        report+= [f"- At cap {cap}, among the {q['eligible']:,} standard full-data statistical declines, the majority released "
                  f"{q['majority_only']:,} when the package withheld and the package released {q['package_only']:,} when the majority withheld. "
                  f"The combined policy returned {all_routes['median']:,} medians, {all_routes['full_first_stage']:,} full first-stage estimates "
                  f"and {all_routes['full_fallback']:,} released fallbacks, with {all_routes['full_nonreturn']:,} nonreturns. "
                  f"{len(ids_changed):,} primary majority decisions changed between seed sets."]
    report += ["", "Each model selection uses the same ordered routine values, anchors and positioning seeds across the two caps. "
               "The movement tables give B-minus-A changes by side in analyte units and absolute changes on the z scale. They include passing-fits-only medians, "
               "the largest movements by Index, and the subset released under both seeds. "
               "[Seed movement summaries](seed-movement-summary.csv) · [Per-Index movement and votes](seed-movement.csv.gz) · [Seed decision changes](seed-decisions.csv) · "
               "[Cap and gate decision changes](release-changes.csv) · [Changed component counts](cap-selection-changes.csv) · [Component-count pairs](cap-component-counts.csv)",
               "", "The vote includes all 50 draws. Technical votes can leave the majority unresolved. A statistically declined fit can contribute valid raw limits "
               "to the primary median. A missing passing-only median takes the saved full-package route. "
               "[Vote and flag records](gate-records.csv.gz) · [Missing contributors](median-contributors-below-50.csv.gz) · [Policy routes](policy-routes.csv) · "
               "[Full versus majority](full-vs-majority.csv) · [Flags by error](flags-by-error.csv)",
               "", "## Endpoint and fit records", "",
               "Each endpoint and fit uses the frozen RIbench 1.0.2 scorer. For individual draws, each cap, seed set and draw number defines one "
               "raw method on the cohort. All scores precede release masking. The extreme lists retain declined raw intervals and identify each gate. "
               "[Endpoint scores](endpoint-scores.csv.gz) · [Fit scores](fit-scores.csv.gz) · [Method decisions](method-decisions.csv.gz) · "
               "[Every extreme endpoint](extreme-endpoints.csv.gz) · [Datasets with mean error above 5](datasets-above-five.csv.gz)",
               "", "The physiological-share analysis pairs each fit with the unweighted fit at the same cap. The optimizer analysis compares iteration-limit "
               "and line-search exits with converged fits. [Component shares](shares/REPORT.md) · [Paired share and error changes](shares/paired-fits.csv.gz) · [Optimizer exits](optimizer/OPTIMIZER.md)", "",
               "RIbench informed the subsampling proposal. Adoption would require independent confirmation and a prespecified sample-size rule for 5,000–50,000 routine results.", ""]
    (out/"REPORT.md").write_text("\n".join(report))
