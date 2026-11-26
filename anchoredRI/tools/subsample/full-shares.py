#!/usr/bin/env python3
"""Pair Task H's full, weighted, draw and median component-1 shares on the scored cohort."""

import argparse
from collections import defaultdict, namedtuple
import csv
import gzip
import math
from pathlib import Path
import runpy
import statistics as stats
import subprocess

HERE = Path(__file__).resolve().parent
old = runpy.run_path(str(HERE / "archived-shares.py"))
finite = old["optional_float"]
mean, quartile, rho = old["mean"], old["quartile"], old["spearman"]
require = old["require"]
# Matched shares are distinct from medians of all valid shares when an endpoint is unavailable.
SummaryRow = namedtuple("SummaryRow", "idx components fraction h share gap realized_gap matched_gap d extreme "
                        "full_gap delta_gap delta_abs_gap delta_d status")
Score = namedtuple("Score", "lower upper z_lower z_upper d")
SCORE_FIELDS = Score._fields


def read(path):
    op = gzip.open if path.suffix == ".gz" else open
    with op(path, "rt", newline="") as stream:
        yield from csv.DictReader(stream)


def indexed(path, key):
    result = {}
    for item in read(path):
        k = item[key]
        require(k not in result, f"Duplicate {key}={k} in {path}")
        result[k] = item
    return result


def score_table(path):
    scores = defaultdict(dict)
    for row in read(path):
        method = row["method"]
        if not (method in {"W3", "W5", "standard_raw", "full3_raw"} or
                method.startswith(("draw_", "matched_", "median_"))):
            continue
        lookup = scores[row["Index"]]
        require(method not in lookup, f"Duplicate score for {row['Index']}, {method}")
        lookup[method] = Score(*(finite(row[name]) for name in SCORE_FIELDS))
    return scores


def sub(left, right):
    return left - right if left is not None and right is not None else None


def gap(share, truth):
    return 100 * (share - truth) if share is not None and truth is not None else None


def absolute(x):
    return abs(x) if x is not None else None


def equal(left, right):
    return (left is None and right is None or left is not None and right is not None and
            math.isclose(left, right, rel_tol=1e-12, abs_tol=1e-13))


def full_method(cap):
    return "standard_raw" if cap == 5 else "full3_raw"


def check_score(s, lower, upper, crp, context):
    require((crp or equal(s.lower, lower)) and equal(s.upper, upper), f"Score endpoints differ: {context}")
    z = (s.z_upper,) if crp else (s.z_lower, s.z_upper)
    expected = stats.mean(abs(v) for v in z) if all(v is not None for v in z) else None
    require(equal(expected, s.d), f"Wrong mean absolute endpoint error: {context}")


def one_case(directory, idx, m, definition, scores):
    fits = list(read(directory / "fits.csv"))
    medians = list(read(directory / "medians.csv"))
    ds = list(read(directory / "decisions.csv"))
    decisions = {(int(r["cap"]), r["seed_set"], r["definition"], r["rule"]): r["status"] for r in ds
                 if r["estimator"] == "median"}
    require(len(fits)==204 and len(medians)==8 and len(ds)==48 and len(decisions)==24 and
            all(r["Index"]==idx for r in fits + medians + ds), f"Incomplete records at {idx}")
    model_keys = {(r["fit_type"],r["cap"],r["seed_set"],r["draw"]) for r in fits}
    expect = {(t,str(c),"","") for t in ("full","W") for c in (3,5)} | {
        ("subsample",str(c),s,str(b)) for c in (3,5) for s in ("A","B") for b in range(1,51)}
    require(model_keys == expect and len({(r["cap"],r["seed_set"],r["definition"]) for r in medians}) == 8,
            f"Incomplete fit or median keys at {idx}")
    full = {int(r["cap"]):r for r in fits if r["fit_type"] == "full"}
    f = float(definition["fractionPathol"])
    bg = float(definition["bg_fraction"])
    h = 1 - f - bg
    counts = old["design_counts"](definition)
    n = int(m["routine_n"])
    realized = counts["phys_n"] / n if sum(counts.values()) == n else None
    crp = definition["Analyte"] == "CRP"
    common = dict(Index=int(idx), analyte=definition["Analyte"], fraction_pathol=f, background_fraction=bg,
                  healthy_share_design=h, background_inclusive_share=1-f, full_routine_n=n, **counts,
                  reconstructed_n=sum(counts.values()), healthy_share_realized_full_vector=realized,
                  routine_hash=m["routine_hash"], anchor_hash=m["anchor_hash"])
    require(definition["Analyte"] == m["Analyte"] and f == float(m["fractionPathol"]) and bg == .001 and
            int(m["anchor_n"]) == 40, f"Design mismatch at {idx}")
    valid_share = lambda x: x is not None and 0 <= x <= 1
    for r in fits + medians:
        is_median = "definition" in r
        cap, bank = int(r["cap"]), r["seed_set"]
        kind = "median" if is_median else r["fit_type"]
        view = r["definition"] if is_median else "all"
        number = "" if is_median else r["draw"]
        draw = kind == "subsample"
        parent = full[cap]
        full_k = int(parent["components"]) if parent["components"] else None
        own_k = None if is_median or not r["components"] else int(r["components"])
        if is_median:
            score_key = f"median_{cap}_{bank}_{view}"
            matched_key = f"matched_{cap}_{bank}_{view}"
            lo, hi = (finite(r[k]) for k in ("lower","upper"))
            matched_lo, matched_hi = (finite(r[k]) for k in ("matched_lower","matched_upper"))
            group_k = full_k
            share = finite(r["share"])
            matched = finite(r["matched_share"])
            share_n, ends_n, paired_n = (int(r[k]) for k in ("share_n","contributors","matched_share_n"))
            require(0 <= paired_n <= min(share_n,ends_n) and max(share_n,ends_n) <= 50 and
                    (valid_share(share) == (share_n>0)) and (valid_share(matched) == (paired_n>0)), f"Median counts at {idx}")
            lrt_flag, lower_flag = "", ""
            status = decisions[cap,bank,view,"majority"]
            full_gate = decisions[cap,bank,view,"full"]
            model_counts = {f"k{k}":int(r[f"k{k}"]) for k in range(1,6)}
            reason = "" if share_n else "no contributing shares"
        else:
            score_key = (full_method(cap) if kind=="full" else f"W{cap}" if kind=="W"
                         else f"draw_{cap}_{bank}_{number}")
            matched_key = score_key
            lo, hi = (finite(r[k]) for k in ("raw_lower", "raw_upper"))
            share = finite(r["share"])
            if share is not None:
                require(valid_share(share) and own_k is not None and 1 <= own_k <= cap and
                        abs(float(r["weight_sum"]) - 1) <= 1e-10 and
                        (own_k != 1 or abs(share-1) <= 1e-10), f"Bad selected routine weights at {idx}")
            else:
                require(bool(r["share_reason"]), f"Missing share reason at {idx}")
            require(r["routine_hash"] == m["routine_hash"] and r["anchor_hash"] == m["anchor_hash"],
                    f"Input hash mismatch at {idx}")
            usable = hi is not None and (crp or lo is not None and lo < hi)
            ends_n, share_n, paired_n = int(usable), int(share is not None), int(usable and share is not None)
            matched = share if paired_n else None
            matched_lo, matched_hi = lo, hi
            group_k = own_k
            status, full_gate = r["first_status"], parent["first_status"]
            lrt_flag, lower_flag, reason = r["lrt_flag"], r["lower_flag"], r["share_reason"]
            model_counts = {f"k{k}":None for k in range(1,6)}
        s, paired_score = scores[score_key], scores[matched_key]
        context = f"{idx}, {score_key}"
        check_score(s, lo, hi, crp, context)
        check_score(paired_score, matched_lo, matched_hi, crp, context)
        d = paired_score.d if paired_n else None
        matched_gap, g = gap(matched,h), gap(share,h)
        zs = (paired_score.z_upper,) if crp else (paired_score.z_lower,paired_score.z_upper)
        extreme = max(abs(z) for z in zs)>5 if d is not None else None
        full_share = finite(parent["share"])
        full_gap = gap(full_share,h)
        full_d = scores[full_method(cap)].d
        abs_change = sub(absolute(g),absolute(full_gap))
        d_change = sub(d,full_d) if full_gap is not None else None
        fields = dict(**common, fit_type=kind, cap=cap, seed_set=bank, draw=number, definition=view,
                      components=own_k, full_components=full_k, group_components=group_k, share=share,
                      fitted_share_n=share_n, endpoint_n=ends_n, matched_share=matched, matched_n=paired_n,
                      lower=lo, upper=hi, raw_d=s.d, matched_lower=matched_lo, matched_upper=matched_hi,
                      z_lower=None if crp else paired_score.z_lower, z_upper=paired_score.z_upper,
                      abs_z_lower=None if crp else absolute(paired_score.z_lower),
                      abs_z_upper=absolute(paired_score.z_upper), matched_d=d, extreme=extreme,
                      gap_pp=g, absolute_gap_pp=absolute(g), matched_gap_pp=matched_gap,
                      realized_full_vector_gap_pp=gap(share,realized), share_reason=reason,
                      full_share=full_share, full_gap_pp=full_gap, full_raw_d=full_d,
                      delta_gap_pp=sub(g,full_gap), delta_abs_gap_pp=abs_change,
                      delta_matched_gap_pp=sub(matched_gap,full_gap),
                      delta_matched_abs_gap_pp=sub(absolute(matched_gap),absolute(full_gap)),
                      delta_matched_d=d_change, status=status, full_gate=full_gate, lrt_flag=lrt_flag,
                      lower_flag=lower_flag, **model_counts)
        summary_row = SummaryRow(int(idx),group_k,f,h,share,g,gap(share,realized),matched_gap,d,extreme,
                                 full_gap,sub(g,full_gap),abs_change,d_change,status)
        yield fields,(kind,cap,bank,view),summary_row


def stats_for(rows):
    gaps = sorted(r.gap for r in rows if r.gap is not None)
    complete = [r for r in rows if r.d is not None and r.matched_gap is not None]
    paired = [r for r in rows if r.delta_abs_gap is not None]
    scored_pairs = [r for r in complete if r.delta_d is not None]
    z = [r.d for r in complete]
    matched_gaps = [r.matched_gap for r in complete]
    revised = [r.realized_gap for r in rows if r.realized_gap is not None]
    extremes = sum(r.extreme for r in complete)
    return dict(n=len(rows),datasets=len({r.idx for r in rows}),valid_shares=len(gaps),
                matched_scorable=len(complete),unscorable=len(rows)-len(complete),released=sum(r.status=="reported" for r in rows),
                mean_fitted_share=mean([r.share for r in rows if r.share is not None]),
                mean_true_share=mean([r.h for r in rows if r.share is not None]), mean_gap_pp=mean(gaps),
                median_gap_pp=quartile(gaps,.5), q25_gap_pp=quartile(gaps,.25),q75_gap_pp=quartile(gaps,.75),
                mean_abs_gap_pp=mean([abs(x) for x in gaps]),mean_matched_d=mean(z),extreme_n=extremes,
                extreme_rate=extremes/len(complete) if complete else None,
                rho_signed=rho(matched_gaps,z),rho_absolute=rho([abs(x) for x in matched_gaps],z),
                paired_shares=len(paired),mean_paired_full_gap_pp=mean([r.full_gap for r in paired]),
                mean_paired_candidate_gap_pp=mean([r.gap for r in paired]),
                mean_delta_gap_pp=mean([r.delta_gap for r in paired]),mean_delta_abs_gap_pp=mean([r.delta_abs_gap for r in paired]),
                paired_scores=len(scored_pairs),mean_paired_delta_matched_d=mean([r.delta_d for r in scored_pairs]),
                mean_paired_score_delta_gap_pp=mean([r.matched_gap-r.full_gap for r in scored_pairs]),
                mean_paired_score_delta_abs_gap_pp=mean([abs(r.matched_gap)-abs(r.full_gap) for r in scored_pairs]),
                realized_full_vector_n=len(revised),mean_realized_full_vector_gap_pp=mean(revised),
                mean_abs_realized_full_vector_gap_pp=mean([abs(x) for x in revised]))


def write(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream,fieldnames=rows[0].keys(),lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def summarize(groups, out, fractions):
    overall, by_fraction, bands = [],[],[]
    for (kind,cap,bank,view), group in sorted(groups.items()):
        key = dict(fit_type=kind,cap=cap,seed_set=bank,definition=view)
        for components in ["all",*range(1,cap+1),"unavailable"]:
            part = group if components == "all" else [r for r in group if r.components == (None if components=="unavailable" else components)]
            dims = dict(**key, group_components=components)
            overall.append(dict(**dims,**stats_for(part)))
            for f in fractions:
                subrows = [r for r in part if r.fraction==f]
                by_fraction.append(dict(**dims,fraction_pathol=f,**stats_for(subrows)))
            for lo,hi,label in old["BANDS"]:
                rr = [r for r in part if r.matched_gap is not None and lo <= abs(r.matched_gap) < hi]
                z = [r for r in rr if r.d is not None]
                n_extreme = sum(r.extreme for r in z)
                bands.append(dict(**dims,abs_matched_gap_band_pp=label,n=len(rr),datasets=len({r.idx for r in rr}),
                                  matched_scorable=len(z),mean_matched_d=mean([r.d for r in z]),extreme_n=n_extreme,
                                  extreme_rate=n_extreme/len(z) if z else None))
    write(out/"by-components.csv",overall)
    write(out/"by-components-and-fraction.csv",by_fraction)
    write(out/"gap-bands.csv",bands)
    return overall


def brief(out, overall, missing_counts):
    body = ["# Component 1 routine shares", "",
            "Each fit is paired with the unweighted full-data fit for the same dataset and component cap. "
            "The table uses the designed physiological fraction, which excludes the 0.1% uniform background. "
            "Positive gap means the model assigned too much routine probability to component 1. "
            "A negative change in absolute gap or endpoint error favors the candidate.", "",
            "| Search | Estimate | Seed set | Definition | Datasets | Mean gap (pp) | Mean change in absolute gap (pp) | Mean change in endpoint error |",
            "| --- | --- | --- | --- | ---: | ---: | ---: | ---: |"]
    def fmt(x): return "NA" if x is None else f"{x:.4f}"
    for r in overall:
        if r["group_components"]!="all" or r["fit_type"]=="subsample" or (r["definition"]=="passing" and r["fit_type"]!="median"):
            continue
        body.append(f"| {r['cap']} | {r['fit_type']} | {r['seed_set'] or '—'} | {r['definition']} | "
                    f"{r['datasets']:,} | {fmt(r['mean_gap_pp'])} | {fmt(r['mean_delta_abs_gap_pp'])} | {fmt(r['mean_paired_delta_matched_d'])} |")
    body.extend(["", "The signed and absolute Spearman correlations use the share and endpoint errors of the same usable fits. "
                 "For a median, the share and limits used for these correlations are medians over the same contributing draws. "
                 "A blank correlation indicates fewer than two distinct values on at least one axis. "
                 "Each draw contributes once within its seed set. Each median contributes once per dataset, cap, bank and definition. "
                 "Individual models are grouped by their selected component count and medians by their full fit’s count. "
                 "The released count for a median uses its majority rule. "
                 "[Component groups](by-components.csv) · [Groups by pathological fraction](by-components-and-fraction.csv) · [Absolute gap bands](gap-bands.csv)", "",
                 f"{missing_counts} datasets had a discrepancy between the reconstructed component counts and observed routine length. "
                 "The realized-fraction sensitivity uses the full routine vector and is reported separately from the primary design-fraction analyses. "
                 "[Count discrepancies](count-discrepancies.csv)", "",
                 "The per-fit export includes requested endpoints, safeguard flags, contribution counts, component counts across draws, "
                 "the parent full fit’s gap and all paired changes. Passing-only individual-model summaries select the draws that passed the first-stage checks. "
                 "[Paired per-fit and per-median records](paired-fits.csv.gz)", "",
                 "Share versus error figures use the component count of the individual fit, or the full fit for a median. "
                 "Columns show the two seed sets separately. [Cap 3, all](share-gap-cap3-all.png) · "
                 "[Cap 3, passing](share-gap-cap3-passing.png) · [Cap 5, all](share-gap-cap5-all.png) · "
                 "[Cap 5, passing](share-gap-cap5-passing.png)", "",
                 "[Prior Task E archived component-share report](https://github.com/amromeo/amromeo.github.io/blob/e5a4577c1fe73869c6a24ce64c8759f190de3197/anchoredRI/tools/subsample/archived-shares/REPORT.md)", ""])
    (out/"REPORT.md").write_text("\n".join(body))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("records","inputs","scores","out"):
        p.add_argument(key,type=Path)
    p.add_argument("--pilot",action="store_true",help="Verify only the twenty previously published pilot datasets")
    args = p.parse_args()
    design = indexed(args.inputs/"design-large.csv","Index")
    defs = indexed(args.inputs/"testsets.csv","Index")
    scores = score_table(args.scores)
    expected = set(indexed(HERE/"pilot-roster.csv","Index")) if args.pilot else set(design)
    require(len(design)==1440 and design.keys()==defs.keys() and expected==scores.keys() and
            {d.name for d in args.records.iterdir() if d.is_dir() and d.name.isdecimal()}==expected,
            "Cohort design, records and score IDs differ")
    args.out.mkdir(parents=True,exist_ok=True)
    groups = defaultdict(list)
    discrepancies = []
    count = 0
    with gzip.open(args.out/"paired-fits.csv.gz","wt",newline="",compresslevel=6) as target:
        for idx in sorted(expected,key=int):
            for fields,key,entry in one_case(args.records/idx,idx,design[idx],defs[idx],scores.pop(idx)):
                if count==0:
                    writer=csv.DictWriter(target,fieldnames=fields.keys(),lineterminator="\n")
                    writer.writeheader()
                writer.writerow(fields)
                groups[key].append(entry)
                if key[0]!="median" and entry.status=="reported":
                    groups[(*key[:3],"passing")].append(entry)
                if fields["fit_type"]=="full" and fields["cap"]==5 and fields["reconstructed_n"]!=fields["full_routine_n"]:
                    discrepancies.append({k:fields[k] for k in ("Index","analyte","full_routine_n","phys_n","left_n","right_n","background_n","reconstructed_n")})
                count+=1
    require(count==len(expected)*(204+8),f"Incomplete paired share record count ({count})")
    with (args.out/"count-discrepancies.csv").open("w",newline="") as stream:
        tab=csv.DictWriter(stream,fieldnames=["Index","analyte","full_routine_n","phys_n","left_n","right_n","background_n","reconstructed_n"],lineterminator="\n")
        tab.writeheader()
        tab.writerows(discrepancies)
    overall=summarize(groups,args.out,sorted({float(d["fractionPathol"]) for d in defs.values()}))
    subprocess.run(["Rscript", str(HERE/"full-share-figure.R"), str(args.out)], check=True)
    brief(args.out,overall,len(discrepancies))
    print(f"Paired {count:,} records across {len(expected):,} distinct datasets")


if __name__=="__main__":
    main()
