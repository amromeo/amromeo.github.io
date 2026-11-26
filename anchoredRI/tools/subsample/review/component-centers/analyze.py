# SPDX-License-Identifier: GPL-3.0-only
"""Classify the extra component centers in archived five-component Task E Arm B fits."""
import argparse
import csv
import hashlib
import io
import math
from collections import defaultdict
from itertools import groupby
from pathlib import Path
from statistics import NormalDist, correlation, fmean, median
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from tail.decompose import healthy_cdf, region_limits

DESIGN_SHA = "177cd392896cb300bf4aaa35a959e46bb3ee48d510beff8cf9957f6643a66e39"
PARAM_SHA = "e77e0fa6e6329f168ba01227069819eedc7bd436e987b3a20437d22cbde7958c"
FIT_SHA = "87dce1dadba2917bd2bf160beab1933bd9d58966291ccae4b8dd6011c42fdef7"
SCORE_SHA = "05badc90db1757901244616f2364fe591926b69e7d7135f6f3b640debda72716"
STD_NORMAL = NormalDist()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def source(path, checksum):
    data = path.read_bytes()
    require(hashlib.sha256(data).hexdigest() == checksum, f"Source SHA-256 differs for {path}")
    return csv.DictReader(io.StringIO(data.decode("utf-8-sig")))


def index_table(rows):
    table = {}
    for row in rows:
        j = int(row["Index"])
        require(j not in table, f"Duplicate Index {j}")
        table[j] = row
    return table


def original_scale(y, transform, shift):
    if transform == "identity":
        # Task E's parameter export has already added one to the internal identity location.
        return y
    if transform == "log" and shift == 0 or transform == "shifted_log":
        return shift + math.exp(y)
    raise ValueError(f"Unexpected fitted transform {transform} with shift {shift}")


def midranks(values):
    ranks = [0.] * len(values)
    order = sorted(range(len(values)), key=values.__getitem__)
    for _, group in groupby(enumerate(order, start=1), key=lambda p: values[p[1]]):
        ties = list(group)
        rank = (ties[0][0] + ties[-1][0]) / 2
        for _, index in ties:
            ranks[index] = rank
    return ranks


def spearman(x, y):
    if len(set(x)) < 2 or len(set(y)) < 2:
        return None
    return correlation(midranks(x), midranks(y))


def records(sources):
    designs = index_table(source(sources / "testsets.csv", DESIGN_SHA))
    fits = index_table(source(sources / "armB-fits.csv", FIT_SHA))
    scores = index_table(row for row in source(sources / "ribench-outcomes.csv", SCORE_SHA)
                         if row["method"] == "armB")
    params = defaultdict(list)
    for row in source(sources / "selected-component-parameters.csv", PARAM_SHA):
        if row["method"] == "armB":
            params[row["dataset_id"]].append(row)
    require(len(fits) == 1440 and designs.keys() == fits.keys() == scores.keys(), "Cohort mismatch")
    require(set(params) == {f"r-{j}" for j in fits}, "Parameter cohort mismatch")
    datasets, components = [], []
    for j, fit in sorted(fits.items()):
        design, score = designs[j], scores[j]
        require(fit["method"] == "armB" and fit["max_components"] == "5" and
                fit["components"] == score["components"], f"Fit/score mismatch at {j}")
        if fit["components"] != "5":
            continue
        dataset = f"r-{j}"
        analyte = design["Analyte"]
        require(fit["dataset_id"] == dataset and fit["Analyte"] == score["Analyte"] == analyte and
                fit["N"] == design["N"] == "50000" and fit["ribench_hashcode"] == design["HashCode"],
                f"Bad fit/design join at {j}")
        for name in ("GT_LRL", "GT_URL"):
            if name == "GT_LRL" and analyte == "CRP":
                continue
            require(float(design[name]) == float(score[name]), f"Truth mismatch {j}, {name}")
        parts = sorted(params[dataset], key=lambda r: int(r["component"]))
        require([int(r["component"]) for r in parts] == list(range(1, 6)), f"Component list at {j}")
        lo, hi = region_limits(design, rounded=False)
        require(healthy_cdf(lo, design) < .25 and healthy_cdf(hi, design) > .75,
                f"Healthy quartiles beyond reference limits at {j}")
        h = 1 - float(design["fractionPathol"]) - float(design["bg_fraction"])
        w1 = float(parts[0]["routine_weight"])
        err = float(score["zzDevAbs_Ov"])
        fields = ("upper",) if analyte == "CRP" else ("lower", "upper")
        signed = []
        for side in fields:
            endpoint = "URL" if side == "upper" else "LRL"
            signed.append(float(score[f"zzDev_{endpoint}"]))
            stored = float(fit[side])
            require(math.isclose(float(score[endpoint]), stored, rel_tol=1e-12, abs_tol=1e-12),
                    f"Pre-safeguard endpoint join at {j}, {side}")
            prob = float(fit[f"{side}_probability"])
            m, sd = (float(parts[0][field]) for field in ("transformed_mean", "transformed_sd"))
            predicted = original_scale(m + sd * STD_NORMAL.inv_cdf(prob), fit["transform"], float(fit["shift"]))
            require(math.isclose(predicted, stored, rel_tol=1e-10, abs_tol=1e-10),
                    f"Inverse transform disagrees with the archived {side} endpoint at {j}")
        require(math.isclose(err, fmean(abs(z) for z in signed), rel_tol=1e-12, abs_tol=1e-14),
                f"Pre-safeguard score mismatch at {j}")
        base = dict(Index=j, analyte=analyte, pathological_fraction=float(design["fractionPathol"]),
                    report_status=fit["report_status"], fitted_transform=fit["transform"],
                    healthy_pct=100*h, component1_pct=100*w1, shortfall_pp=100*(h-w1), d=err)
        fitrows = []
        for k, part in enumerate(parts, 1):
            w, shift, m, sd = (float(part[key]) for key in
                               ("routine_weight", "shift", "transformed_mean", "transformed_sd"))
            require(part["transform"] == fit["transform"] and shift == float(fit["shift"]) and
                    part["setting"] == "pooled" and (part["screened_component"] == "TRUE") == (k == 1)
                    and 0 <= w <= 1 and math.isfinite(sd) and sd > 0, f"Component fields at {j}:{k}")
            if k == 1:
                continue
            center = original_scale(m, part["transform"], shift)
            q = healthy_cdf(center, design)
            require(math.isfinite(center) and 0 <= q <= 1, f"Center/CDF outside range at {j}:{k}")
            fitrows.append(dict(Index=j, component=k, shift=shift, transformed_mean=m, transformed_sd=sd,
                                routine_weight=w, center_original=center, true_healthy_percentile=100*q,
                                in_reference=lo <= center <= hi, in_central50=.25 <= q <= .75))
        require(abs(w1 + sum(r["routine_weight"] for r in fitrows) - 1) <= 1e-10, f"Weights at {j}")
        require(all(not r["in_central50"] or r["in_reference"] for r in fitrows), f"Region mismatch {j}")
        components.extend(fitrows)
        for key in ("reference", "central50"):
            inside = [r for r in fitrows if r["in_" + key]]
            base[key + "_n"] = len(inside)
            base[key + "_weight_pct"] = 100 * sum(r["routine_weight"] for r in inside)
        datasets.append(base)
    require(len(datasets) == 941 and len(components) == 3764, "Incomplete five-component cohort")
    return datasets, components


def summaries(datasets):
    groups = defaultdict(list)
    for row in datasets:
        for axis, value in (("all", "all"), ("analyte", row["analyte"]),
                            ("pathological_fraction", row["pathological_fraction"]),
                            ("first_stage", row["report_status"])):
            groups[axis, value].append(row)
    out, counts = [], []
    for (axis, value), rows in sorted(groups.items()):
        for region in ("reference", "central50"):
            with_center = [r for r in rows if r[region + "_n"] > 0]
            without = [r for r in rows if r[region + "_n"] == 0]
            shares = [r[region + "_weight_pct"] for r in rows]
            cond = [r[region + "_weight_pct"] for r in with_center]
            rec = dict(axis=axis, value=value, region=region, fits=len(rows), fits_with_center=len(with_center),
                       percent_fits_with_center=100*len(with_center)/len(rows),
                       centers_inside=sum(r[region + "_n"] for r in rows), eligible_centers=4*len(rows),
                       mean_sum_weight_pct=fmean(shares), median_sum_weight_pct=median(shares),
                       mean_sum_weight_pct_if_present=fmean(cond) if cond else None,
                       median_sum_weight_pct_if_present=median(cond) if cond else None,
                       rho_weight_shortfall=spearman(shares, [r["shortfall_pp"] for r in rows]),
                       rho_weight_d=spearman(shares, [r["d"] for r in rows]))
            for label, data in (("with", with_center), ("without", without)):
                rec[label + "_mean_shortfall_pp"] = fmean(r["shortfall_pp"] for r in data) if data else None
                rec[label + "_mean_d"] = fmean(r["d"] for r in data) if data else None
            out.append(rec)
            if axis == "all":
                for n in range(5):
                    data = [r for r in rows if r[region + "_n"] == n]
                    counts.append(dict(region=region, centers_inside=n, fits=len(data),
                                       mean_weight_pct=fmean(r[region + "_weight_pct"] for r in data) if data else None,
                                       mean_shortfall_pp=fmean(r["shortfall_pp"] for r in data) if data else None,
                                       mean_d=fmean(r["d"] for r in data) if data else None))
    return out, counts


def write(path, rows):
    with path.open("w", newline="") as target:
        table = csv.DictWriter(target, fieldnames=rows[0].keys(), lineterminator="\n")
        table.writeheader()
        table.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sources", type=Path, help="Directory containing the four pinned CSV inputs")
    parser.add_argument("--out", type=Path, default=HERE)
    args = parser.parse_args()
    datasets, components = records(args.sources)
    summary, counts = summaries(datasets)
    args.out.mkdir(parents=True, exist_ok=True)
    for name, rows in (("datasets", datasets), ("components", components), ("summary", summary),
                       ("by-count", counts)):
        write(args.out / (name + ".csv"), rows)


if __name__ == "__main__":
    main()
