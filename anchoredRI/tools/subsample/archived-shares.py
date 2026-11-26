"""Summarize the saved Task E mixtures at N = 50,000. Requires Python >= 3.10.

Run with --anchorrefine-repo PATH, where PATH is a checkout of amromeo/anchorRefine.
Both Git repositories must contain the pinned revisions below.
"""
import argparse
import csv
import hashlib
import io
import itertools
import json
import math
import statistics as stats
import subprocess
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
TASK_E = "8bcd49c3c7977c50982524b9dbc48d16bdfa5914"
STAGE_3 = "638e2d3500faa01d4e2816895a95d7a6685567a0"
TASK_E_DIR = "anchoredRI/tools/boxcox/accuracy/large/results/"
TASK_E_FILES = {
    "selected-component-parameters.csv": "e77e0fa6e6329f168ba01227069819eedc7bd436e987b3a20437d22cbde7958c",
    "armB-fits.csv": "87dce1dadba2917bd2bf160beab1933bd9d58966291ccae4b8dd6011c42fdef7",
    "baseline-fits.csv": "cab045a054d0651dbe0b20c8e5d8602957e521284d2b98b21c4256bc26fe960a",
    "ribench-outcomes.csv": "05badc90db1757901244616f2364fe591926b69e7d7135f6f3b640debda72716",
}
INPUT_FILES = {
    "stage1/large/testsets.csv": "177cd392896cb300bf4aaa35a959e46bb3ee48d510beff8cf9957f6643a66e39",
    "stage2/ribench/design-large.csv": "e7b3b9b0673211f84534c11f2125a00eaad3295b8180772e419d6de7c0faa267",
    "vendor/RIbench_1.0.2.tar.gz": "5703783160e2ed6bef1891d44c1401cfa0303820c46eb0f3e56f213c2ae3dd35",
}
CAPS = {"armB": 5, "anchoredRI": 3}
BANDS = ((0, 5, "<5"), (5, 10, "5 to <10"), (10, 20, "10 to <20"), (20, math.inf, ">=20"))


def require(condition, detail):
    if not condition:
        raise ValueError(detail)


def load_files(repo, revision, prefix, hashes):
    result = {}
    for name, expected in hashes.items():
        data = subprocess.check_output(["git", "-C", str(repo), "show", revision + ":" + prefix + name])
        require(hashlib.sha256(data).hexdigest() == expected, f"Source SHA-256 differs for {name}")
        if name.endswith(".csv"):
            result[name] = list(csv.DictReader(io.StringIO(data.decode("utf-8-sig"))))
    return result


def indexed(rows, fields):
    lookup = {}
    for row in rows:
        key = tuple(row[k] for k in fields)
        require(key not in lookup, f"Duplicate key {key}")
        lookup[key] = row
    return lookup


def optional_float(value):
    if value in ("", "NA", "NaN"):
        return None
    number = float(value)
    return number if math.isfinite(number) else None


def design_counts(design):
    """Mirror R's positive size-to-integer conversion, keeping the generator's operation order."""
    n = int(design["N"])
    fraction, background = float(design["fractionPathol"]), float(design["bg_fraction"])
    left, right, total = map(float, (design["left_ratio"], design["right_ratio"], design["ratio_sum"]))
    return dict(phys_n=int(n * (1 - fraction - background)),
                left_n=0 if design["Analyte"] == "CRP" else int(n * (fraction * (left / total))),
                right_n=int(n * (fraction * (right / total))), background_n=int(n * background))


def make_records(tables, inputs):
    definition = indexed(inputs["stage1/large/testsets.csv"], ["Index"])
    design = indexed(inputs["stage2/ribench/design-large.csv"], ["Index"])
    require(design.keys() == definition.keys() and len(design) == 1440, "Incomplete design")
    outcomes = indexed((r for r in tables["ribench-outcomes.csv"] if r["method"] in CAPS), ["Index", "method"])
    weights = defaultdict(list)
    for row in tables["selected-component-parameters.csv"]:
        weights[row["dataset_id"], row["method"]].append(row)
    fits = indexed(tables["armB-fits.csv"] + tables["baseline-fits.csv"], ["Index", "method"])
    expected = {(i, method) for (i,) in design for method in CAPS}
    require(fits.keys() == outcomes.keys() == expected, "Incomplete fit/outcome table")
    require(set(weights) == {(f"r-{i}", m) for i, m in expected}, "Incomplete parameter table")
    result = []
    for i, method in sorted(fits, key=lambda key: (int(key[0]), key[1])):
        fit, score = fits[i, method], outcomes[i, method]
        defn, planned = definition[i,], design[i,]
        dataset = fit["dataset_id"]
        context = f"{method}, Index {i}"
        require(dataset == "r-" + i and fit["N"] == planned["N"] == defn["N"] == "50000", context)
        require(fit["Analyte"] == score["Analyte"] == planned["Analyte"] == defn["Analyte"], context)
        require(fit["routine_hash"] == planned["routine_hash"] and fit["anchor_hash"] == planned["anchor_hash"]
                and fit["ribench_hashcode"] == planned["ribench_hashcode"] == defn["HashCode"]
                and fit["n_anchor"] == planned["anchor_n"] == "40"
                and fit["n_routine"] == planned["routine_n"], context)
        k, n = int(fit["components"]), int(planned["routine_n"])
        parts = sorted(weights[dataset, method], key=lambda r: int(r["component"]))
        require(int(fit["max_components"]) == CAPS[method] and 1 <= k <= CAPS[method]
                and score["components"] == fit["components"]
                and [int(r["component"]) for r in parts] == list(range(1, k + 1)), context)
        require(all(r["setting"] == "pooled" and r["transform"] == fit["transform"] == score["transform"]
                    and (r["screened_component"] == "TRUE") == (pos == 0) for pos, r in enumerate(parts)), context)
        shares = [float(r["routine_weight"]) for r in parts]
        require(all(0 <= w <= 1 for w in shares) and abs(sum(shares) - 1) <= 1e-10, context)
        fraction, bg = float(defn["fractionPathol"]), float(defn["bg_fraction"])
        require(fraction == float(fit["pathological_fraction"]) == float(score["fraction"])
                == float(planned["fractionPathol"]) and bg == .001, context)
        healthy = 1 - fraction - bg
        counts = design_counts(defn)
        reconstructed = sum(counts.values())
        realized = counts["phys_n"] / n if reconstructed == n else None
        crp = defn["Analyte"] == "CRP"
        zl = optional_float(score["zzDev_LRL"]) if not crp else None
        zu, d = optional_float(score["zzDev_URL"]), optional_float(score["zzDevAbs_Ov"])
        endpoints = [zu] if crp else [zl, zu]
        scorable = all(v is not None for v in endpoints)
        require(scorable == (d is not None), context)
        if scorable:
            require(math.isclose(d, stats.mean(abs(v) for v in endpoints), rel_tol=1e-12, abs_tol=1e-14), context)
        for side in ("upper",) if crp else ("lower", "upper"):
            raw, stored = optional_float(fit[side]), optional_float(score[{"lower":"LRL","upper":"URL"}[side]])
            require((raw is None and stored is None) or (raw is not None and stored is not None and
                    math.isclose(raw, stored, rel_tol=1e-12, abs_tol=1e-14)), context)
        require(all(fit[flag] == score[flag] for flag in ["declined", "technical"]), context)
        gap = 100 * (shares[0] - healthy)
        result.append(dict(Index=int(i), dataset_id=dataset, method=method, analyte=defn["Analyte"],
                           components=k, transform=fit["transform"], fraction_pathol=fraction,
                           background_fraction=bg, background_inclusive_share=1-fraction,
                           healthy_share_design=healthy, routine_n=n, **counts, reconstructed_n=reconstructed,
                           healthy_share_realized=realized, fitted_share=shares[0], gap_pp=gap,
                           absolute_gap_pp=abs(gap), realized_gap_pp=None if realized is None else 100*(shares[0]-realized),
                           z_lower=zl, z_upper=zu, abs_z_lower=None if zl is None else abs(zl),
                           abs_z_upper=None if zu is None else abs(zu), d=d,
                           extreme=any(abs(v) > 5 for v in endpoints) if scorable else None,
                           report_status=fit["report_status"], routine_hash=fit["routine_hash"], anchor_hash=fit["anchor_hash"]))
    return result


def ranks(values):
    result = [0.] * len(values)
    ordered = sorted(range(len(values)), key=values.__getitem__)
    for _, items in itertools.groupby(enumerate(ordered, start=1), key=lambda entry: values[entry[1]]):
        block = list(items)
        rank = (block[0][0] + block[-1][0]) / 2
        for _, index in block:
            result[index] = rank
    return result


def spearman(x, y):
    if len(set(x)) < 2 or len(set(y)) < 2:
        return None
    return stats.correlation(ranks(x), ranks(y))


def mean(values):
    return stats.fmean(values) if values else None


def quartile(ordered, p):
    """Linear interpolation between adjacent ordered values (R type 7)."""
    if not ordered:
        return None
    pos = (len(ordered) - 1) * p
    low = int(pos)
    return ordered[low] + (pos-low)*(ordered[math.ceil(pos)]-ordered[low])


def summary(rows):
    scored = [r for r in rows if r["d"] is not None]
    errors = [r["d"] for r in scored]
    gaps = sorted(r["gap_pp"] for r in rows)
    revised = [r["realized_gap_pp"] for r in rows if r["realized_gap_pp"] is not None]
    extreme = sum(r["extreme"] for r in scored)
    return dict(n=len(rows), scorable=len(scored), unscorable=len(rows)-len(scored),
                declined=sum(r["report_status"] == "declined" for r in rows),
                mean_fitted_share=mean([r["fitted_share"] for r in rows]),
                mean_true_share=mean([r["healthy_share_design"] for r in rows]), mean_gap_pp=mean(gaps),
                median_gap_pp=quartile(gaps, .5), q25_gap_pp=quartile(gaps,.25), q75_gap_pp=quartile(gaps,.75),
                mean_abs_gap_pp=mean([abs(g) for g in gaps]), mean_d=mean(errors),
                rho_signed=spearman([r["gap_pp"] for r in scored], errors),
                rho_absolute=spearman([r["absolute_gap_pp"] for r in scored], errors),
                extreme_n=extreme, extreme_rate=extreme/len(scored) if scored else None,
                realized_n=len(revised), mean_realized_gap_pp=mean(revised),
                mean_abs_realized_gap_pp=mean([abs(v) for v in revised]))


def write_csv(path, rows, fields=None):
    with path.open("w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=fields or list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task-e-repo", type=Path, default=HERE.parents[2])
    parser.add_argument("--anchorrefine-repo", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=HERE / "archived-shares")
    args = parser.parse_args()
    tables = load_files(args.task_e_repo, TASK_E, TASK_E_DIR, TASK_E_FILES)
    inputs = load_files(args.anchorrefine_repo, STAGE_3, "", INPUT_FILES)
    records = make_records(tables, inputs)
    components, fractions, bands = [], [], []
    fraction_levels = sorted({r["fraction_pathol"] for r in records})
    for method, cap in CAPS.items():
        arm = [r for r in records if r["method"] == method]
        for k in ["all", *range(1, cap + 1)]:
            rows = arm if k == "all" else [r for r in arm if r["components"] == k]
            group = dict(method=method, components=k)
            components.append(dict(**group, **summary(rows)))
            for f in fraction_levels:
                subset = [r for r in rows if r["fraction_pathol"] == f]
                fractions.append(dict(**group, fraction_pathol=f, **summary(subset)))
            for lo, hi, label in BANDS:
                subset = [r for r in rows if lo <= r["absolute_gap_pp"] < hi]
                out = summary(subset)
                bands.append(dict(**group, absolute_gap_band_pp=label,
                                  **{key:out[key] for key in ["n","scorable","unscorable","mean_d","extreme_n","extreme_rate"]}))
    args.out.mkdir(parents=True, exist_ok=True)
    write_csv(args.out/"fits.csv", records)
    write_csv(args.out/"by-components.csv", components)
    write_csv(args.out/"by-components-and-fraction.csv", fractions)
    write_csv(args.out/"gap-bands.csv", bands)
    discrepancy_fields = ["Index", "analyte", "routine_n", "phys_n", "left_n", "right_n", "background_n", "reconstructed_n"]
    discrepancies = [{k:r[k] for k in discrepancy_fields} for r in records
                     if r["method"] == "armB" and r["routine_n"] != r["reconstructed_n"]]
    write_csv(args.out/"count-discrepancies.csv", discrepancies, discrepancy_fields)
    sources = [
        dict(repository="amromeo/amromeo.github.io", revision=TASK_E,
             files={TASK_E_DIR + n:h for n,h in TASK_E_FILES.items()}),
        dict(repository="amromeo/anchorRefine", revision=STAGE_3, files=INPUT_FILES),
    ]
    (args.out/"sources.json").write_text(json.dumps(sources, indent=2)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
