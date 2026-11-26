# SPDX-License-Identifier: GPL-3.0-only
"""Localize the signed healthy-share deficit in the archived Task H full fits.

Run with Python 3.10+ and the verified per-Index release files. See REPORT.md.
"""
import argparse
import csv
from collections import defaultdict
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR
from hashlib import sha256
import io
import math
from pathlib import Path
from statistics import fmean

HERE = Path(__file__).resolve().parent
STUDY = HERE.parent.parent
FULL = STUDY / "results/full-45860591866"
DESIGN_SHA = "177cd392896cb300bf4aaa35a959e46bb3ee48d510beff8cf9957f6643a66e39"
MANIFEST_SHA = "f09924475a9b4410962d0f7c64b2c54705902e516feb7f3f2dfbb48bc022ad8d"


def check(ok, detail):
    if not ok:
        raise ValueError(detail)


def read_verified(path, expected):
    data = path.read_bytes()
    check(sha256(data).hexdigest() == expected, f"SHA-256 mismatch at {path}")
    return csv.DictReader(io.StringIO(data.decode()))


def phi(z):
    return .5 * math.erfc(-z / math.sqrt(2))


def healthy_cdf(x, spec):
    lam, mu, sd, shift = (float(spec["nonp_" + k]) for k in ("lambda", "mu", "sigma", "shift"))
    check(lam >= 0 and shift >= 0 and sd > 0, "Unexpected RIbench distribution")
    if x <= 0:
        return 0.0
    if lam == 1:
        # invBoxCox(z,1) + shift = z + 1 + shift. The generator removes x < 0.
        base = -1 - shift
        latent = x - 1 - shift
    else:
        if x <= shift:
            return 0.0
        base = -1 / lam if lam else -math.inf
        latent = math.expm1(lam * math.log(x - shift)) / lam if lam else math.log(x - shift)
    excluded = phi((base - mu) / sd)
    return (phi((latent - mu) / sd) - excluded) / (1 - excluded)


def component_cdf(x, fit):
    family = fit["transform"]
    mu, sd, shift = (float(fit[key]) for key in ("mu", "sigma", "shift"))
    check(sd > 0, "Invalid component spread")
    if family == "identity":
        latent = x
    else:
        check(family in ("log", "shifted_log") and (family != "log" or shift == 0),
              f"Unexpected fitted transform {family} with shift {shift}")
        if x <= shift:
            return 0.0
        latent = math.log(x - shift)
    return phi((latent - mu) / sd)


def region_limits(spec, rounded):
    one_sided = spec["Analyte"] == "CRP"
    lo = None if one_sided else Decimal(spec["GT_LRL"])
    hi = Decimal(spec["GT_URL"])
    if rounded:
        unit = Decimal(1).scaleb(-int(spec["decimals"]))
        if lo is not None:
            lo = ((lo / unit).to_integral_value(rounding=ROUND_CEILING) - Decimal(".5")) * unit
        hi = ((hi / unit).to_integral_value(rounding=ROUND_FLOOR) + Decimal(".5")) * unit
    return float(lo) if lo is not None else -math.inf, float(hi)


def regional_record(spec, fit, h, w, rounded):
    lo, hi = region_limits(spec, rounded)
    out = dict(Index=int(spec["Index"]), analyte=spec["Analyte"], fraction_pathol=float(spec["fractionPathol"]),
               scale="rounded observations" if rounded else "continuous", first_stage=fit["first_status"],
               transform=fit["transform"], mu=fit["mu"], sigma=fit["sigma"], shift=fit["shift"],
               healthy_share=h, component1_share=w, lower_threshold=lo, upper_threshold=hi)
    for group, weight, cdf in (("healthy", h, lambda x: healthy_cdf(x, spec)),
                               ("component1", w, lambda x: component_cdf(x, fit))):
        left, through_upper = cdf(lo), cdf(hi)
        check(0 <= left <= through_upper <= 1, f"Invalid CDF at {spec['Index']}")
        for region, prob in (("below", left), ("inside", through_upper-left), ("above", 1-through_upper),
                             ("outside", left+1-through_upper), ("total", 1)):
            out[f"{group}_{region}_pp"] = 100 * weight * prob
    for region in ("below", "inside", "above", "outside", "total"):
        out[f"deficit_{region}_pp"] = out[f"healthy_{region}_pp"] - out[f"component1_{region}_pp"]
    check(abs(out["deficit_inside_pp"] + out["deficit_outside_pp"] - 100 * (h-w)) < 2e-12,
          f"Mass failed to reconcile at {spec['Index']}")
    return out


def summarize(key, rows):
    scale, group, value = key
    out = dict(scale=scale, group=group, value=value, fits=len(rows))
    for region in ("total", "inside", "outside", "below", "above"):
        for measure in ("healthy", "component1", "deficit"):
            field = f"{measure}_{region}_pp"
            out[field] = fmean(r[field] for r in rows)
        deficit = out[f"deficit_{region}_pp"]
        healthy = out[f"healthy_{region}_pp"]
        out[f"pct_of_total_deficit_{region}"] = 100*deficit/out["deficit_total_pp"] if out["deficit_total_pp"] else None
        out[f"pct_of_true_mass_missing_{region}"] = 100*deficit/healthy if healthy else None
    return out


def write(path, rows):
    with path.open("w", newline="") as target:
        writer = csv.DictWriter(target, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("records", type=Path, help="verified release tree, <Index>/fits.csv")
    parser.add_argument("testsets", type=Path, help="pinned anchorRefine stage1/large/testsets.csv")
    parser.add_argument("--out", type=Path, default=HERE)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    manifest = {r["Index"]: r["sha256"] for r in read_verified(FULL / "verification/file-hashes.csv", MANIFEST_SHA)
                if r["file"] == "fits.csv"}
    specs = list(read_verified(args.testsets, DESIGN_SHA))
    check(len(specs) == 1440 and {s["Index"] for s in specs} == set(manifest), "Incorrect cohort")
    rows = []
    for spec in sorted(specs, key=lambda s: int(s["Index"])):
        idx = spec["Index"]
        fit, = (r for r in read_verified(args.records / idx / "fits.csv", manifest[idx])
                if r["cap"] == "5" and r["fit_type"] == "full")
        check(fit["Index"] == idx and spec["N"] == "50000", f"Index mismatch {idx}")
        if fit["components"] != "5":
            continue
        h = 1 - float(spec["fractionPathol"]) - float(spec["bg_fraction"])
        w = float(fit["share"])
        check(0 < w <= 1 and 0 < h <= 1, f"Invalid routine weight at {idx}")
        # Confirm the stored parameters are for the unmasked first-stage model.
        one_sided = spec["Analyte"] == "CRP"
        check(abs(component_cdf(float(fit["raw_upper"]), fit) - (.95 if one_sided else .975)) < 1e-9,
              f"Upper quantile mismatch at {idx}")
        if not one_sided:
            check(abs(component_cdf(float(fit["raw_lower"]), fit) - .025) < 1e-9,
                  f"Lower quantile mismatch at {idx}")
        rows.extend(regional_record(spec, fit, h, w, rounded) for rounded in (False, True))
    check(len(rows) == 2 * 934, "Incorrect five-component cohort")
    groups = defaultdict(list)
    for r in rows:
        for group, key in (("all", "all"), ("analyte", r["analyte"]), ("first stage", r["first_stage"]),
                           ("net deficit", "positive" if r["deficit_total_pp"] > 0 else "zero or negative")):
            groups[r["scale"], group, key].append(r)
    summaries = [summarize(key, members) for key, members in sorted(groups.items())]
    write(args.out / "per-dataset.csv", rows)
    write(args.out / "summary.csv", summaries)


if __name__ == "__main__":
    main()
