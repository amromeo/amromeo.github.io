# SPDX-License-Identifier: GPL-3.0-only
"""Reconstruct the four-step policy from the committed Task H paired records."""

import csv
import hashlib
import math
from pathlib import Path
from statistics import mean, median

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parents[1] / "results/full-45860591866/paired-index.csv"
SOURCE_SHA256 = "d7a3fe0958cdc73deabbb578483596655c156b1b04d2785a5763139f97570d3c"
FIELDS = ("lower", "upper", "z_lower", "z_upper", "d")


def interval(row, prefix):
    values = {}
    for field in FIELDS:
        value = float(row[f"{prefix}_{field}"])
        values[field] = value if math.isfinite(value) else None
    return values


def write_table(filename, rows):
    with (HERE / filename).open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def apply_policy(row):
    if row["panel"] not in ("A", "B"):
        raise ValueError(f"Unknown first-stage panel at {row['Index']}")
    if any(row[key] not in ("True", "False") for key in ("full_released", "median_released")):
        raise ValueError(f"Unknown release decision at {row['Index']}")
    full_passes = row["panel"] == "A"
    package_returns = row["full_released"] == "True"
    majority_passes = row["median_released"] == "True"
    package = interval(row, "ref")
    subsample = interval(row, "median")
    ends = ("upper",) if row["analyte"] == "CRP" else ("lower", "upper")

    if full_passes and not package_returns:
        raise ValueError(f"Full fit passed without a package return at {row['Index']}")
    for values, expected in [(package, package_returns), (subsample, True)]:
        checked = ["d", *ends, *("z_" + e for e in ends)]
        if any((values[key] is not None) != expected for key in checked):
            raise ValueError(f"Inconsistent saved limits or score at {row['Index']}")
        if expected and not math.isclose(values["d"], mean(abs(values["z_" + e]) for e in ends),
                                         rel_tol=1e-13, abs_tol=1e-15):
            raise ValueError(f"Saved score disagrees with endpoints at {row['Index']}")

    if full_passes:
        step, result = "a", subsample
    elif package_returns:
        step, result = "b", package
    elif majority_passes:
        step, result = "c", subsample
    else:
        step, result = "d", package

    output = dict(bank=row["bank"], Index=int(row["Index"]), analyte=row["analyte"],
                  full_first_stage_passes=full_passes, package_returns=package_returns,
                  majority_approves=majority_passes, step=step)
    for label, values in [("package", package), ("policy", result)]:
        output.update({f"{label}_{key}": value for key, value in values.items()})
        released = values["d"] is not None
        output[f"{label}_score_all"] = min(5, values["d"]) if released else 5
        output[f"{label}_extreme"] = any(abs(values["z_" + e]) > 5 for e in ends) if released else None
    output["paired_delta_d"] = result["d"] - package["d"] if package_returns else None
    return output


def summarize(records, label):
    returned = [row for row in records if row[f"{label}_d"] is not None]
    return dict(assigned=len(records), returns=len(returned),
                nonreturns=len(records) - len(returned),
                released_mean_d=mean(row[f"{label}_d"] for row in returned) if returned else None,
                released_max_d=max((row[f"{label}_d"] for row in returned), default=None),
                all_assignment_mean=mean(row[f"{label}_score_all"] for row in records),
                released_intervals_with_extreme_endpoint=sum(row[f"{label}_extreme"] for row in returned))


def main():
    with SOURCE.open("rb") as handle:
        if hashlib.file_digest(handle, "sha256").hexdigest() != SOURCE_SHA256:
            raise ValueError("The input hash differs from the published Task H paired records")
    with SOURCE.open(newline="") as handle:
        records = [apply_policy(row) for row in csv.DictReader(handle)
                   if row["cap"] == "5" and row["definition"] == "all"]
    records.sort(key=lambda r: (r["bank"], r["Index"]))
    by_bank = {bank: [r for r in records if r["bank"] == bank] for bank in ("A", "B")}
    indices = [r["Index"] for r in by_bank["A"]]
    if (len(records) != 2880 or len(set(indices)) != 1440 or
            indices != [r["Index"] for r in by_bank["B"]]):
        raise ValueError("Both seed sets must contain the same 1,440 distinct datasets")
    for a, b in zip(by_bank["A"], by_bank["B"]):
        shared = ["analyte", "full_first_stage_passes", "package_returns", *("package_" + x for x in FIELDS)]
        if any(a[k] != b[k] for k in shared):
            raise ValueError(f"Full fit or package changed across seeds at {a['Index']}")

    summary = [dict(bank="", method="package", **summarize(by_bank["A"], "package"))]
    routes, paired = [], []
    for bank, group in by_bank.items():
        summary.append(dict(bank=bank, method="policy", **summarize(group, "policy")))
        for step in "abcd":
            cases = [r for r in group if r["step"] == step]
            routes.append(dict(bank=bank, step=step, **summarize(cases, "policy")))
        pairs = [r for r in group if r["paired_delta_d"] is not None]
        delta = [r["paired_delta_d"] for r in pairs]
        paired.append(dict(bank=bank, common_returns=len(pairs),
                           policy_mean_d=mean(r["policy_d"] for r in pairs),
                           package_mean_d=mean(r["package_d"] for r in pairs),
                           mean_policy_minus_package=mean(delta), median_policy_minus_package=median(delta),
                           policy_wins=sum(x < 0 for x in delta), exact_ties=sum(x == 0 for x in delta),
                           policy_losses=sum(x > 0 for x in delta)))
    write_table("by-index.csv", records)
    write_table("summary.csv", summary)
    write_table("routes.csv", routes)
    write_table("common-returns.csv", paired)


if __name__ == "__main__":
    main()
