#!/usr/bin/env python3
"""Summarize selected optimizer exits in saved, scored Task H first-stage fits."""
import argparse
import csv
import gzip
import math
from collections import Counter, defaultdict
from pathlib import Path

# anchoredRI/R/mixture_fit.R::.ri_optimizer_status labels the selected candidate.
# An unavailable selection has no optimizer status; a safeguard decline can have one.
STATUSES = ("converged", "iteration_limit", "line_search_failure", "other_exit", "unavailable")
STAGES = ("before_check", "first_stage_released")
METRICS = ("z_lower", "abs_z_lower", "z_upper", "abs_z_upper", "d", "any_extreme")
ARMS = [(fit_type, cap, bank) for fit_type in ("full", "W", "subsample")
        for cap in (3, 5) for bank in (("A", "B") if fit_type == "subsample" else ("",))]


class Scores:
    def __init__(self):
        self.records = 0
        self.n = [0] * len(METRICS)
        self.sums = [0.0] * len(METRICS)

    def add(self, values):
        self.records += 1
        for i, value in enumerate(values):
            if value is not None:
                self.n[i] += 1
                self.sums[i] += value

    def means(self):
        return [s / n if n else None for s, n in zip(self.sums, self.n)]


class Group(Scores):
    def __init__(self):
        super().__init__()
        self.indices = set()
        self.status_counts = Counter()

    def record(self, index, first, scores):
        self.indices.add(index)
        self.status_counts[first] += 1
        self.add(scores)

    def size(self):
        return sum(self.status_counts.values())


def number(x):
    if x in ("", "NA", "NaN"):
        return None
    value = float(x)
    return value if math.isfinite(value) else None


def read_scores(row):
    lower, upper, d = (number(row[field]) for field in ("z_lower", "z_upper", "d"))
    if row["Analyte"] == "CRP" and lower is not None:
        raise ValueError(f"One-sided CRP acquired a lower score at Index {row['Index']}")
    requested = [upper] if row["Analyte"] == "CRP" else [lower, upper]
    # Count a known extreme even with a missing companion score. A non-extreme
    # call requires every requested score. Keep that rate's own denominator.
    extreme = None
    if any(x is not None and abs(x) > 5 for x in requested):
        extreme = 1
    elif all(x is not None for x in requested):
        extreme = 0
    return [lower, abs(lower) if lower is not None else None,
            upper, abs(upper) if upper is not None else None, d, extreme]


def subtract(left, right):
    return [a - b if a is not None and b is not None else None for a, b in zip(left, right)]


def fields():
    return ["fit_type", "cap", "seed_set", "selected_k", "phase", "optimizer_status"]


def write_table(path, rows):
    with path.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def analyze(scores_path, out):
    groups = defaultdict(Group)
    within_index = defaultdict(Scores)
    seen = set()
    per_index = Counter()
    opener = gzip.open if scores_path.suffix == ".gz" else open
    with opener(scores_path, "rt", newline="") as source:
        for row in csv.DictReader(source):
            index = int(row["Index"])
            arm = row["fit_type"], int(row["cap"]), row["seed_set"]
            kind, cap, bank = arm
            draw = row["draw"]
            if arm not in ARMS or (kind == "subsample" and int(draw) not in range(1, 51)) or (
                    kind != "subsample" and draw):
                raise ValueError(f"Unexpected fit assignment at Index {index}: {arm}, draw {draw}")
            key = (index, *arm, int(draw) if draw else None)
            if key in seen:
                raise ValueError(f"Duplicate fit: {key}")
            seen.add(key)
            per_index[index] += 1
            raw_k = number(row["components"])
            selected = int(raw_k) if raw_k is not None else "unavailable"
            if raw_k is not None and (raw_k != selected or selected not in range(1, cap + 1)):
                raise ValueError(f"Invalid selected component count at {key}: {raw_k}")
            status = row["optimizer_status"]
            if status in ("", "NA"):
                status = "unavailable"
            first = row["first_status"]
            if status not in STATUSES or first not in ("reported", "declined", "withheld"):
                raise ValueError(f"Unrecognized optimizer or first-stage status at {key}: {status}, {first}")
            values = read_scores(row)
            for phase in STAGES if first == "reported" else STAGES[:1]:
                for selected_k in ("all", selected):
                    group = (*arm, selected_k, phase)
                    groups[(*group, status)].record(index, first, values)
                    if kind == "subsample":
                        within_index[(*group, index, status)].add(values)
    if not per_index or any(n != 204 for n in per_index.values()):
        raise ValueError("Expected 204 unique first-stage fit records for every Index.")
    summary = []
    paired = []
    for arm in ARMS:
        kind, cap, bank = arm
        for k in ("all", *range(1, cap + 1), "unavailable"):
            for phase in STAGES:
                group = (*arm, k, phase)
                total = sum(groups[(*group, s)].size() for s in STATUSES)
                before_count = sum(groups[(*arm, k, STAGES[0], s)].size() for s in STATUSES)
                converged = groups[(*group, "converged")]
                for status in STATUSES:
                    current = groups[(*group, status)]
                    base = dict(zip(fields(), (*group, status)))
                    record = dict(base, all_attempts=before_count, phase_fits=total,
                                  fits=current.size(), fraction_of_phase=current.size() / total if total else None,
                                  distinct_datasets=len(current.indices), reported=current.status_counts["reported"],
                                  declined=current.status_counts["declined"], withheld=current.status_counts["withheld"],
                                  extreme_fits=round(current.sums[-1]), converged_fits=converged.size())
                    means = current.means()
                    diffs = subtract(means, converged.means())
                    for i, name in enumerate(METRICS):
                        record[f"{name}_n"] = current.n[i]
                        record[f"{name}_mean"] = means[i]
                        record[f"{name}_converged_n"] = converged.n[i]
                        record[f"{name}_minus_converged"] = diffs[i]
                    summary.append(record)
                    if kind != "subsample" or status == "converged":
                        continue
                    contrasts = Scores()
                    jointly_observed = current.indices & converged.indices
                    for index in sorted(jointly_observed):
                        lhs = within_index[(*group, index, status)]
                        rhs = within_index[(*group, index, "converged")]
                        contrasts.add(subtract(lhs.means(), rhs.means()))
                    pairs = dict(base, datasets_with_both_exits=len(jointly_observed),
                                 exit_fits_in_pairs=sum(
                                     within_index[(*group, i, status)].records for i in jointly_observed),
                                 converged_fits_in_pairs=sum(
                                     within_index[(*group, i, "converged")].records for i in jointly_observed))
                    for i, (name, diff) in enumerate(zip(METRICS, contrasts.means())):
                        pairs[f"{name}_paired_datasets"] = contrasts.n[i]
                        pairs[f"{name}_exit_minus_converged"] = diff
                    paired.append(pairs)
    out.mkdir(exist_ok=True, parents=True)
    write_table(out / "optimizer-summary.csv", summary)
    write_table(out / "optimizer-within-index.csv", paired)
    write_report(out / "OPTIMIZER.md", summary, paired, len(per_index))


def display(x):
    return "—" if x is None else f"{x:.4f}"


def arm_label(row):
    kind = {"full": "Full", "W": "W", "subsample": "Subsamples"}[row["fit_type"]]
    return kind + (" " + row["seed_set"] if row["seed_set"] else "")


def write_report(path, summary, paired, n):
    text = [
        "# Selected optimizer exits and endpoint error",
        "",
        f"The analysis covers {n:,} datasets and {204 * n:,} first-stage fit records. "
        "Each subsample contributes a record under each cap. The full and W fits enter once per cap. "
        "The denominator for an exit rate includes every attempted fit in that group. A failed selection appears as unavailable.",
        "",
        "The frozen fitter labels code 0 as converged. It labels code 1 or an iteration-limit message as iteration_limit, "
        "a line-search or LNSRCH message as line_search_failure, and other codes as other_exit. "
        "The label belongs to the selected mixture candidate. Either check can decline a fit after any optimizer exit. "
        "[Status mapping in anchoredRI 0.4.0](https://github.com/amromeo/amromeo.github.io/blob/fe8105ac6961f1b45bbfbe034648938f227c4726/anchoredRI/R/mixture_fit.R)",
        "",
        "## Before the safeguards",
        "",
        "| Fit | Cap | Attempts | Converged | Iteration limit | Line-search failure | Other exit | Unavailable |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    raw = [r for r in summary if r["selected_k"] == "all" and r["phase"] == STAGES[0]]
    for arm in ARMS:
        rows = [r for r in raw if tuple(r[k] for k in ("fit_type", "cap", "seed_set")) == arm]
        text.append(f"| {arm_label(rows[0])} | {arm[1]} | {rows[0]['phase_fits']:,} | " +
                    " | ".join(f"{r['fits']:,}" for r in rows) + " |")
    text.extend(["", "Each cell below gives iteration-limit exits / line-search exits / fits selecting that K. "
                 "The values use all attempted fits before the safeguards.", "",
                 "| Fit | Cap | K = 1 | K = 2 | K = 3 | K = 4 | K = 5 |",
                 "|---|---:|---:|---:|---:|---:|---:|"])
    for arm in ARMS:
        arm_rows = [r for r in summary if tuple(r[k] for k in ("fit_type", "cap", "seed_set")) == arm
                    and r["phase"] == STAGES[0]]
        cells = []
        for k in range(1, 6):
            if k > arm[1]:
                cells.append("—")
                continue
            exits = {r["optimizer_status"]: r for r in arm_rows if r["selected_k"] == k}
            cells.append(f"{exits['iteration_limit']['fits']:,} / {exits['line_search_failure']['fits']:,} / "
                         f"{exits['converged']['phase_fits']:,}")
        text.append(f"| {arm_label(arm_rows[0])} | {arm[1]} | " + " | ".join(cells) + " |")
    text.extend([
        "",
        "Signed lower and upper z scores measure direction. Absolute z measures endpoint error. "
        "The dataset score d is the uncapped mean absolute z across its requested endpoints. CRP uses only its upper endpoint. "
        "An extreme fit has at least one requested endpoint with absolute z above 5. Each metric has its own denominator.",
        "",
        "| Fit | Cap | Exit | Scored fits | Mean d | Converged mean d | Difference | Fits with an extreme endpoint / assessed |",
        "|---|---:|---|---:|---:|---:|---:|---:|",
    ])
    for row in raw:
        if row["optimizer_status"] == "converged" or not row["fits"]:
            continue
        ref = next(other for other in raw if other["fit_type"] == row["fit_type"] and
                   other["cap"] == row["cap"] and other["seed_set"] == row["seed_set"] and
                   other["optimizer_status"] == "converged")
        text.append(f"| {arm_label(row)} | {row['cap']} | {row['optimizer_status']} | {row['d_n']:,} | "
                    f"{display(row['d_mean'])} | {display(ref['d_mean'])} | "
                    f"{display(row['d_minus_converged'])} | {row['extreme_fits']:,} / {row['any_extreme_n']:,} |")
    text.extend([
        "",
        "The detailed table gives counts, denominators, exit rates, signed and absolute endpoint errors, and differences "
        "from converged fits for each cap and selected component count. The released analysis uses the same first-stage "
        "scores on fits passing both applicable safeguards. [All exit and error summaries](optimizer-summary.csv)",
        "",
        "## Subsamples compared within the same dataset",
        "",
        "For each dataset with both exit types, we average each metric over its draws separately by exit and subtract the "
        "converged mean from the other exit mean. Each dataset with both means contributes once. The selected-K rows "
        "also hold the selected component count fixed. Full and W have one fit per cap per dataset.",
        "",
        "| Seed set | Cap | Exit | Datasets with paired d | Mean within-dataset difference in d |",
        "|---|---:|---|---:|---:|",
    ])
    for row in paired:
        if row["selected_k"] == "all" and row["phase"] == STAGES[0] and row["datasets_with_both_exits"]:
            text.append(f"| {row['seed_set']} | {row['cap']} | {row['optimizer_status']} | "
                        f"{row['d_paired_datasets']:,} | {display(row['d_exit_minus_converged'])} |")
    text.extend(["", "These comparisons describe the saved fits. Optimizer exit and model selection both depend "
                 "on the draw. [Within-dataset contrasts before and after the checks](optimizer-within-index.csv)", ""])
    path.write_text("\n".join(text))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("scores", type=Path, help="fit-scores.csv.gz from the Task H scorer")
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    analyze(args.scores, args.output)
