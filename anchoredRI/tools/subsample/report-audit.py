"""Aggregate cap-3 archive checks. Input contains the 48 shards, possibly split into parts."""
import argparse
import csv
import gzip
import hashlib
import json
import math
import statistics
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXACT = ("transform", "components", "lrt_flag", "lower_flag", "status")
ENDPOINTS = ("raw_lower", "raw_upper", "lower", "upper")
CHECKS = (*EXACT, *(f"{k}_missingness" for k in ENDPOINTS), *ENDPOINTS, "p")
ARCHIVE_URL = ("https://github.com/amromeo/amromeo.github.io/blob/"
               "8bcd49c3c7977c50982524b9dbc48d16bdfa5914/"
               "anchoredRI/tools/boxcox/accuracy/large/results/baseline-fits.csv")


def read(path):
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def write(path, rows, columns):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def number(value):
    return float(value) if value else None


def evaluate(row):
    """Use saved 17-digit values to independently reconstruct the R worker's assertions."""
    results = {}
    crp = row["Analyte"] == "CRP"
    width = float(row["archive_width"])
    if not math.isfinite(width) or width <= 0:
        raise ValueError(f"Invalid archived span at {row['Index']}")
    for key in EXACT:
        old, new = row[f"archive_{key}"], row[f"observed_{key}"]
        if key == "lower_flag" and crp:
            status = "inapplicable"
        elif not old or not new:
            status = "unavailable"
        else:
            status = "pass" if new == old else "fail"
        results[key] = (status, None)
    for key in ENDPOINTS:
        old, new = (number(row[f"{source}_{key}"]) for source in ("archive", "observed"))
        missing_ok = (old is None) == (new is None)
        results[f"{key}_missingness"] = ("pass" if missing_ok else "fail", None)
        if crp and "lower" in key:
            if old is not None:
                raise ValueError(f"Archived CRP lower limit on {row['Index']}")
            results[key] = ("inapplicable", None)
        elif old is None and new is None and key in ("lower", "upper"):
            results[key] = ("inapplicable", None)
        elif old is None or new is None:
            results[key] = ("unavailable", None)
        else:
            delta = abs(new - old)
            results[key] = ("pass" if delta <= .001 * width else "fail", delta)
    p_old, p_new = (number(row[f"{source}_p"]) for source in ("archive", "observed"))
    if p_old is None or p_new is None:
        results["p"] = ("unavailable", None)
    else:
        delta = abs(p_new - p_old)
        results["p"] = ("pass" if delta <= .0001 else "fail", delta)
    passed = not any(status in ("unavailable", "fail") for status, _ in results.values())
    if passed != (row["pass"] == "TRUE"):
        raise ValueError(f"Independent decision differs at {row['Index']} ({row['attempt']})")
    return results


def check_archive_identity(rows):
    # The archived-shares export independently pairs these 1,440 IDs and hashes.
    references = {r["Index"]: r for r in read(HERE / "archived-shares" / "fits.csv")
                  if r["method"] == "anchoredRI"}
    if len(references) != 1440 or {r["Index"] for r in rows} != set(references):
        raise ValueError("The audit must cover all 1,440 Task E baseline datasets")
    for r in rows:
        ref = references[r["Index"]]
        for field in ("routine_hash", "anchor_hash"):
            if r[field] != ref[field]:
                raise ValueError(f"Archived {field} mismatch at {r['Index']}")
        if r["archive_seed"] != r["Index"] or int(r["anchor_seed"]) != 17000000 + int(r["Index"]):
            raise ValueError(f"Invalid seeds at {r['Index']}")
        if abs(float(r["archive_share"]) - float(ref["fitted_share"])) > 1e-14:
            raise ValueError(f"Component-1 archive share mismatch at {r['Index']}")
        if r["observed_share"] and not math.isclose(float(r["weight_sum"]), 1, rel_tol=0, abs_tol=1e-10):
            raise ValueError(f"Invalid recorded routine shares at {r['Index']}")


def summarize(rows, checks):
    result = []
    for field in CHECKS:
        records = [(r, checks[(r["Index"], r["attempt"])][field]) for r in rows]
        counts = Counter(status for _, (status, _) in records)
        differences = [(d, d / float(r["archive_width"]), r["Index"])
                       for r, (_, d) in records if d is not None]
        worst = max(differences, default=(None, None, None), key=lambda x: x[0])
        width_max = max((v[1] for v in differences), default=None) if field in ENDPOINTS else None
        result.append(dict(check=field, comparable=counts["pass"] + counts["fail"],
                           passed=counts["pass"], failed=counts["fail"],
                           unavailable=counts["unavailable"], inapplicable=counts["inapplicable"],
                           max_abs=worst[0], max_abs_Index=worst[2], max_fraction_width=width_max))
    return result


def show(value):
    return "—" if value is None else f"{value:.7g}"


def table(summary):
    lines = ["| Check | Comparable | Pass | Fail | Unavailable | Inapplicable | Max abs. difference | Max fraction of width |",
             "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for item in summary:
        lines.append("| {check} | {comparable} | {passed} | {failed} | {unavailable} | {inapplicable} | ".format(**item)
                     + f"{show(item['max_abs'])} | {show(item['max_fraction_width'])} |")
    return lines


def load(root):
    rows = []
    for path in sorted(root.rglob("cap3-shard-*-part-*.csv")):
        current = read(path)
        fit_file = path.with_suffix(".rds")
        if not current or not fit_file.is_file() or fit_file.stat().st_size == 0:
            raise ValueError(f"Missing fit evidence or empty shard at {path}")
        rows.extend(current)
    if not rows:
        raise ValueError(f"No cap-3 audit outputs in {root}")
    originals = sorted((r for r in rows if r["attempt"] == "study"), key=lambda r: int(r["Index"]))
    repeats = sorted((r for r in rows if r["attempt"] == "archived_seed"), key=lambda r: int(r["Index"]))
    checks = {(r["Index"], r["attempt"]): evaluate(r) for r in rows}
    if len(checks) != len(rows) or len(originals) + len(repeats) != len(rows):
        raise ValueError("Unexpected or duplicate diagnostic attempt")
    check_archive_identity(originals)
    for r in originals:
        if int(r["fit_seed"]) != 1:
            raise ValueError(f"Wrong primary seed at {r['Index']}")
    failed = {r["Index"]: r for r in originals if r["pass"] != "TRUE"}
    if set(failed) != {r["Index"] for r in repeats}:
        raise ValueError("Missing or spurious archived-seed reruns")
    for r in repeats:
        primary = failed[r["Index"]]
        if r["fit_seed"] != r["Index"] or any(r[k] != primary[k] for k in r
                if k.startswith("archive_") or k in ("routine_hash", "anchor_hash", "Analyte")):
            raise ValueError(f"Changed reference or incorrect archived seed at {r['Index']}")
    return originals, repeats, checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workers", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    original, repeats, checks = load(args.workers)
    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    exceptions = [r for r in original if r["pass"] != "TRUE"]
    repeat_by_id = {r["Index"]: r for r in repeats}
    columns = list(original[0])
    write(output / "archive-audit.csv.gz", original, columns)
    write(output / "exceptions.csv", exceptions, columns)
    write(output / "seed-repeats.csv", repeats, columns)
    report = ["# Full-data cap-3 agreement with the archived Task E baseline", "",
              f"The study-seed fit passed every required check on **{len(original) - len(exceptions):,} "
              f"of {len(original):,}** matched datasets. The archived seed resolved "
              f"**{sum(r['pass'] == 'TRUE' for r in repeats)} of {len(repeats)}** original exceptions.", "",
              "The study seed was 1. Each diagnostic rerun used the original RIbench Index as its fitting seed. "
              "Both runs used the same archived routine values, anchors, rounding and requested percentiles. "
              f"[Task E baseline]({ARCHIVE_URL})", ""]
    diagnostics = {}
    summary_csv = []
    for label, records, title in (("study", original, "Study seed"), ("archived_seed", repeats, "Archived-seed diagnostic")):
        summary = summarize(records, checks)
        summary_csv.extend(dict(attempt=label, **entry) for entry in summary)
        share_deltas = [float(r["share_difference_pp"]) for r in records if r["share_difference_pp"]]
        data = dict(attempts=len(records), passed=sum(r["pass"] == "TRUE" for r in records),
                    cpu_s=sum(float(r["cpu_s"]) for r in records), fit_elapsed_s=sum(float(r["elapsed_s"]) for r in records),
                    comparable_shares=len(share_deltas), mean_share_diff_pp=statistics.mean(share_deltas) if share_deltas else None,
                    max_abs_share_diff_pp=max((abs(x) for x in share_deltas), default=None))
        diagnostics[label] = data
        report.extend([f"## {title}", ""])
        if records:
            report.extend(table(summary))
            report.extend(["", "Endpoint differences are in analyte units. For two-sided fits the denominator is the "
                           "archived raw interval width. For CRP it is the archived upper limit minus zero. The endpoint "
                           "tolerance is 0.001 of that span. The absolute LRT p-value tolerance is 0.0001. Lower-endpoint "
                           "missingness is also checked for CRP. Two matching declines have no released endpoints to compare.",
                           "", f"Component 1's mean study-minus-archive share difference was "
                           f"{show(data['mean_share_diff_pp'])} percentage points across {len(share_deltas):,} matched "
                           f"shares. The maximum absolute difference was {show(data['max_abs_share_diff_pp'])} percentage "
                           "points. Share differences are descriptive.", ""])
        else:
            report.extend(["No original exception required a rerun.", ""])
    report.extend(["## Exceptions by original Index", ""])
    if exceptions:
        report.extend(["| Index | Analyte | Study-seed differences | Archived-seed differences |",
                       "|---:|---|---|---|"])
        for row in exceptions:
            second = repeat_by_id[row["Index"]]
            report.append(f"| {row['Index']} | {row['Analyte']} | {row['differences'].replace('|', ', ')} | "
                          f"{second['differences'].replace('|', ', ') or 'resolved'} |")
    else:
        report.append("No study-seed comparison exceeded a tolerance or differed on a required exact check.")
    report.extend(["", "The files retain the Task E and study values on the same row, the two input hashes, "
                   "fit seed, first-stage flags and decisions, component-1 shares and fit times. "
                   "[All original comparisons](archive-audit.csv.gz) · [Original exceptions](exceptions.csv) · "
                   "[Archived-seed reruns](seed-repeats.csv) · [Check summaries](audit-checks.csv)", ""])
    (output / "AUDIT.md").write_text("\n".join(report), encoding="utf-8")
    write(output / "audit-checks.csv", summary_csv, list(summary_csv[0]))
    (output / "audit.json").write_text(json.dumps(diagnostics, indent=2, allow_nan=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
