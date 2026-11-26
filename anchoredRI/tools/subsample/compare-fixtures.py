#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-only
"""Verify exact repeat and separate-host Task H fixture evidence."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def read_csv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for label in ("first", "same", "other", "out"):
        parser.add_argument(label, type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    if any(args.out.iterdir()):
        parser.error("Use an empty audit output directory to retain prior attempts.")

    records = []

    def record(kind, case, field, a, b, ok, difference="", limit=""):
        records.append(dict(comparison=kind, case=case, field=field, reference=a,
                            observed=b, difference=difference, tolerance=limit, passed=bool(ok)))

    runs = (args.first, args.same, args.other)
    hosts = []
    for label, run in zip(("first", "same", "other"), runs):
        rows = read_csv(run / "checks.csv")
        record("fixtures", label, "unique checks", len(rows), len({r["test"] for r in rows}),
               bool(rows) and len(rows) == len({r["test"] for r in rows}))
        for row in rows:
            record("fixtures", label, row["test"], "TRUE", row["pass"], row["pass"] == "TRUE",
                   row["detail"])
        host, = read_csv(run / "host.csv")
        hosts.append(host)
    for kind, host in zip(("same-host", "cross-host"), hosts[1:]):
        boot1 = hosts[0]["boot_id"]
        boot2 = host["boot_id"]
        record(kind, "host", "boot_id", boot1, boot2,
               bool(boot1 and boot2) and ((boot1 == boot2) == (kind == "same-host")))
        for field in ("R", "package"):
            record(kind, "runtime", field, hosts[0][field], host[field], hosts[0][field] == host[field])
    unique_runs = len({h["invocation"] for h in hosts if h["invocation"]})
    record("fixtures", "independent invocations", "invocation UUID count", 3, unique_runs, unique_runs == 3)
    # science.rds retains double-precision selected parameters, candidate
    # scores, safeguards and both reporting decisions; CSVs are human-readable.
    for file in ("science.rds", "fits.csv", "positions.csv", "checks.csv"):
        a = digest(args.first / file)
        b = digest(args.same / file)
        record("same-host", "all", file, a, b, a == b)
    a = digest(args.first / "positions.csv")
    b = digest(args.other / "positions.csv")
    record("cross-host", "all", "positions.csv", a, b, a == b)
    required = {f"{label} {cap}" for label in ("full", "W", "CRP subsample") for cap in (3, 5)}
    a_rows, b_rows = [], []
    for label, run, rows in (("first", args.first, a_rows), ("other", args.other, b_rows)):
        rows.extend(read_csv(run / "fits.csv"))
        labels = [row["case"] for row in rows]
        record("cross-host", label, "complete fit roster", sorted(required), sorted(labels),
               len(labels) == len(required) and set(labels) == required)

    a_rows = {r["case"]: r for r in a_rows}
    b_rows = {r["case"]: r for r in b_rows}
    endpoint_diffs, p_diffs, share_diffs = [], [], []
    for case in sorted(required.intersection(a_rows, b_rows)):
        a, b = a_rows[case], b_rows[case]
        for field in ("shape", "transform", "components", "first_status", "final_status", "method",
                      "lrt_flag", "lower_flag", "lower_skipped"):
            record("cross-host", case, field, a[field], b[field], a[field] == b[field])
        # These synthetic CRP fixtures use a zero lower reporting bound.
        lower = 0.0 if a["shape"] == "upper" else float(a["raw_lower"])
        width = float(a["raw_upper"]) - lower
        valid_span = math.isfinite(width) and width > 0
        record("cross-host", case, "reference span > 0", "", width, valid_span)
        for field in ("raw_lower", "raw_upper", "lower", "upper", "lrt_p"):
            left, right = a[field], b[field]
            if not left or not right:
                record("cross-host", case, field + " missingness", left, right, left == right)
                continue
            diff = abs(float(left) - float(right))
            is_p = field == "lrt_p"
            limit = .0001 if is_p else width * .001
            record("cross-host", case, field, left, right,
                   math.isfinite(diff) and valid_span and diff <= limit, diff, limit)
            if is_p:
                p_diffs.append(diff)
            elif valid_span:
                endpoint_diffs.append(diff / width)
        if a["share"] and b["share"]:
            share_diffs.append(100 * abs(float(a["share"]) - float(b["share"])))

    failures = [r for r in records if not r["passed"]]
    for name, rows in (("comparisons.csv", records), ("exceptions.csv", failures)):
        with (args.out / name).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(rows)
    summary = {"comparisons": len(records), "failures": len(failures),
               "max_cross_host_endpoint_difference_per_reference_width": max(endpoint_diffs, default=None),
               "max_cross_host_absolute_lrt_p_difference": max(p_diffs, default=None),
               "max_cross_host_share_difference_percentage_points": max(share_diffs, default=None),
               "hosts": hosts}
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return int(bool(failures))


if __name__ == "__main__":
    raise SystemExit(main())
