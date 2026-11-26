#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-only
"""Verify downloaded confirmation fit records against the independent readback."""
import argparse
import csv
from hashlib import sha256
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODELS = {"fits.csv", "models.rds", "positions.rds", "receipt.csv", "routes.csv"}
COUNTS = {"n_models": 101, "n_draws": 100, "n_positions": 500000, "n_routes": 2}


def rows(path):
    with path.open(newline="") as source:
        return list(csv.DictReader(source))


def digest(path):
    with path.open("rb") as source:
        h = sha256()
        for chunk in iter(lambda: source.read(1 << 20), b""):
            h.update(chunk)
        return h.hexdigest()


def require(ok, detail):
    if not ok:
        raise ValueError(detail)


def csv_write(path, items):
    with path.open("w", newline="") as dest:
        writer = csv.DictWriter(dest, fieldnames=list(items[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(items)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("records", "verification", "out"):
        parser.add_argument(key, type=Path)
    a = parser.parse_args()
    roster = rows(HERE / "roster.csv")
    ids = {int(row["j"]) for row in roster}
    require(len(ids) == len(roster) == 1440, "Expected 1,440 unique roster IDs")
    for folder in (a.records, a.verification):
        children = list(folder.iterdir())
        require(all(p.is_dir() and not p.is_symlink() for p in children)
                and {p.name for p in children} == {str(i) for i in ids}, f"Incomplete directory: {folder}")
    records, hashes = [], []
    for j in sorted(ids):
        r, v = a.records / str(j), a.verification / str(j)
        required = MODELS | {"SHA256SUMS"}
        require({f.name for f in r.iterdir()} == required and
                {f.name for f in v.iterdir()} == {"records.csv", "hashes.csv"},
                f"Unexpected inventory for {j}")
        checksum_lines = (r / "SHA256SUMS").read_text().splitlines()
        saved = {name: h for h, name in (line.split("  ", 1) for line in checksum_lines)}
        require(len(checksum_lines) == len(MODELS) and saved.keys() == MODELS, f"Invalid SHA256SUMS: {j}")
        case_hashes = rows(v / "hashes.csv")
        require(len(case_hashes) == len(required) and {x["file"] for x in case_hashes} == required,
                f"Invalid independent file hash list: {j}")
        for h in case_hashes:
            f = r / h["file"]
            require(int(h["j"]) == j and not f.is_symlink() and f.is_file()
                    and 0 < f.stat().st_size == int(h["bytes"])
                    and h["sha256"] == digest(f) and (f.name == "SHA256SUMS" or h["sha256"] == saved[f.name]),
                    f"Downloaded file differs from verification: {f}")
            hashes.append(h)
        case, = rows(v / "records.csv")
        receipt, = rows(r / "receipt.csv")
        require(int(receipt["j"]) == int(case["j"]) == j and case["verified"] == "TRUE"
                and all(int(case[k]) == n for k, n in COUNTS.items())
                and len(case["routine_hash"]) == len(case["anchor_hash"]) == 64
                and all(receipt[k] == case[k] for k in ("routine_hash", "anchor_hash"))
                and all(math.isfinite(float(case[k])) and float(case[k]) >= 0 for k in
                        ("max_abs_lrt_p_difference","max_abs_limit_difference","max_abs_component_sum_error")),
                f"Invalid or incomplete verification counts for {j}")
        records.append(case)
    a.out.mkdir(parents=True, exist_ok=True)
    csv_write(a.out / "verified-datasets.csv", records)
    csv_write(a.out / "file-hashes.csv", hashes)
    summary = {"datasets": len(records), "record_files": len(hashes),
               "counts": {k: sum(int(r[k]) for r in records) for k in COUNTS},
               "technical_first_stage_fits": sum(int(r["technical_fits"]) for r in records),
               "max_abs_lrt_p_difference": max(float(r["max_abs_lrt_p_difference"]) for r in records),
               "max_abs_limit_difference": max(float(r["max_abs_limit_difference"]) for r in records),
               "max_abs_component_sum_error": max(float(r["max_abs_component_sum_error"]) for r in records),
               "output_sha256": {k: digest(a.out / k) for k in ("verified-datasets.csv", "file-hashes.csv")}}
    (a.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
