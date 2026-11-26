#!/usr/bin/env python3
"""Publish or verify the ordered 1,420-dataset complement of the frozen pilot."""

import argparse
import csv
import hashlib
import io
from pathlib import Path

HERE = Path(__file__).resolve().parent
PILOT_SHA256 = "61f6d0c8aa46aaf55b8d2d19d88b7399ae68a1f923a096a060b3f84ad40d091d"
INDEX_SHA256 = "69a4307507cbd22a18941d2587f74513dc3d67fe7ad0de6f7b1c97e83b26d274"


def pinned_csv(path, sha):
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != sha:
        raise ValueError(f"Wrong SHA-256 for {path}")
    return list(csv.DictReader(io.StringIO(raw.decode("utf-8"))))


def roster(input_dir):
    index = pinned_csv(input_dir / "task-g-index.csv", INDEX_SHA256)
    pilot = pinned_csv(HERE / "pilot-roster.csv", PILOT_SHA256)
    indexed = {int(r["Index"]): r for r in index}
    seen = {int(r["Index"]) for r in pilot}
    if len(index) != 1440 or len(indexed) != 1440 or len(seen) != len(pilot) or len(seen) != 20:
        raise ValueError("Expected 1,440 unique source indices and 20 unique pilot indices")
    for r in pilot:
        orig = indexed[int(r["Index"])]
        if any(r[k] != orig[k] for k in ("Analyte", "fractionPathol", "shard", "routine_n", "routine_hash", "anchor_hash")):
            raise ValueError(f"Pilot/source mismatch at Index {r['Index']}")
    remaining = sorted(indexed.keys() - seen)
    if len(remaining) != 1420 or set(remaining) | seen != indexed.keys():
        raise ValueError("Incomplete full cohort")
    stream = io.StringIO()
    writer = csv.writer(stream, lineterminator="\n")
    writer.writerow(["Index", "block", "slot", "shard", "Analyte", "pathology_stratum", "fractionPathol",
                     "routine_n", "routine_hash", "anchor_hash"])
    for pos, idx in enumerate(remaining):
        r = indexed[idx]
        writer.writerow([idx, 1 + pos // 8, 1 + pos % 8, r["shard"], r["Analyte"],
                         "<=20%" if float(r["fractionPathol"]) <= .2 else ">20%", r["fractionPathol"],
                         r["routine_n"], r["routine_hash"], r["anchor_hash"]])
    return stream.getvalue().encode("utf-8")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input_dir", type=Path)
    p.add_argument("--write", action="store_true", help="write full-roster.csv; default verifies the published bytes")
    args = p.parse_args()
    data = roster(args.input_dir)
    target = HERE / "full-roster.csv"
    if args.write:
        target.write_bytes(data)
    elif data != target.read_bytes():
        raise ValueError("Published block roster differs from the sorted pilot complement")
    print(f"1420 remaining indices, 177 blocks of 8 and one block of 4. SHA-256 {hashlib.sha256(data).hexdigest()}")


if __name__ == "__main__":
    main()
