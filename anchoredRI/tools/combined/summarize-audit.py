"""Require all expected Task G release routes and endpoints to reproduce."""

from collections import Counter
import csv
from pathlib import Path
import sys

rows = []
for file in Path(sys.argv[1]).rglob("audit-*.csv"):
    with file.open() as handle:
        rows.extend(csv.DictReader(handle))
if len(rows) != 749 or len({row["j"] for row in rows}) != 749:
    raise ValueError("Expected 749 unique audited datasets")
counts = Counter(row["method"] for row in rows)
if counts != {"anchoredRI": 714, "anchorRefine": 28, "": 7}:
    raise ValueError(f"Unexpected route counts: {counts}")
errors = [row for row in rows if row["error"]]
if errors:
    raise ValueError(f"{len(errors)} audit failures: {errors[:5]}")
maximum = max(float(row["max_endpoint_error"]) for row in rows)
if not 0 <= maximum <= 1e-10:
    raise ValueError(f"Endpoint difference exceeded tolerance: {maximum}")
print("749 of 749 routes and endpoint records matched.")
print("714 anchoredRI intervals, 28 anchorRefine intervals, 7 nonreturns.")
print(f"Largest absolute endpoint difference: {maximum:.6g}.")
