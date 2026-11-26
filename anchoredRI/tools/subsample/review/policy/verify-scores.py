# SPDX-License-Identifier: GPL-3.0-only
"""Run the pinned Task H scorer on the routed endpoints and compare the saved scores."""

import argparse
import csv
import hashlib
import math
import subprocess
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCORE = HERE.parents[1] / "pilot-score.R"
DESIGN_SHA256 = "177cd392896cb300bf4aaa35a959e46bb3ee48d510beff8cf9957f6643a66e39"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("design", type=Path, help="Frozen Task H testsets.csv")
    design = parser.parse_args().design.resolve()
    with design.open("rb") as handle:
        if hashlib.file_digest(handle, "sha256").hexdigest() != DESIGN_SHA256:
            raise ValueError("The design differs from the frozen RIbench cohort")
    with (HERE / "by-index.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    expected = {}
    with tempfile.TemporaryDirectory() as temp:
        source, scored = (Path(temp) / name for name in ("endpoints.csv", "scores.csv"))
        with source.open("w", newline="") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(["Index", "method", "lower", "upper"])
            for row in rows:
                targets = [("policy_" + row["bank"], "policy")]
                if row["bank"] == "A":
                    targets.append(("package", "package"))
                for method, prefix in targets:
                    writer.writerow([row["Index"], method, row[prefix + "_lower"], row[prefix + "_upper"]])
                    expected[row["Index"], method] = row, prefix
        subprocess.run(["Rscript", str(SCORE), str(source), str(design), str(scored)], check=True)
        stats = {}
        with scored.open(newline="") as handle:
            for actual in csv.DictReader(handle):
                method = actual["method"]
                row, prefix = expected.pop((actual["Index"], method))
                key = method, row["analyte"]
                group = stats.setdefault(key, dict(method=method, analyte=row["analyte"], assigned=0,
                                                   returns=0, zero_upper_limits=0, max_abs_z_difference=0,
                                                   max_abs_d_difference=0))
                group["assigned"] += 1
                if row[prefix + "_upper"]:
                    group["returns"] += 1
                    group["zero_upper_limits"] += int(float(row[prefix + "_upper"]) == 0)
                for field in ("z_lower", "z_upper", "d"):
                    saved = row[prefix + "_" + field]
                    if not saved and not actual[field]:
                        continue
                    if not saved or not actual[field]:
                        raise ValueError(f"Score availability differs: {key}, {row['Index']}, {field}")
                    delta = abs(float(saved) - float(actual[field]))
                    if not math.isfinite(delta) or delta > 1e-12:
                        raise ValueError(f"Score differs: {key}, {row['Index']}, {field} by {delta}")
                    col = "max_abs_d_difference" if field == "d" else "max_abs_z_difference"
                    group[col] = max(group[col], delta)
        if expected:
            raise ValueError("R scoring omitted an assignment")
    with (HERE / "score-readback.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=next(iter(stats.values())), lineterminator="\n")
        writer.writeheader()
        writer.writerows(stats[key] for key in sorted(stats))


if __name__ == "__main__":
    main()
