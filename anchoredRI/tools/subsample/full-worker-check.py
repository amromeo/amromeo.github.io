#!/usr/bin/env python3
"""Seal or verify the durably written parts of one Task H case, without fitting."""
import argparse
import csv
import hashlib
import os
import re
import tempfile
from pathlib import Path

RECEIPT = "case-receipt.csv"
COMPLETE = ("COMPLETE", "decisions.csv", "draws.csv", "fits.csv", "input-id.rds",
            "medians.csv", "models.rds", "positions.rds", "timing.csv")
ATTEMPT = re.compile(r"(full_public|full_cap3|W|[AB](0[1-9]|[1-4][0-9]|50))[.]rds")
FIELDS = ("Index", "state", "path", "bytes", "sha256")


def require(value, message):
    if not value:
        raise ValueError(message)


def inspect(case, index):
    require(case.name == str(index) and index > 0 and case.is_dir(), "Wrong case directory or Index")
    if (case / "COMPLETE").is_file():
        state = "complete"
        require((case / "COMPLETE").read_text().splitlines() == [f"Index {index}; 204 fits, 100 draws"],
                "Invalid completion marker")
        names = list(COMPLETE)
        for name, count in (("fits.csv", 204), ("draws.csv", 100), ("medians.csv", 8), ("decisions.csv", 48)):
            with (case / name).open(newline="") as stream:
                ids = [r["Index"] for r in csv.DictReader(stream)]
            require(len(ids) == count and set(ids) == {str(index)}, f"Incomplete or mis-keyed {name}")
    else:
        state = "partial"
        # The R writer stages each attempt under a temporary name before rename.
        # Unfinished summaries and temporary writes are rebuilt from the attempts.
        names = ["input-id.rds"]
        attempts = case / "attempts"
        if attempts.exists():
            require(attempts.is_dir() and not attempts.is_symlink(), "Invalid attempts directory")
            for p in attempts.glob("*.rds"):
                require(ATTEMPT.fullmatch(p.name), f"Unexpected attempt {p.name}")
                names.append(f"attempts/{p.name}")

    rows = []
    for name in sorted(names):
        path = case / name
        require(path.is_file() and not path.is_symlink() and path.stat().st_size > 0, f"Missing payload: {name}")
        digest = hashlib.sha256()
        with path.open("rb") as stream:
            for part in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(part)
        rows.append(dict(Index=str(index), state=state, path=name, bytes=str(path.stat().st_size),
                         sha256=digest.hexdigest()))
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=("seal", "verify"))
    p.add_argument("case", type=Path)
    p.add_argument("index", type=int)
    args = p.parse_args()
    rows = inspect(args.case, args.index)
    path = args.case / RECEIPT
    if args.action == "verify" or path.exists():
        with path.open(newline="") as stream:
            reader = csv.DictReader(stream)
            require(tuple(reader.fieldnames or ()) == FIELDS, "Unexpected receipt header")
            previous = list(reader)
        require(bool(previous), "Empty receipt")
        if args.action == "verify" or previous[0]["state"] == "complete":
            require(previous == rows, "Saved files differ from the checkpoint receipt")
        elif rows[0]["state"] == "partial":
            require(all(row in rows for row in previous), "A previously saved attempt changed")
        else:
            original = next(row for row in previous if row["path"] == "input-id.rds")
            current = next(row for row in rows if row["path"] == "input-id.rds")
            require(original == dict(current, state="partial"), "The case identity changed")
    if args.action == "seal":
        # A previous partial receipt is replaced only after the new state exists.
        with tempfile.NamedTemporaryFile(mode="w", newline="", dir=args.case, prefix=".case-receipt.",
                                         delete=False) as stream:
            tmp = Path(stream.name)
            try:
                writer = csv.DictWriter(stream, FIELDS, lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
                stream.flush()
                os.fsync(stream.fileno())
            except BaseException:
                tmp.unlink(missing_ok=True)
                raise
        try:
            tmp.replace(path)
        finally:
            tmp.unlink(missing_ok=True)
    print(rows[0]["state"])


if __name__ == "__main__":
    main()
