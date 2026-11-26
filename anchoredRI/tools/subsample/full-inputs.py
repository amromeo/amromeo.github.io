#!/usr/bin/env python3
"""Restore the study's original fit inputs once per block, omitting comparator and scoring files."""

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import runpy
import zipfile

HERE = Path(__file__).resolve().parent
frozen = runpy.run_path(str(HERE / "prepare-inputs.py"))
block_roster = runpy.run_path(str(HERE / "full-roster.py"))["roster"]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("out", type=Path)
    p.add_argument("--block", type=int, required=True, choices=range(1, 179), metavar="1..178")
    p.add_argument("--anchorrefine-repo", type=Path, help="local git object store at the original pinned revision")
    args = p.parse_args()
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    repo = frozen["REFINE"]
    clones = {repo: args.anchorrefine_repo} if args.anchorrefine_repo else {}
    cached, read, check = (frozen[k] for k in ("cached","read_pinned","check"))
    archive = frozen["ARCHIVE"]

    def restore(path, name, sha):
        cached(out / name, sha, lambda: read(repo, path, clones))

    for source, path, name, sha in frozen["FIXED"]:
        if source == repo and name in ("task-g-manifest.json", "taskg-SHA256SUMS", "design-large.csv", "testsets.csv"):
            restore(path, name, sha)
    manifest = json.loads((out / "task-g-manifest.json").read_bytes())
    check(manifest["run_id"] == "45541935930" and manifest["datasets"] == 1440, "Wrong Task G archive")
    restore(archive + "task-g-index.csv", "task-g-index.csv", manifest["sha256"]["task-g-index.csv"])
    index = frozen["metadata"](out)
    roster = block_roster(out)
    check(roster == (HERE / "full-roster.csv").read_bytes(), "Block roster differs from frozen design")
    members = [r for r in csv.DictReader(io.StringIO(roster.decode())) if int(r["block"]) == args.block]
    needed = {int(r["shard"]) for r in members}
    hashes = {}
    for line in (out / "taskg-SHA256SUMS").read_text().splitlines():
        sha, name = line.split(None, 1)
        hashes[name] = sha
    objects = {m["shard"]:m for m in manifest["primary_rds"]}
    check(len(members) == (4 if args.block == 178 else 8) and needed <= objects.keys(), "Bad block/shard mapping")
    for shard in sorted(needed):
        obj = objects[shard]
        cases = [row for row in index.values() if int(row["shard"]) == shard]
        check(all(r["archive"] == obj["archive"] and r["member"] == obj["member"] for r in cases) and
              sorted(int(r["rds_position"]) for r in cases) == list(range(1, obj["datasets"] + 1)), "Invalid RDS member positions")
        def get_member():
            data = read(repo, archive + obj["archive"], clones)
            check(hashlib.sha256(data).hexdigest() == hashes[obj["archive"]], f"Wrong ZIP digest, shard {shard}")
            with zipfile.ZipFile(io.BytesIO(data)) as z:
                raw = z.read(obj["member"])
            check(len(raw) == obj["bytes"], f"Wrong RDS length, shard {shard}")
            return raw
        cached(out / f"taskg-shard-{shard}.rds", obj["rds_sha256"], get_member)
    print(f"Block {args.block}: verified {len(needed)} original RDS members for {len(members)} indices")


if __name__ == "__main__":
    main()
