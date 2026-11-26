# SPDX-License-Identifier: GPL-3.0-only
"""Verify and unpack the eight fit archives and the fixture/input archive."""
import argparse
import csv
import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil
import tarfile

BASE = "task-h-confirmation-"
HERE = Path(__file__).resolve().parent


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def unpack(source, dest, allowed, cases, seen):
    with tarfile.open(source, "r:") as bundle:
        inventory = bundle.getmembers()
        members = {item.name: item for item in inventory}
        if len(members) != len(inventory) or "SHA256SUMS" not in members:
            raise ValueError(f"Duplicate entry or missing manifest in {source.name}")
        manifest = bundle.extractfile("SHA256SUMS").read().decode().splitlines()
        entries = dict(line.split("  ", 1)[::-1] for line in manifest)
        if len(entries) != len(manifest) or entries.keys() != members.keys() - {"SHA256SUMS"}:
            raise ValueError(f"Archive membership differs from checksums in {source.name}")
        used = set()
        for name, expected in entries.items():
            relative = PurePosixPath(name)
            item = members[name]
            if (relative.is_absolute() or str(relative) != name or ".." in relative.parts
                    or len(relative.parts) < 2 or relative.parts[0] not in allowed
                    or not item.isfile() or name in seen):
                raise ValueError(f"Unexpected archive entry {name}")
            if cases is not None:
                if len(relative.parts) != 3 or int(relative.parts[1]) not in cases:
                    raise ValueError(f"Case outside frozen block: {name}")
                used.add(int(relative.parts[1]))
            seen.add(name)
            dst = dest / name
            dst.parent.mkdir(parents=True, exist_ok=True)
            if not dst.exists():
                tmp = dst.with_suffix(dst.suffix + ".partial")
                with bundle.extractfile(item) as inp, tmp.open("wb") as out:
                    shutil.copyfileobj(inp, out)
                if sha(tmp) != expected:
                    raise ValueError(f"Archived member digest differs at {name}")
                tmp.replace(dst)
            elif sha(dst) != expected:
                raise ValueError(f"Previously extracted file differs at {name}")
        if cases is not None and used != cases:
            raise ValueError(f"Missing planned datasets in {source.name}")
        return len(entries)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("assets", type=Path)
    parser.add_argument("extracted", type=Path)
    args = parser.parse_args()
    args.extracted.mkdir(parents=True, exist_ok=True)
    with (HERE / "roster.csv").open() as handle:
        roster = list(csv.DictReader(handle))
    labels = ["verification"] + [f"records-{g:02d}" for g in range(1, 9)]
    summary, seen = [], set()
    for ix, label in enumerate(labels):
        source = args.assets / (BASE + label + ".tar")
        digest = sha(source)
        if source.with_suffix(".tar.sha256").read_text().strip() != f"{digest}  {source.name}":
            raise ValueError(f"Archive digest differs at {source.name}")
        cases = {int(r["j"]) for r in roster if 1+(int(r["block"])-1)//23 == ix} if ix else None
        allowed = {"records","verification","jobs"} if ix else {"inputs","fixtures"}
        files = unpack(source, args.extracted, allowed, cases, seen)
        summary.append(dict(asset=source.name, bytes=source.stat().st_size, sha256=digest,
                            member_files=files, cases=len(cases) if ix else 0))
    expected_ids = {int(r["j"]) for r in roster}
    by_run = set()
    provenance = []
    for j in sorted(expected_ids):
        env, = (args.extracted / "jobs" / str(j)).iterdir()
        if env.name != "job.json":
            raise ValueError(f"Missing job provenance for {j}")
        source = json.loads(env.read_text())
        if (source["j"] != j or source["attempt"] < 1 or not source["runner"]
                or not str(source["run"]).isdigit() or len(source["commit"]) != 40):
            raise ValueError(f"Invalid run provenance for {j}")
        by_run.add((str(source["run"]), source["commit"]))
        provenance.append(source)
    if len(by_run) != 1:
        raise ValueError("Release combines different run IDs or source commits")
    with (args.extracted/"job-provenance.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=provenance[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(provenance)
    if sum(r["cases"] for r in summary) != 1440:
        raise ValueError("Release membership is incomplete")
    with (args.extracted/"release-assets.csv").open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=summary[0],lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary)
    print(f"Read back {len(summary)} bundles, {len(seen)} member digests and 1,440 case IDs")


if __name__ == "__main__":
    main()
