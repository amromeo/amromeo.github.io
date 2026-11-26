#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-only
"""Build a slim, checksummed Stage 2 reference bundle without fitting or generating cases."""
import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import runpy
import subprocess
import tempfile
import zipfile

HERE = Path(__file__).resolve().parent
FROZEN = runpy.run_path(str(HERE.parent / "prepare-inputs.py"))
CODE = ("anchorRefine", "c043917226b20929296b80cda53e0ae293968b42")
EVIDENCE = ("anchorRefine", "43b97a85df790c1d6b8e01466796a6c72f145251")
ARCHIVE = "results/stage2-confirmation-45303419260/"
MANIFEST_SHA = "8080f57ddb77b522107f99b03fa8c1881e76d30d3d742a7963c774586f5dec2f"
SOURCE_HASHES = {
    "generator.R": "b6b754278afb630ff707b154eff48e8517dd60b52f36a410db558b0401215964",
    "confirmation-design.R": "dccefeab45d0af8c07f06cbadeb511d373d8657d74a1422988894c91b920e945",
    "confirmation-generator.R": "410ae1bed94e09ad2b705440e26133eaecfcb67e800a6048f370e3462f4fc720",
    "confirmation-design.csv": "5f4a78f982c3a882c8cadb742ada0f9b191e7a4ba658651fa5a69dddc48652c7",
}


def verified(payload, expected, name):
    if hashlib.sha256(payload).hexdigest() != expected:
        raise ValueError(f"Incorrect SHA-256 at {name}")
    return payload


def roster(design):
    original = list(csv.DictReader(io.StringIO(design.decode())))
    selected = [r for r in original if r["N"] == "50000"]
    expected = [j for first in (481, 1441, 2401) for j in range(first, first + 480)]
    if len(original) != 2880 or [int(r["j"]) for r in selected] != expected:
        raise ValueError("Unexpected Stage 2 cohort")
    fields = ["j", "block", "slot", "cell_id", *(k for k in selected[0] if k != "j")]
    out = io.StringIO()
    writer = csv.DictWriter(out, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for pos, record in enumerate(selected):
        writer.writerow(dict(record, block=pos // 8 + 1, slot=pos % 8 + 1, cell_id=(int(record["j"]) - 1) // 8 + 1))
    return out.getvalue().encode()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--repo", type=Path, help="Read the pinned commits from this anchorRefine Git object store")
    parser.add_argument("--hosted", action="store_true", help="Run R in the fixed Task H container from repo root")
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    clones = {CODE: args.repo, EVIDENCE: args.repo} if args.repo else {}

    def read(revision, name, sha):
        return verified(FROZEN["read_pinned"](revision, name, clones), sha, name)

    for name, sha in SOURCE_HASHES.items():
        (out / name).write_bytes(read(CODE, "stage2/" + name, sha))
    expected = roster((out / "confirmation-design.csv").read_bytes())
    if (HERE / "roster.csv").read_bytes() != expected:
        raise ValueError("Study roster differs from the archived N=50,000 design")
    manifest_bytes = read(EVIDENCE, ARCHIVE + "SHA256SUMS", MANIFEST_SHA)
    (out / "archive-SHA256SUMS").write_bytes(manifest_bytes)
    manifest = {f: h for h, f in (line.split("  ", 1) for line in manifest_bytes.decode().splitlines())}

    receipt = {"source": CODE[1], "evidence": EVIDENCE[1], "archive_manifest_sha256": MANIFEST_SHA,
               "source_sha256": SOURCE_HASHES, "primary_members": []}
    # Stage 2 archived all 48 primary RDS files in ZIPs alongside replicate/audit data.
    # Only the original inputs' identity fields and guidance panel leave this temporary directory.
    with tempfile.TemporaryDirectory(dir=out.parent, prefix="stage2-input-members-") as work:
        work = Path(work)
        for shard in range(1, 49):
            path = f"artifacts/selected-confirmation-{shard}.zip"
            raw = read(EVIDENCE, ARCHIVE + path, manifest[path])
            member = f"primary-{shard}/results.rds"
            with zipfile.ZipFile(io.BytesIO(raw)) as archive:
                if archive.namelist().count(member) != 1:
                    raise ValueError(f"Missing or duplicated original data member {member}")
                content = archive.read(member)
            (work / f"{shard}.rds").write_bytes(content)
            receipt["primary_members"].append(dict(shard=shard, archive=path, sha256=manifest[path],
                                                   member=member, member_sha256=hashlib.sha256(content).hexdigest()))
        cmd = ["Rscript", "--vanilla"]
        path_args = [HERE / "extract-reference.R", work, out]
        if args.hosted:
            cmd = ["bash", "anchoredRI/tools/subsample/full-hosted.sh", "r"]
            path_args = [p.relative_to(Path.cwd().resolve()) for p in path_args]
        subprocess.run([*cmd, *(str(p) for p in path_args)], check=True)
    (out / "input-sources.json").write_text(json.dumps(receipt, indent=2) + "\n")
    files = sorted([*SOURCE_HASHES, "references.rds", "reference-hashes.csv", "input-sources.json", "archive-SHA256SUMS"])
    if {p.name for p in out.iterdir()} - set(files) - {"INPUT-SHA256SUMS"}:
        raise ValueError(f"Unexpected file in input bundle at {out}")
    manifest_path = out / "INPUT-SHA256SUMS"
    with tempfile.NamedTemporaryFile("w", dir=out, delete=False) as tmp:
        for name in files:
            tmp.write(hashlib.sha256((out / name).read_bytes()).hexdigest() + "  " + name + "\n")
    os.replace(tmp.name, manifest_path)
    print("Sealed original metadata and panels for 1,440 cases across 48 verified primary shards")


if __name__ == "__main__":
    main()
