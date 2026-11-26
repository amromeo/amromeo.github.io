"""Freeze Task H sources and authorize the original sorted block assignments."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess

STUDY = Path("anchoredRI/tools/subsample")
WORKFLOW = Path(".github/workflows/anchoredri-task-h-full.yml")
PILOT = STUDY / "results/pilot-45854822638"
MANIFEST = STUDY / "FULL-SOURCE.sha256"
RECEIPT = STUDY / "launch-full.json"
PILOT_COMMIT = "beac1a10d6451f006a58f3ef3f97417c9838da2d"
ROSTER_DIGEST = "17e29b5a50419bfc42de3c7be5ac664640f695ed3eb1e7b43041daf1747523bc"
RUNTIME_DIGEST = "4ae1bb45031aca048062684b0700fcc9ba723906141508c75796e8868afb98ce"
IMAGE = "rocker/r-ver@sha256:5e8963d73b87865e0ad779df4202555fb1b436ee95dc7d7442450308b8ea77f3"
RUNTIME = PILOT / "runtime/runtime.tgz"
RELEASE = "task-h-subsample-50000"
FULL_APPROVAL = "9a1cf84efe3cffc40d13bcd6eaf3a6ba776242ef"


def git(*args):
    return subprocess.check_output(["git", *args], text=True).strip()


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def sources():
    files = {p for p in STUDY.iterdir() if p.is_file() and p not in (RECEIPT, MANIFEST)}
    for line in (STUDY / "SOURCE.sha256").read_text().splitlines():
        original, name = line.split("  ", 1)
        path = Path(name)
        # The full launch changes only the protocol and pilot cohort entry point
        # among the original scientific and verification sources.
        if path != STUDY / "PROTOCOL.md" and path != STUDY / "pilot.R" and sha(path) != original:
            raise ValueError(f"Frozen pilot source changed: {path}")
        files.add(path)
    files.update({WORKFLOW, RUNTIME})
    files.update(p for p in (STUDY / "full-slot").rglob("*") if p.is_file())
    return sorted(files)


def check_cohort():
    paths = [STUDY / name for name in ("full-roster.csv", "pilot-roster.csv")]
    with paths[0].open(newline="") as f, paths[1].open(newline="") as g:
        rows = list(csv.DictReader(f))
        ids = [int(row["Index"]) for row in rows]
        pilots = {int(row["Index"]) for row in csv.DictReader(g)}
    if (len(rows) != 1420 or len(set(ids) | pilots) != 1440 or len(pilots) != 20 or ids != sorted(ids)
            or sha(paths[0]) != ROSTER_DIGEST):
        raise ValueError("Unexpected run cohort or ordering")
    if any((int(row["block"]), int(row["slot"])) != (1 + pos // 8, 1 + pos % 8) for pos, row in enumerate(rows)):
        raise ValueError("Blocks do not match ordered positions")
    if sha(RUNTIME) != RUNTIME_DIGEST:
        raise ValueError("Incorrect compiled pilot runtime")
    for line in (STUDY / "package-trees.txt").read_text().splitlines():
        old, filename = line.split(" ", 1)
        if git("rev-parse", f"HEAD:{filename}") != old:
            raise ValueError(f"Package source differs from the pilot at {filename}")


def launch(commit):
    return dict(implementation_commit=commit, pilot_commit=PILOT_COMMIT, protocol_approval=FULL_APPROVAL,
                source_manifest_sha256=sha(MANIFEST), roster="anchoredRI/tools/subsample/full-roster.csv",
                roster_sha256=ROSTER_DIGEST, reused_pilot=20, remaining=1420, order="ascending original numeric Index",
                split="177 consecutive blocks of 8, followed by one block of 4", fitting_jobs=178,
                preflight_jobs=1, bundle_jobs=8, total_jobs=187, maximum_concurrent_fitting_jobs=20,
                bundle_groups="group = 1 + floor((block - 1) / 23)", r_version="4.5.0",
                container=IMAGE, runtime_sha256=RUNTIME_DIGEST, release_tag=RELEASE,
                pilot_asset_sha256="ed6f798c91aae5ed930c628420b3a6ce66bd78987ab6044e9c1b05df5fbab15e",
                storage="GitHub release assets for records, branch summaries and SHA-256 checksums",
                stop_rule="Retry only a documented runner interruption. Stop and report any other failure.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["create", "receipt", "verify"])
    mode = parser.parse_args().mode
    check_cohort()
    manifest = "".join(f"{sha(p)}  {p}\n" for p in sources())
    if mode == "create":
        MANIFEST.write_text(manifest)
        return
    if MANIFEST.read_text() != manifest:
        raise ValueError("Study source hash mismatch")
    if mode == "receipt":
        if git("status", "--porcelain"):
            raise ValueError("Create launch receipt in clean, committed implementation checkout")
        RECEIPT.write_text(json.dumps(launch(git("rev-parse", "HEAD")), indent=2) + "\n")
        return
    parent = git("rev-parse", "HEAD^")
    if json.loads(RECEIPT.read_text()) != launch(parent) or git("diff", "--name-only", parent, "HEAD") != str(RECEIPT):
        raise ValueError("Receipt does not identify the frozen parent")
    print(f"Verified source, runtime, protocol and 178 ordered blocks at {parent}")


if __name__ == "__main__":
    main()
