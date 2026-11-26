# SPDX-License-Identifier: GPL-3.0-only
"""Create the confirmation source manifest and a receipt for its committed revision."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

study = Path("anchoredRI/tools/subsample")
directory = study / "confirmation"
manifest = directory / "SOURCE.sha256"
receipt = directory / "launch.json"
runtime = study / "results/pilot-45854822638/runtime/runtime.tgz"
approved = "140bd4372d9741c70c76999d2b47e998fac17e9d"
image = "rocker/r-ver@sha256:5e8963d73b87865e0ad779df4202555fb1b436ee95dc7d7442450308b8ea77f3"


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def git(*args):
    return subprocess.check_output(["git", *args], text=True).strip()


def freeze():
    # Package code, compiled library, and dependencies are identical to the
    # original Task H fit. The workflow, generator adapter and report are new.
    pilot_files = (study / "FULL-SOURCE.sha256").read_text().splitlines()
    files = set()
    for line in pilot_files:
        original, name = line.split("  ", 1)
        files.add(Path(name))
        if sha(Path(name)) != original:
            raise ValueError(f"Original Task H source changed at {name}")
    files.update({study / "FULL-SOURCE.sha256", runtime,
                  Path(".github/workflows/anchoredri-task-h-confirm.yml")})
    files.update(path for path in directory.rglob("*") if path.is_file() and
                 path not in (manifest, receipt) and path.suffix != ".pyc"
                 and "__pycache__" not in path.parts)
    return "".join(f"{sha(p)}  {p}\n" for p in sorted(files))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["manifest", "receipt"])
    parser.add_argument("--release-id", type=int)
    args = parser.parse_args()
    if sha(runtime) != "4ae1bb45031aca048062684b0700fcc9ba723906141508c75796e8868afb98ce":
        raise ValueError("Compiled Task H runtime changed")
    if git("rev-parse", f"{approved}:{directory}/PROTOCOL.md") != git("hash-object", str(directory / "PROTOCOL.md")):
        raise ValueError("Confirmation protocol differs from the approved revision")
    for line in (study / "package-trees.txt").read_text().splitlines():
        expected, name = line.split(" ", 1)
        if git("rev-parse", f"HEAD:{name}") != expected:
            raise ValueError(f"Package source changed at {name}")
    lines = freeze()
    if args.command == "manifest":
        manifest.write_text(lines)
        return
    if git("status", "--porcelain") or manifest.read_text() != lines:
        raise ValueError("Launch requires a clean, committed source freeze")
    if not args.release_id or args.release_id < 1:
        raise ValueError("Provide the reserved GitHub data release ID")
    receipt.write_text(json.dumps(dict(
        implementation_commit=git("rev-parse", "HEAD"), protocol_commit=approved,
        source_manifest_sha256=sha(manifest), roster_sha256=sha(directory/"roster.csv"),
        container=image, runtime_sha256=sha(runtime), r_version="4.5.0",
        library_lock_sha256=sha(study/"library-lock.csv"), anchoredRI_version="0.4.0",
        generator_commit="c043917226b20929296b80cda53e0ae293968b42",
        evidence_commit="43b97a85df790c1d6b8e01466796a6c72f145251",
        release_tag="task-h-confirmation-50000", release_id=args.release_id,
        datasets=1440, N=50000, draws_per_bank=50, banks=["A","B"],
        blocks=180, cases_per_block=8, maximum_parallel_fit_jobs=20,
        fitting_jobs=180, fixture_hosts=2, release_groups=8,
        order="ascending original j in signed design at N = 50000",
        stop_after="confirmation report"), indent=2)+"\n")


if __name__ == "__main__":
    main()
