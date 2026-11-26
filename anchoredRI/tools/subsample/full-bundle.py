"""Assemble durable Task H release assets from checked slot transfers."""
import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tarfile

STUDY = Path(__file__).resolve().parent
PILOT = STUDY / "results/pilot-45854822638"
NAME = re.compile(r"task-h-g(\d\d)-b(\d{3})-s([1-8])-a([1-9]\d*)")


def digest(file):
    with file.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def csv_rows(file):
    with file.open(newline="") as stream:
        return list(csv.DictReader(stream))


def make_asset(name, destination, entries):
    destination.mkdir(parents=True, exist_ok=True)
    asset = destination / name
    hashes = []
    staging = asset.with_suffix(".tmp")
    with tarfile.open(staging, "w", format=tarfile.USTAR_FORMAT) as archive:
        for relative, source in sorted(entries.items()):
            if not source.is_file() or source.is_symlink():
                raise ValueError(f"Expected a regular file: {source}")
            hashes.append(f"{digest(source)}  {relative}\n")
            info = tarfile.TarInfo(relative)
            info.size = source.stat().st_size
            info.mode = 0o644
            with source.open("rb") as stream:
                archive.addfile(info, stream)
        value = "".join(hashes).encode()
        info = tarfile.TarInfo("SHA256SUMS")
        info.size = len(value)
        info.mode = 0o644
        archive.addfile(info, io.BytesIO(value))
    checksum = digest(staging)
    if asset.exists() and digest(asset) != checksum:
        raise ValueError(f"Existing bundle has different contents: {asset}")
    staging.replace(asset)
    asset.with_suffix(".tar.sha256").write_text(f"{checksum}  {name}\n")
    print(f"{name}: {len(entries)} files, {asset.stat().st_size} bytes, SHA-256 {checksum}")


def pilot_bundle(incoming, destination):
    known = dict((path, sha) for sha, path in (s.split("  ", 1) for s in (PILOT / "SHA256SUMS").read_text().splitlines()))
    entries = {}
    for row in csv_rows(STUDY / "pilot-roster.csv"):
        key = f"pilot/data/{row['Index']}.tar.gz"
        path = (incoming if (incoming / key).is_file() else PILOT) / key
        if digest(path) != known[key]:
            raise ValueError(f"Pilot archive differs from published checksum: {key}")
        entries[key] = path
    if len(entries) != 20:
        raise ValueError("Expected the 20 archived timing pilots")
    make_asset("task-h-records-pilot.tar", destination, entries)


def fitted_bundle(group, incoming, destination):
    roster = {}
    for row in csv_rows(STUDY / "full-roster.csv"):
        batch, slot = int(row["block"]), int(row["slot"])
        if 1 + (batch - 1) // 23 == group:
            roster[batch, slot] = int(row["Index"])
    if not roster:
        raise ValueError(f"No planned data for group {group}")
    selected = {}
    for entry in incoming.iterdir():
        match = NAME.fullmatch(entry.name)
        if not match or not entry.is_dir():
            raise ValueError(f"Unexpected slot artifact: {entry}")
        found_group, batch, slot, attempt = map(int, match.groups())
        key = batch, slot
        if found_group != group or key not in roster:
            raise ValueError(f"Unexpected slot outside planned group: {entry}")
        receipt = json.loads((entry / "slot.json").read_text())
        if ((int(receipt["Index"]), int(receipt["group"]), int(receipt["block"]), int(receipt["slot"]), int(receipt["attempt"]))
                != (roster[key], group, batch, slot, attempt)):
            raise ValueError(f"Slot identity mismatch: {entry}")
        if os.getenv("GITHUB_RUN_ID") and str(receipt["run_id"]) != os.environ["GITHUB_RUN_ID"]:
            raise ValueError(f"Incorrect workflow run: {entry}")
        if os.getenv("GITHUB_SHA") and receipt["commit"] != os.environ["GITHUB_SHA"]:
            raise ValueError(f"Incorrect implementation: {entry}")
        if receipt["state"] == "partial":
            continue
        if (receipt["state"] != "complete" or receipt.get("worker_exit") != 0
                or receipt.get("verification_exit") != 0 or key in selected):
            raise ValueError(f"Error or duplicate successful attempt for {entry}")
        case_receipt = entry / "records" / str(roster[key]) / "case-receipt.csv"
        if receipt.get("case_receipt_sha256") != digest(case_receipt):
            raise ValueError(f"Slot checksum differs for {entry}")
        selected[key] = attempt, entry
    if selected.keys() != roster.keys():
        raise ValueError(f"Missing completed slots for {sorted(roster.keys() - selected.keys())}")
    entries = {}
    for key, idx in sorted(roster.items()):
        entry = selected[key][1]
        case = entry / "records" / str(idx)
        subprocess.run([sys.executable, str(STUDY / "full-worker-check.py"), "verify", str(case), str(idx)], check=True)
        if not (case / "COMPLETE").is_file():
            raise ValueError(f"Unfinished dataset {idx}")
        for file in sorted(case.rglob("*")):
            if file.is_file():
                entries[f"records/{idx}/{file.relative_to(case)}"] = file
        entries[f"verification/{idx}/slot.json"] = entry / "slot.json"
        verification = entry / "verification"
        required = {"verified-datasets.csv", "file-hashes.csv", "observations.csv", "failures.csv"}
        if {p.name for p in verification.glob("*.csv")} != required or csv_rows(verification / "failures.csv"):
            raise ValueError(f"Incomplete or failed verification for {idx}")
        verified = csv_rows(verification / "verified-datasets.csv")
        if len(verified) != 1 or (verified[0]["Index"], verified[0]["source"]) != (str(idx), "full"):
            raise ValueError(f"Independent verification refers to a different dataset: {idx}")
        hashes = csv_rows(verification / "file-hashes.csv")
        expected = {p.name for p in case.iterdir() if p.is_file()} - {"case-receipt.csv"}
        if len(hashes) != len(expected) or {r["file"] for r in hashes} != expected:
            raise ValueError(f"Missing independently verified file for {idx}")
        for r in hashes:
            if (r["Index"] != str(idx) or r["source"] != "full" or r["sha256"] != digest(case / r["file"])
                    or int(r["bytes"]) != (case / r["file"]).stat().st_size):
                raise ValueError(f"Independent file digest mismatch for {idx}: {r['file']}")
        for file in sorted(verification.rglob("*")):
            if file.is_file():
                entries[f"verification/{idx}/{file.relative_to(verification)}"] = file
    make_asset(f"task-h-records-{group:02}.tar", destination, entries)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("group", help="Group number or pilot")
    parser.add_argument("incoming", type=Path, help="Slot artifact directory or root of a restored pilot bundle")
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    if args.group == "pilot":
        pilot_bundle(args.incoming, args.destination)
    else:
        fitted_bundle(int(args.group), args.incoming, args.destination)


if __name__ == "__main__":
    main()
