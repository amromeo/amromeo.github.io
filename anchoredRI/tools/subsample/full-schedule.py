"""GitHub Actions slot checkpoints, record verification and immutable release assets."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

STUDY = Path("anchoredRI/tools/subsample")
RELEASE = "task-h-subsample-50000"


def sha(path):
    with path.open("rb") as source:
        return hashlib.file_digest(source, "sha256").hexdigest()


def roster():
    with (STUDY / "full-roster.csv").open() as source:
        return list(csv.DictReader(source))


def groups(block):
    return 1 + (block - 1) // 23


def matrix():
    rows = roster()
    ids = [int(row["Index"]) for row in rows]
    with (STUDY / "pilot-roster.csv").open() as source:
        pilots = {int(row["Index"]) for row in csv.DictReader(source)}
    if len(ids) != 1420 or ids != sorted(set(ids)) or len(pilots) != 20 or pilots.intersection(ids):
        raise ValueError("Expected the 1,420 distinct original nonpilot indices, in order")
    for i, row in enumerate(rows):
        if (int(row["block"]), int(row["slot"])) != (i // 8 + 1, i % 8 + 1):
            raise ValueError("Invalid block or slot order")
    return {"include": [dict(block=i, padded=f"{i:03d}", group=f"{groups(i):02d}")
                        for i in range(1, 179)]}


def atomic_json(path, obj):
    temp = path.with_suffix(".tmp")
    temp.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n")
    temp.replace(path)


def utc():
    return datetime.now(timezone.utc).isoformat(timespec="milliseconds")


def run(*args, capture=False):
    return subprocess.run([str(v) for v in args], capture_output=capture, text=True, check=True)


def receipt(mode, case, idx):
    return run(sys.executable, STUDY / "full-worker-check.py", mode, case, idx, capture=True).stdout.strip()


def key(block, slot):
    hits = [r for r in roster() if (int(r["block"]), int(r["slot"])) == (block, slot)]
    if len(hits) != 1:
        raise ValueError(f"No unique dataset at block {block}, slot {slot}")
    return dict(schema=1, run_id=int(os.environ["GITHUB_RUN_ID"]), commit=os.environ["GITHUB_SHA"],
                group=groups(block), block=block, slot=slot, Index=int(hits[0]["Index"]))


def prior_state(fields):
    block, group, slot = fields["block"], fields["group"], fields["slot"]
    prefix = f"task-h-g{group:02d}-b{block:03d}-s{slot}-a"
    found = []
    for folder in Path("prior").glob(f"{prefix}*"):
        info = json.loads((folder / "slot.json").read_text())
        if (any(info.get(k) != v for k, v in fields.items())
                or folder.name != prefix + str(info["attempt"])
                or not 0 < info["attempt"] < int(os.environ["GITHUB_RUN_ATTEMPT"])):
            raise ValueError(f"Checkpoint identity mismatch in {folder}")
        case = folder / "records" / str(fields["Index"])
        if info["state"] == "error":
            raise ValueError(f"Study or verification error in {folder}; review before resumption")
        if info["state"] not in ("complete", "partial"):
            raise ValueError(f"Invalid state in {folder}")
        if info.get("case_receipt_sha256"):
            if sha(case / "case-receipt.csv") != info["case_receipt_sha256"]:
                raise ValueError(f"Receipt digest mismatch in {folder}")
            state = receipt("verify", case, fields["Index"])
            if state != ("complete" if info.get("awaiting_verification") else info["state"]):
                raise ValueError(f"Checkpoint state differs in {folder}")
        elif info["state"] != "partial" or case.exists():
            raise ValueError(f"Missing case receipt in {folder}")
        if info["state"] == "complete" and not (folder / "verification" / "verified-datasets.csv").is_file():
            raise ValueError(f"Missing independent verification for {folder}")
        found.append((info["attempt"], info["state"], folder))
    if sum(state == "complete" for _, state, _ in found) > 1:
        raise ValueError(f"Duplicate completed attempts at block {block}, slot {slot}")
    found.sort(reverse=True)
    if any(state == "complete" for _, state, _ in found[1:]):
        raise ValueError("A later attempt followed a completed slot")
    return found[0] if found else (None, None, None)


def fit_slot(block, slot):
    identity = key(block, slot)
    idx = identity["Index"]
    _, prior, directory = prior_state(identity)
    if prior == "complete":
        print(f"Index {idx} already completed and verified in {directory.name}")
        return 0
    root = Path("full-slots") / f"s{slot}"
    root.mkdir(parents=True)
    record = root / "records" / str(idx)
    if prior == "partial" and (directory / "records" / str(idx)).exists():
        shutil.copytree(directory / "records" / str(idx), record)
    meta = dict(identity, attempt=int(os.environ["GITHUB_RUN_ATTEMPT"]), state="partial",
                restored_from=directory.name if directory else None, started_at=utc(), finished_at=None,
                worker_exit=None, verification_exit=None, wall_seconds=None, case_receipt_sha256=None)
    status = root / "slot.json"
    atomic_json(status, meta)
    begin = time.monotonic()
    try:
        worker = ["bash", STUDY / "full-hosted.sh", "r", STUDY / "full-worker.R", "inputs", root / "records", block, slot]
        meta["worker_exit"] = subprocess.run(list(map(str, worker))).returncode
        if meta["worker_exit"] != 0:
            # Only host interruptions are resumable from saved draws. An R or
            # audit error remains visible and stops all remaining matrix work.
            if meta["worker_exit"] not in (-15, -9, -2, 130, 137, 143):
                meta["state"] = "error"
            return 1
        if receipt("seal", record, idx) != "complete":
            raise ValueError(f"Worker returned without complete records for Index {idx}")
        meta["case_receipt_sha256"] = sha(record / "case-receipt.csv")
        verifier = ["bash", STUDY / "full-hosted.sh", "r", STUDY / "full-verify.R",
                    root / "records", "inputs", root / "verification", "full", idx]
        meta["verification_exit"] = subprocess.run(list(map(str, verifier))).returncode
        if meta["verification_exit"] != 0:
            meta["state"] = "error"
            return 1
        meta["state"] = "complete"
        return 0
    except Exception:
        meta["state"] = "error"
        raise
    finally:
        meta["finished_at"] = utc()
        meta["wall_seconds"] = time.monotonic() - begin
        atomic_json(status, meta)


def checkpoint(block, slot):
    # Also called by the composite action on cancellation. Previously uploaded
    # slots already live in Actions even if the runner itself disappears.
    root = Path("full-slots") / f"s{slot}"
    status = root / "slot.json"
    if not status.is_file():
        return
    info = json.loads(status.read_text())
    idx = info["Index"]
    if info["state"] == "complete":
        return
    case = root / "records" / str(idx)
    if case.is_dir() and (case / "input-id.rds").exists():
        state = receipt("seal", case, idx)
        info["case_receipt_sha256"] = sha(case / "case-receipt.csv")
        # A completed but unaudited interrupted slot is checked on resumption.
        if info["state"] == "partial" and state == "complete":
            info["awaiting_verification"] = True
    elif case.exists():
        # Contains no validated reference to the ordered inputs; on a later
        # run start this case afresh without using any unverified output.
        if any(f.is_file() for f in case.rglob("*")):
            raise ValueError(f"Files without input identity at {case}")
        case.rmdir()
    atomic_json(status, info)


def gh(*args):
    return run("gh", *args, capture=True).stdout


def release_assets(repo, release_id):
    return {a["name"]: a for a in json.loads(gh("api", f"repos/{repo}/releases/{release_id}/assets?per_page=100"))}


def upload_group(group):
    repo = os.environ["GITHUB_REPOSITORY"]
    launch = json.loads((STUDY / "launch-full.json").read_text())
    release = json.loads(gh("release", "view", RELEASE, "--repo", repo, "--json", "apiUrl,isDraft,tagName,targetCommitish"))
    if release["targetCommitish"] != launch["implementation_commit"]:
        raise ValueError("Draft release target differs from the frozen implementation")
    if release["tagName"] != RELEASE or not release["isDraft"]:
        raise ValueError("Expected the designated draft Task H release")
    release_id = int(release["apiUrl"].rsplit("/", 1)[1])
    if release_id != 459667772:
        raise ValueError("Unexpected Task H release ID")
    bucket = Path("bundles")
    tar = bucket / f"task-h-records-{group}.tar"
    digest_file = tar.with_suffix(".tar.sha256")
    digest, filename = digest_file.read_text().strip().split()
    if group == "pilot" and digest != launch["pilot_asset_sha256"]:
        raise ValueError("Pilot bundle differs from the pre-launch digest")
    if digest != sha(tar) or filename.lstrip("*") != tar.name:
        raise ValueError(f"Invalid tar digest file {digest_file}")
    for path in (tar, digest_file):
        expected = sha(path)
        assets = release_assets(repo, release_id)
        if path.name not in assets:
            run("gh", "release", "upload", RELEASE, path, "--repo", repo)
            assets = release_assets(repo, release_id)
        asset = assets[path.name]
        if int(asset["size"]) != path.stat().st_size or asset["state"] != "uploaded":
            raise ValueError(f"Release asset size or state differs for {path.name}")
        if asset.get("digest"):
            if asset["digest"] != "sha256:" + expected:
                raise ValueError(f"Existing release asset differs for {path.name}")
        else:
            with tempfile.TemporaryDirectory() as temp:
                run("gh", "release", "download", RELEASE, "--repo", repo, "--pattern", path.name, "--dir", temp)
                if sha(Path(temp) / path.name) != expected:
                    raise ValueError(f"Downloaded release asset differs for {path.name}")
        print(f"{path.name}: sha256:{expected}")


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    commands = cli.add_subparsers(required=True, dest="command")
    commands.add_parser("matrix")
    for name in ("slot", "checkpoint"):
        p = commands.add_parser(name)
        p.add_argument("block", type=int)
        p.add_argument("slot", type=int)
    commands.add_parser("upload").add_argument("group", choices=["pilot", *(f"{i:02d}" for i in range(1, 9))])
    args = cli.parse_args()
    if args.command == "matrix":
        print(json.dumps(matrix(), separators=(",", ":")))
    elif args.command == "upload":
        upload_group(args.group)
    elif args.command == "checkpoint":
        checkpoint(args.block, args.slot)
    else:
        sys.exit(fit_slot(args.block, args.slot))
