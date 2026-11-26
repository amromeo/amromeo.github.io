"""Transfer verified Actions checkpoints without the toolkit's 1,000-artifact listing limit."""
import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import stat
import subprocess
import tempfile
import zipfile

BLOCK = re.compile(r"task-h-block-g(\d\d)-b(\d{3})-a([1-9]\d*)")
SLOT = re.compile(r"task-h-g(\d\d)-b(\d{3})-s([1-8])-a([1-9]\d*)")


def artifacts():
    repo = os.environ["GITHUB_REPOSITORY"]
    run_id = int(os.environ["GITHUB_RUN_ID"])
    params = ["gh", "api", "--paginate", "--slurp", f"repos/{repo}/actions/runs/{run_id}/artifacts?per_page=100"]
    pages = json.loads(subprocess.check_output(params))
    by_id = {}
    for page in pages:
        for artifact in page["artifacts"]:
            by_id[artifact["id"]] = artifact
    # Jobs may add artifacts during an interrupted-run restore. Previously
    # published checkpoints are immutable; only earlier attempts are eligible.
    return list(by_id.values())


def download(asset, dest, allowed=None):
    if asset["expired"]:
        raise ValueError(f"Expired Actions artifact {asset['name']}")
    repo = os.environ["GITHUB_REPOSITORY"]
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "payload.zip"
        with path.open("wb") as out:
            subprocess.run(["gh", "api", f"repos/{repo}/actions/artifacts/{int(asset['id'])}/zip"],
                           stdout=out, check=True)
        with path.open("rb") as result:
            digest = "sha256:" + hashlib.file_digest(result, "sha256").hexdigest()
        if digest != asset["digest"] or path.stat().st_size != asset["size_in_bytes"]:
            raise ValueError(f"GitHub ZIP checksum or length differs for {asset['name']}")
        dest.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(path) as archive:
            seen = set()
            for member in archive.infolist():
                name = PurePosixPath(member.filename)
                if (not name.parts or ".." in name.parts or name.is_absolute() or name in seen
                        or stat.S_ISLNK(member.external_attr >> 16) or (allowed and name.parts[0] not in allowed)):
                    raise ValueError(f"Unexpected member of {asset['name']}: {name}")
                seen.add(name)
                if (dest / name).exists() and not member.is_dir():
                    raise ValueError(f"Duplicate payload file {name}")
            archive.extractall(dest)


def prior(block):
    group = (block - 1) // 23 + 1
    attempt = int(os.environ["GITHUB_RUN_ATTEMPT"])
    for asset in artifacts():
        match = SLOT.fullmatch(asset["name"])
        if match:
            g, b, _, a = map(int, match.groups())
            if (g, b) == (group, block) and a < attempt:
                download(asset, Path("prior") / asset["name"], {"records", "verification", "slot.json"})


def snapshot(block):
    group = (block - 1) // 23 + 1
    attempt = int(os.environ["GITHUB_RUN_ATTEMPT"])
    dst = Path("block-snapshot")
    dst.mkdir()
    for slot in range(1, 5 if block == 178 else 9):
        prefix = f"task-h-g{group:02d}-b{block:03d}-s{slot}-a"
        current = Path("full-slots") / f"s{slot}"
        options = [(attempt, current)] if current.is_dir() else []
        options.extend((int(p.name.removeprefix(prefix)), p) for p in Path("prior").glob(prefix + "*"))
        completed = [(a, p) for a, p in options if json.loads((p / "slot.json").read_text())["state"] == "complete"]
        if len(completed) != 1:
            raise ValueError(f"Expected one completed checkpoint for block {block}, slot {slot}")
        a, folder = completed[0]
        shutil.copytree(folder, dst / (prefix + str(a)), copy_function=os.link)


def group_files(group):
    latest = {}
    for item in artifacts():
        match = BLOCK.fullmatch(item["name"])
        if not match:
            continue
        g, block, attempt = map(int, match.groups())
        if g == group and (block not in latest or attempt > latest[block][0]):
            latest[block] = (attempt, item)
    expected = set(range((group-1)*23+1, min(group*23, 178)+1))
    if set(latest) != expected:
        raise ValueError(f"Wrong block snapshots in group {group}, missing {sorted(expected-set(latest))}")
    for block in sorted(latest):
        attempt, item = latest[block]
        slots = range(1, 5 if block == 178 else 9)
        # A resumed block may contain checkpoints from different run attempts.
        allowed = {f"task-h-g{group:02d}-b{block:03d}-s{s}-a{a}"
                   for s in slots for a in range(1, attempt+1)}
        download(item, Path("incoming"), allowed)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)
    for action in ("prior", "snapshot"):
        sub.add_parser(action).add_argument("block", type=int, choices=range(1, 179))
    sub.add_parser("group").add_argument("group", type=int, choices=range(1, 9))
    sub.add_parser("smoke").add_argument("name")
    args = parser.parse_args()
    if args.cmd == "smoke":
        found = [a for a in artifacts() if a["name"] == args.name]
        if len(found) != 1:
            raise ValueError("Missing or duplicated preflight artifact")
        download(found[0], Path("preflight/returned"))
    elif args.cmd == "group":
        group_files(args.group)
    elif args.cmd == "prior":
        prior(args.block)
    else:
        snapshot(args.block)
