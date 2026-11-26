# SPDX-License-Identifier: GPL-3.0-only
"""Frozen block roster, verified Actions transfers, and immutable release bundles."""
import argparse
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess

HERE = Path("anchoredRI/tools/subsample/confirmation")
TAG = "task-h-confirmation-50000"
PARENT = HERE.parent
RECORDS = {"fits.csv", "routes.csv", "receipt.csv", "models.rds", "positions.rds", "SHA256SUMS"}


def original(script):
    path = PARENT / (script + ".py")
    spec = importlib.util.spec_from_file_location(script, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


transfer = original("full-transfer")
bundle = original("full-bundle")


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def run(*args):
    return subprocess.check_output([str(a) for a in args], text=True).strip()


def roster():
    with (HERE / "roster.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    ids = [int(row["j"]) for row in rows]
    expected = [*range(481, 961), *range(1441, 1921), *range(2401, 2881)]
    blocks = [int(row["block"]) for row in rows]
    if ids != expected or blocks != [i // 8 + 1 for i in range(1440)]:
        raise ValueError("The roster must contain the original 1,440 eligible j values in 180 consecutive blocks")
    return [ids[i:i + 8] for i in range(0, len(ids), 8)]


def group(block):
    return 1 + (block - 1) // 23


def fixture_name(host, attempt):
    return f"confirm-fixtures-{host}-a{attempt}"


def asset_name(block, attempt):
    return f"confirm-block-g{group(block):02d}-b{block:03d}-a{attempt}"


def earlier(assets, make_name):
    # Only the current workflow run is eligible. The REST inventory is
    # paginated, including when 1,440 individual checkpoints are present.
    names = {asset["name"]: asset for asset in assets}
    for attempt in range(int(os.environ["GITHUB_RUN_ATTEMPT"]), 0, -1):
        if make_name(attempt) in names:
            return names[make_name(attempt)]
    return None


def record_check(directory, j):
    records = directory / "records" / str(j)
    with (records / "SHA256SUMS").open() as handle:
        sums = dict((name, digest) for digest, name in (s.rstrip("\n").split("  ", 1) for s in handle))
    if set(sums) != RECORDS - {"SHA256SUMS"} or {p.name for p in records.iterdir()} != RECORDS:
        raise ValueError(f"Incomplete saved model files at j={j}")
    for name, digest in sums.items():
        if sha(records / name) != digest:
            raise ValueError(f"Saved case file checksum mismatch at j={j}, {name}")
    verified = directory / "verification" / str(j)
    if {p.name for p in verified.iterdir()} != {"records.csv", "hashes.csv"}:
        raise ValueError(f"Incomplete independent readback for j={j}")
    with (verified / "records.csv").open(newline="") as handle:
        verification, = csv.DictReader(handle)
    if (int(verification["j"]) != j or verification["verified"] != "TRUE"
            or [int(verification[k]) for k in ("n_models", "n_draws", "n_positions", "n_routes")]
            != [101, 100, 500000, 2]):
        raise ValueError(f"Wrong independently verified counts for j={j}")
    with (verified / "hashes.csv").open(newline="") as handle:
        reference = list(csv.DictReader(handle))
    if len(reference) != len(RECORDS) or {r["file"] for r in reference} != RECORDS:
        raise ValueError(f"Missing independent file digests for j={j}")
    for row in reference:
        file = records / row["file"]
        if int(row["j"]) != j or int(row["bytes"]) != file.stat().st_size or row["sha256"] != sha(file):
            raise ValueError(f"Independent digest or length differs for j={j}: {file.name}")
    with (records / "receipt.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 1 or int(rows[0]["j"]) != j:
        raise ValueError(f"Case receipt identifies a different j={j}")
    envelope = json.loads((directory / "job.json").read_text())
    if (envelope["commit"] != os.environ["GITHUB_SHA"] or str(envelope["run"]) != os.environ["GITHUB_RUN_ID"]
            or envelope["j"] != j):
        raise ValueError(f"Checkpoint did not originate in this run at j={j}")
    return records, verified


def restore(block):
    assets = transfer.artifacts()
    for j in roster()[block - 1]:
        asset = earlier(assets, lambda a: f"confirm-case-j{j}-a{a}")
        if asset:
            dst = Path("blockdata") / str(j)
            transfer.download(asset, dst, {"records", "verification", "job.json"})
            record_check(dst, j)
            print(f"Restored verified j={j} from {asset['name']}")


def release():
    receipt = json.loads((HERE / "launch.json").read_text())
    repo = os.environ["GITHUB_REPOSITORY"]
    endpoint = f"repos/{repo}/releases/{int(receipt['release_id'])}"
    remote = json.loads(run("gh", "api", endpoint))
    if (remote["tag_name"] != TAG or receipt["release_tag"] != TAG or not remote["draft"]
            or remote["target_commitish"] != receipt["implementation_commit"]):
        raise ValueError("Release identity, frozen implementation target, or draft state differs")
    return repo, endpoint


def validate():
    receipt = json.loads((HERE / "launch.json").read_text())
    if run("git", "rev-parse", "HEAD^") != receipt["implementation_commit"]:
        raise ValueError("Launch must directly follow the frozen implementation")
    changed = run("git", "diff", "--name-only", "HEAD^", "HEAD").splitlines()
    if changed != [str(HERE / "launch.json")]:
        raise ValueError("Only the launch receipt may change after the implementation freeze")
    manifest = HERE / "SOURCE.sha256"
    if sha(manifest) != receipt["source_manifest_sha256"]:
        raise ValueError("Frozen source manifest checksum differs")
    subprocess.run(["sha256sum", "--quiet", "-c", str(manifest)], check=True)
    release()
    blocks = roster()
    matrix = {"include": [dict(block=b, padded=f"{b:03d}", group=f"{group(b):02d}")
                          for b in range(1, len(blocks) + 1)]}
    with open(os.environ["GITHUB_OUTPUT"], "a") as out:
        print("matrix=" + json.dumps(matrix, separators=(",", ":")), file=out)


def download_inputs():
    repo = os.environ["GITHUB_REPOSITORY"]
    asset_id = int(os.environ["INPUT_ARTIFACT_ID"])
    asset = json.loads(run("gh", "api", f"repos/{repo}/actions/artifacts/{asset_id}"))
    if (int(asset["workflow_run"]["id"]) != int(os.environ["GITHUB_RUN_ID"])
            or asset["workflow_run"]["head_sha"] != os.environ["GITHUB_SHA"]
            or not asset["name"].startswith("confirm-inputs-a")):
        raise ValueError("Input metadata belongs to another launch")
    transfer.download(asset, Path("inputs"))
    subprocess.run(["sha256sum", "--quiet", "-c", "INPUT-SHA256SUMS"], cwd="inputs", check=True)


def fixtures():
    assets = transfer.artifacts()
    for host in ("primary", "other"):
        asset = earlier(assets, lambda a: fixture_name(host, a))
        if not asset:
            raise ValueError("Missing fixture host " + host)
        transfer.download(asset, Path("fixtures") / host)
    identities = [(Path("fixtures") / host / "host.json").read_text() for host in ("primary", "other")]
    if json.loads(identities[0])["boot_id"] == json.loads(identities[1])["boot_id"]:
        raise ValueError("Cross-host fixtures were executed on the same kernel boot")


def make_group(g):
    assets = transfer.artifacts()
    planned = roster()
    entries = {}
    for block, ids in enumerate(planned, 1):
        if group(block) != g:
            continue
        asset = earlier(assets, lambda a: asset_name(block, a))
        if not asset:
            raise ValueError(f"Missing completed snapshot for block {block}")
        root = Path("incoming") / str(block)
        transfer.download(asset, root, {str(j) for j in ids})
        if {p.name for p in root.iterdir()} != {str(j) for j in ids}:
            raise ValueError(f"Unexpected case membership in block {block}")
        for j in ids:
            case = root / str(j)
            records, verified = record_check(case, j)
            for kind, folder in (("records", records), ("verification", verified)):
                for path in folder.iterdir():
                    entries[f"{kind}/{j}/{path.name}"] = path
            entries[f"jobs/{j}/job.json"] = case / "job.json"
    return entries


def publish(label):
    if label == "verification":
        entries = {str(p): p for directory in ("inputs", "fixtures")
                   for p in Path(directory).rglob("*") if p.is_file()}
        if not entries or not Path("fixtures/checks.txt").is_file() or not any(Path("inputs").iterdir()):
            raise ValueError("Missing fixture audit or original input metadata")
        name = "task-h-confirmation-verification.tar"
    else:
        entries = make_group(int(label))
        name = f"task-h-confirmation-records-{int(label):02d}.tar"
    dst = Path("bundles")
    bundle.make_asset(name, dst, entries)
    repo, endpoint = release()
    for path in (dst / name, dst / (name + ".sha256")):
        def remote_asset():
            rows = json.loads(run("gh", "api", endpoint + "/assets?per_page=100"))
            return next((item for item in rows if item["name"] == path.name), None)
        item = remote_asset()
        if item is None:
            subprocess.run(["gh", "release", "upload", TAG, str(path), "--repo", repo], check=True)
            item = remote_asset()
        if (item is None or item.get("state") != "uploaded" or item.get("size") != path.stat().st_size
                or item.get("digest") != "sha256:" + sha(path)):
            raise ValueError(f"Public asset digest, length or upload status differs for {path.name}")
        print(f"{path.name}: {item['digest']}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_subparsers(dest="action", required=True)
    for command in ("validate", "inputs", "fixtures"):
        actions.add_parser(command)
    actions.add_parser("restore").add_argument("block", type=int, choices=range(1, 181))
    case = actions.add_parser("case")
    case.add_argument("block", type=int, choices=range(1, 181))
    case.add_argument("slot", type=int, choices=range(1, 9))
    actions.add_parser("seal").add_argument("j", type=int)
    actions.add_parser("publish").add_argument("label", choices=["verification", *[f"{i:02d}" for i in range(1,9)]])
    args = parser.parse_args()
    if args.action == "validate":
        validate()
    elif args.action == "inputs":
        download_inputs()
    elif args.action == "fixtures":
        fixtures()
    elif args.action == "restore":
        restore(args.block)
    elif args.action == "case":
        print(roster()[args.block - 1][args.slot - 1])
    elif args.action == "seal":
        j = args.j
        root = Path("blockdata") / str(j)
        job = dict(j=j, run=os.environ["GITHUB_RUN_ID"], commit=os.environ["GITHUB_SHA"],
                   attempt=int(os.environ["GITHUB_RUN_ATTEMPT"]), runner=os.environ["RUNNER_NAME"])
        (root / "job.json").write_text(json.dumps(job, indent=2) + "\n")
        record_check(root, j)
    elif args.action == "publish":
        publish(args.label)


if __name__ == "__main__":
    main()
