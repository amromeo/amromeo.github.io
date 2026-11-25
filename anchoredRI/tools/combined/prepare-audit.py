"""Restore one shard of the published Task G audit and its expected results."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
from hashlib import sha256
from io import BytesIO, TextIOWrapper
import json
from pathlib import Path
import time
from urllib.request import urlopen
from zipfile import ZipFile

GITHUB = "https://raw.githubusercontent.com/amromeo/"
STAGE3 = GITHUB + "anchorRefine/638e2d3500faa01d4e2816895a95d7a6685567a0/"
CONFIRMATION = STAGE3 + "results/stage3-confirmation-45463873162/"
TASK_G = GITHUB + "amromeo.github.io/655a007f6d5ec67d41262c45f21293b3d0b945b5/anchoredRI/tools/combined/results/task-g-evidence.zip"


def download(url, checksum=None):
    for attempt in range(4):
        try:
            with urlopen(url, timeout=90) as response:
                data = response.read()
            if checksum and sha256(data).hexdigest() != checksum:
                raise ValueError(f"SHA-256 differs for {url}")
            return data
        except OSError:
            if attempt == 3:
                raise
            time.sleep(2 ** attempt)


def index_by_id(rows):
    index = {}
    for row in rows:
        case_id = int(row["j"])
        if case_id in index:
            raise ValueError(f"Duplicate ID {case_id}")
        index[case_id] = row
    return index


def evidence_rows(archive, prefix):
    for name in archive.namelist():
        if name.startswith("art/taskg-audit-") and Path(name).name.startswith(prefix):
            with archive.open(name) as handle:
                yield from csv.DictReader(TextIOWrapper(handle, encoding="utf-8"))


def prepare(folder, group, groups):
    folder.mkdir(parents=True, exist_ok=True)
    manifest = json.loads(download(CONFIRMATION + "task-g-manifest.json",
                                   "53c6bab239a025c6d8c25039ba55f48cfff8070accf89a941d98c50679ecc784"))
    sums = dict(line.split()[::-1] for line in download(CONFIRMATION + "SHA256SUMS").decode().splitlines())
    index = download(CONFIRMATION + "task-g-index.csv", manifest["sha256"]["task-g-index.csv"])
    index = index_by_id(csv.DictReader(index.decode().splitlines()))
    design = download(STAGE3 + "stage3/design-confirmation.csv",
                      "f3f6f4e70dd1d06b50fb408bef03a41f0a2cf3a0cb18537a274688c3caf4b505")
    design = index_by_id(csv.DictReader(design.decode().splitlines()))
    with ZipFile(BytesIO(download(TASK_G, "989328da9c8cbc7c4d5772d986420e32dca97e7f2800856c4ca1927da244785c"))) as evidence:
        ids = index_by_id(csv.DictReader(evidence.read("audit.csv").decode().splitlines()))
        anchored = index_by_id(evidence_rows(evidence, "audit-"))
        blend = index_by_id(row for row in evidence_rows(evidence, "checks-")
                            if row["arm"] == "t0.25")
    if len(ids) != 749 or ids.keys() != anchored.keys() or ids.keys() != blend.keys():
        raise ValueError("The published audit must have exactly 749 complete fit and check records")

    archives = sorted({index[j]["archive"] for j in ids})[group::groups]
    members = {row["archive"]: row for row in manifest["primary_rds"]}

    def restore(name):
        member = members[name]
        with ZipFile(BytesIO(download(CONFIRMATION + name, sums[name]))) as archive:
            data = archive.read(member["member"])
        if sha256(data).hexdigest() != member["rds_sha256"]:
            raise ValueError(f"Primary RDS checksum differs for {name}")
        (folder / Path(name).with_suffix(".rds").name).write_bytes(data)

    with ThreadPoolExecutor(max_workers=8) as pool:
        list(pool.map(restore, archives))

    fixtures = []
    for case_id in sorted(ids):
        case = index[case_id]
        if case["archive"] not in archives:
            continue
        arm_b, post = anchored[case_id], blend[case_id]
        if any(case[key] != arm_b[key] for key in ["routine_hash", "anchor_hash"]):
            raise ValueError(f"Task G hashes disagree for dataset {case_id}")
        if arm_b["status"] == "reported":
            method = "anchoredRI"
            limits = (arm_b["lower"], arm_b["upper"])
        elif post["status"] == "keep":
            method = "anchorRefine"
            limits = (post["released_native_lower"], post["released_native_upper"])
        else:
            method = ""
            limits = ("", "")
        fixtures.append(dict(
            j=case_id, rds=Path(case["archive"]).with_suffix(".rds").name,
            position=case["rds_position"], N=case["N"], h=design[case_id]["h"],
            routine_hash=case["routine_hash"], anchor_hash=case["anchor_hash"],
            fallback_seed=design[case_id]["fitter_seed"], pit_seed=post["seed"],
            first_status=arm_b["status"], first_reason=arm_b["reason"],
            first_transform=arm_b["family"], first_components=arm_b["components"],
            check_status=post["status"], check_reason=post["reason"], check_p=post["p"],
            method=method, lower=limits[0], upper=limits[1],
        ))
    with (folder / "fixtures.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fixtures[0]))
        writer.writeheader()
        writer.writerows(fixtures)
    print(f"Restored {len(fixtures)} audit datasets from {len(archives)} verified shards")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path)
    parser.add_argument("group", type=int)
    parser.add_argument("groups", type=int)
    args = parser.parse_args()
    if not 0 <= args.group < args.groups <= 96:
        parser.error("require 0 <= group < groups <= 96")
    prepare(args.folder, args.group, args.groups)
