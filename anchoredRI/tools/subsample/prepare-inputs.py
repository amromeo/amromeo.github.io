#!/usr/bin/env python3
"""Restore pinned, checksummed Task H inputs. Does not fit or sample any data."""

import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import tempfile
import time
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen
import zipfile

REFINE = ("anchorRefine", "638e2d3500faa01d4e2816895a95d7a6685567a0")
SITE = ("amromeo.github.io", "8bcd49c3c7977c50982524b9dbc48d16bdfa5914")
ARCHIVE = "results/stage3-ribench-large-45541935930/"
E = "anchoredRI/tools/boxcox/accuracy/"
FIXED = [
    (REFINE, ARCHIVE + "task-g-manifest.json", "task-g-manifest.json", "c13112de13d49229b0846e704a7f455d0b02d15346280b651da145528662a300"),
    (REFINE, ARCHIVE + "SHA256SUMS", "taskg-SHA256SUMS", "82f439cca6e89ce93bca0e9be887a60efe30408c6f903784d46f7064b39c3978"),
    (REFINE, "stage2/ribench/design-large.csv", "design-large.csv", "e7b3b9b0673211f84534c11f2125a00eaad3295b8180772e419d6de7c0faa267"),
    (REFINE, "stage1/large/testsets.csv", "testsets.csv", "177cd392896cb300bf4aaa35a959e46bb3ee48d510beff8cf9957f6643a66e39"),
    (SITE, E + "large/results/baseline-fits.csv", "baseline-fits.csv", "cab045a054d0651dbe0b20c8e5d8602957e521284d2b98b21c4256bc26fe960a"),
    (SITE, E + "large/results/armB-fits.csv", "armB-fits.csv", "87dce1dadba2917bd2bf160beab1933bd9d58966291ccae4b8dd6011c42fdef7"),
    (SITE, E + "large/results/selected-component-parameters.csv", "selected-component-parameters.csv", "e77e0fa6e6329f168ba01227069819eedc7bd436e987b3a20437d22cbde7958c"),
    (SITE, E + "large/results/ribench-outcomes.csv", "ribench-outcomes.csv", "05badc90db1757901244616f2364fe591926b69e7d7135f6f3b640debda72716"),
    (SITE, E + "score.R", "task-e-score.R", "29c3ba640a2731dd15bb23309442baf3c0548396964aecdd1eb3b74ead91eb8a"),
    (SITE, "anchoredRI/tools/boxcox/native-scores.R", "native-scores.R", "ad46624c72bc3ef2790e47b26e193d88df8eadccbc3555fb502b797e72a087b9"),
    (REFINE, "vendor/RIbench_1.0.2.tar.gz", "RIbench_1.0.2.tar.gz", "5703783160e2ed6bef1891d44c1401cfa0303820c46eb0f3e56f213c2ae3dd35"),
]


def check(ok, message):
    if not ok:
        raise ValueError(message)


def read_pinned(repo, path, clones):
    name, commit = repo
    if repo in clones:
        return subprocess.check_output(["git", "-C", str(clones[repo]), "show", f"{commit}:{path}"])
    url = f"https://raw.githubusercontent.com/amromeo/{name}/{commit}/{path}"
    for attempt in range(4):
        try:
            request = Request(url, headers={"User-Agent": "Task-H-archive-restore"})
            with urlopen(request, timeout=120) as response:
                return response.read()
        except (URLError, TimeoutError, ConnectionError) as error:
            if isinstance(error, HTTPError) and error.code < 500 and error.code != 429:
                raise
            if attempt == 3:
                raise
            time.sleep(2 ** attempt)


def cached(path, sha, supply):
    if path.exists():
        check(hashlib.sha256(path.read_bytes()).hexdigest() == sha, f"Existing file digest differs: {path}")
        return
    payload = supply()
    check(hashlib.sha256(payload).hexdigest() == sha, f"Downloaded bytes differ: {path.name}")
    # Several workers can safely share the directory. A partial download is never visible.
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as staging:
        staging.write(payload)
    os.replace(staging.name, path)


def rows(path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def metadata(out):
    defs, design, index = [rows(out / name) for name in ("testsets.csv", "design-large.csv", "task-g-index.csv")]
    keyed = [{int(row["Index"]): row for row in data} for data in (defs, design, index)]
    check(all(len(data) == len(mapping) == 1440 for data, mapping in zip((defs, design, index), keyed)), "Cohort size or duplicate Index")
    check(keyed[0].keys() == keyed[1].keys() == keyed[2].keys(), "Index mismatch")
    for i, entry in keyed[2].items():
        ref, named = keyed[1][i], keyed[0][i]
        check(all(entry[key] == ref[key] for key in ("Analyte", "Scenario", "N", "routine_n", "routine_hash", "anchor_hash"))
              and entry["N"] == "50000" and ref["anchor_n"] == "40" and int(ref["anchor_seed"]) == 17000000 + i
              and int(entry["j"]) == i and int(entry["routine_n"]) >= 5000
              and all(named[k] == ref[k] for k in ("Analyte", "Scenario", "N"))
              and float(entry["fractionPathol"]) == float(ref["fractionPathol"]) == float(named["fractionPathol"]),
              f"Frozen design mismatch at Index {i}")
    return keyed[2]


def numbers(value):
    values = [int(i) for i in value.split(",")]
    if len(values) != len(set(values)):
        raise argparse.ArgumentTypeError("Duplicate numeric ID")
    return set(values)


def roster_bytes(index):
    groups = {}
    for i, row in index.items():
        key = (row["Analyte"], "<=20%" if float(row["fractionPathol"]) <= 0.2 else ">20%")
        rank = hashlib.sha256(f"TaskH-pilot-v1|{i}".encode("utf-8")).hexdigest()
        groups.setdefault(key, []).append((rank, i))
    check(len(groups) == 20 and all(len(group) == 72 for group in groups.values()), "Unexpected pilot stratum sizes")
    stream = io.StringIO()
    writer = csv.writer(stream, lineterminator="\n")
    writer.writerow(["Index", "Analyte", "pathology_stratum", "fractionPathol", "routine_n", "shard",
                     "eligible_datasets", "selection_sha256", "routine_hash", "anchor_hash"])
    for (analyte, stratum), group in sorted(groups.items()):
        rank, i = min(group)
        row = index[i]
        writer.writerow([i, analyte, stratum, row["fractionPathol"], row["routine_n"], row["shard"],
                         len(group), rank, row["routine_hash"], row["anchor_hash"]])
    return stream.getvalue().encode("utf-8")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("out", type=Path)
    ap.add_argument("--anchorrefine-repo", type=Path, help="local clone containing the pinned anchorRefine revision")
    ap.add_argument("--site-repo", type=Path, help="local clone containing the pinned Task E revision")
    selector = ap.add_mutually_exclusive_group()
    selector.add_argument("--metadata-only", action="store_true")
    selector.add_argument("--ids", type=numbers, help="comma-separated original RIbench indices")
    selector.add_argument("--shards", type=numbers, help="comma-separated original Task G shard numbers")
    ap.add_argument("--write-roster", type=Path, help="write the roster derived entirely from design fields")
    args = ap.parse_args()
    out = args.out
    out.mkdir(exist_ok=True, parents=True)
    clones = {r: p for r, p in ((REFINE, args.anchorrefine_repo), (SITE, args.site_repo)) if p}

    def restore(repo, path, name, sha):
        cached(out / name, sha, lambda: read_pinned(repo, path, clones))

    for source in FIXED:
        restore(*source)
    manifest = json.loads((out / "task-g-manifest.json").read_bytes())
    check(manifest["run_id"] == "45541935930" and manifest["datasets"] == 1440, "Unexpected Task G manifest")
    for path, name in (("task-g-index.csv", "task-g-index.csv"), ("evidence/fits.csv.gz", "anchorRefine-fits.csv.gz"),
                       ("evidence/checks.csv.gz", "anchorRefine-checks.csv.gz")):
        restore(REFINE, ARCHIVE + path, name, manifest["sha256"][path])
    index = metadata(out)
    roster = roster_bytes(index)
    if args.write_roster:
        args.write_roster.write_bytes(roster)
    else:
        check(Path(__file__).with_name("pilot-roster.csv").read_bytes() == roster, "Pilot roster does not match frozen design")
    archive_hash = {}
    for line in (out / "taskg-SHA256SUMS").read_text().splitlines():
        sha, name = line.split(None, 1)
        archive_hash[name] = sha
    shard_info = {item["shard"]: item for item in manifest["primary_rds"]}
    check(set(shard_info) == set(range(1, 49)), "Missing or duplicate Task G shards")
    shard_ids = set()
    if not args.metadata_only:
        if args.ids:
            check(args.ids <= index.keys(), f"Unknown indices: {args.ids - index.keys()}")
            shard_ids = {int(index[i]["shard"]) for i in args.ids}
        else:
            shard_ids = args.shards or set(shard_info)
    check(shard_ids <= shard_info.keys(), "Unknown shard")
    for s in sorted(shard_ids):
        info = shard_info[s]
        members = [row for row in index.values() if int(row["shard"]) == s]
        check(all(row["archive"] == info["archive"] and row["member"] == info["member"] for row in members)
              and sorted(int(row["rds_position"]) for row in members) == list(range(1, info["datasets"] + 1)),
              f"Task G shard {s} index differs")
        def read_member():
            raw = read_pinned(REFINE, ARCHIVE + info["archive"], clones)
            check(hashlib.sha256(raw).hexdigest() == archive_hash[info["archive"]], f"ZIP checksum differs, shard {s}")
            with zipfile.ZipFile(io.BytesIO(raw)) as archive:
                payload = archive.read(info["member"])
            check(len(payload) == info["bytes"], f"RDS length differs, shard {s}")
            return payload
        cached(out / f"taskg-shard-{s}.rds", info["rds_sha256"], read_member)
    print(f"Verified 1440 paired design records and {len(shard_ids)} selected Task G RDS members in {out}")


if __name__ == "__main__":
    main()
