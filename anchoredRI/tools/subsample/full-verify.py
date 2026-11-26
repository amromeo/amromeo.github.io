#!/usr/bin/env python3
"""Check the joined 20 archived + 1,420 new Task H readback receipts."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import tarfile

HERE = Path(__file__).resolve().parent
PILOT = HERE / "results/pilot-45854822638"
PILOT_MANIFEST_HASH = "e7e44be0bc2cf2364a52ecda691063f4d465ba9dfdf8d7b522bb74eaf4dc8f06"
INPUT_HASHES = {
    "task-g-manifest.json": "c13112de13d49229b0846e704a7f455d0b02d15346280b651da145528662a300",
    "design-large.csv": "e7b3b9b0673211f84534c11f2125a00eaad3295b8180772e419d6de7c0faa267",
}
FILES = {"COMPLETE", "decisions.csv", "draws.csv", "fits.csv", "input-id.rds",
         "medians.csv", "models.rds", "positions.rds", "timing.csv"}
COUNTS = {"models": 204, "draws": 100, "position_indices": 500000, "sample_seed_checks": 100,
          "medians": 8, "decisions": 48}


def require(ok, detail):
    if not ok:
        raise ValueError(detail)


def fingerprint(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows, fields=None):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields or rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def by_index(rows):
    indexed = {int(r["Index"]): r for r in rows}
    require(len(rows) == len(indexed), "Duplicate Index in dataset roster")
    return indexed


def verify_original_pilots(hashes, archives):
    manifest_path = PILOT / "SHA256SUMS"
    require(fingerprint(manifest_path) == PILOT_MANIFEST_HASH, "Pilot manifest changed")
    manifest = dict((p, s) for s, p in (line.split("  ", 1) for line in manifest_path.read_text().splitlines()))
    image = PILOT / "runtime/image.json"
    require(fingerprint(image) == manifest["runtime/image.json"], "Archived pilot image receipt changed")
    for idx, files in hashes.items():
        key = f"pilot/data/{idx}.tar.gz"
        path = archives / f"{idx}.tar.gz"
        require(fingerprint(path) == manifest[key], f"Archived pilot bundle changed: {idx}")
        with tarfile.open(path) as archive:
            for member in archive:
                require(member.isfile() and member.name in files, f"Invalid pilot member: {idx}/{member.name}")
                row = files[member.name]
                require(member.size == int(row["bytes"]), f"Pilot byte count changed: {idx}/{member.name}")
                with archive.extractfile(member) as stream:
                    actual = hashlib.file_digest(stream, "sha256").hexdigest()
                require(actual == row["sha256"], f"Original pilot bytes differ: {idx}/{member.name}")
            names = archive.getnames()
            require(len(names) == len(FILES) and set(names) == FILES, f"Incomplete pilot record: {idx}")


def verify_worker_receipt(path, idx, hashes):
    rows = read_csv(path)
    require(len(rows) == len(FILES) and {r["path"] for r in rows} == FILES, f"Incomplete worker receipt for {idx}")
    for r in rows:
        row = hashes[r["path"]]
        require(int(r["Index"]) == idx and r["state"] == "complete" and r["bytes"] == row["bytes"]
                and r["sha256"] == row["sha256"], f"Worker and independent checksums differ: {idx}/{r['path']}")


def collect(root, source, expected, design, raw):
    datasets, files, observations = {}, {}, []
    for receipt in sorted(root.rglob("file-hashes.csv")):
        group = by_index(read_csv(receipt.parent / "verified-datasets.csv"))
        require(not read_csv(receipt.parent / "failures.csv"), f"Failed scientific checks in {receipt.parent}")
        group_files = {idx: {} for idx in group}
        for r in read_csv(receipt):
            idx, name = int(r["Index"]), r["file"]
            require(idx in group_files and name in FILES and name not in group_files[idx]
                    and r["source"] == source and int(r["bytes"]) > 0 and re.fullmatch("[0-9a-f]{64}", r["sha256"]),
                    f"Invalid file record for {idx}/{name}")
            group_files[idx][name] = r
        for idx, r in group.items():
            require(idx in expected and idx not in datasets and r["source"] == source,
                    f"Duplicate, misrouted or unexpected {source} dataset {idx}")
            d = design[idx]
            require(group_files[idx].keys() == FILES and all(int(r[k]) == n for k, n in COUNTS.items()),
                    f"Incomplete {source} dataset {idx}")
            require(r["Analyte"] == d["Analyte"] and r["stratum"] == ("<=20%" if float(d["fractionPathol"]) <= .2 else ">20%")
                    and all(r[k] == d[k] for k in ("routine_hash", "anchor_hash", "routine_n")), f"Design mismatch for {idx}")
            if source == "full" and raw:
                verify_worker_receipt(raw / str(idx) / "case-receipt.csv", idx, group_files[idx])
                for name, file in group_files[idx].items():
                    path = raw / str(idx) / name
                    require(path.stat().st_size == int(file["bytes"]) and fingerprint(path) == file["sha256"],
                            f"Transferred data differs from scientific readback: {idx}/{name}")
            files[idx] = group_files[idx]
            datasets[idx] = r
        obs = read_csv(receipt.parent / "observations.csv")
        require(all(int(r["Index"]) in group for r in obs), f"Mis-keyed observations in {receipt.parent}")
        observations.extend(obs)
    require(set(datasets) == expected,
            f"{source} coverage mismatch: missing {sorted(expected - datasets.keys())}; unexpected {sorted(datasets.keys() - expected)}")
    return datasets, files, observations


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("inputs", "pilot_verification", "full_verification", "out"):
        parser.add_argument(name, type=Path)
    parser.add_argument("--pilot-archives", type=Path, default=PILOT / "pilot/data",
                        help="Directory containing the original <Index>.tar.gz files from the pilot release")
    parser.add_argument("--records-root", type=Path,
                        help="Read every transferred full-run file under RECORDS_ROOT/<Index>/ and reconcile worker receipts")
    parser.add_argument("--pilot-only", action="store_true", help="Dry check the 20 archived pilots before collection")
    a = parser.parse_args()
    for name, sha in INPUT_HASHES.items():
        require(fingerprint(a.inputs / name) == sha, f"Wrong {name}")
    manifest = json.loads((a.inputs / "task-g-manifest.json").read_bytes())
    require(fingerprint(a.inputs / "task-g-index.csv") == manifest["sha256"]["task-g-index.csv"], "Task G dataset index changed")
    design = by_index(read_csv(a.inputs / "design-large.csv"))
    indexed = by_index(read_csv(a.inputs / "task-g-index.csv"))
    require(len(indexed) == 1440 and indexed.keys() == design.keys(), "Expected exactly 1,440 original indices")
    require(all(all(d[key] == indexed[i][key] for key in ("routine_hash", "anchor_hash", "routine_n"))
                for i, d in design.items()), "Task G index and original design disagree")
    pilot = HERE / "pilot-roster.csv"
    require(fingerprint(pilot) == "61f6d0c8aa46aaf55b8d2d19d88b7399ae68a1f923a096a060b3f84ad40d091d",
            "Wrong frozen pilot roster")
    pilot_ids = set(by_index(read_csv(pilot)))
    full_ids = set(by_index(read_csv(HERE / "full-roster.csv")))
    require(not pilot_ids & full_ids and pilot_ids | full_ids == set(design)
            and len(full_ids) == 1420 and len(pilot_ids) == 20, "Pilot and full roster do not partition the cohort")

    rows, files, obs = collect(a.pilot_verification, "pilot", pilot_ids, design, None)
    verify_original_pilots(files, a.pilot_archives)
    if not a.pilot_only:
        others, more_files, more_obs = collect(a.full_verification, "full", full_ids, design, a.records_root)
        rows.update(others)
        files.update(more_files)
        obs.extend(more_obs)
    a.out.mkdir(parents=True, exist_ok=True)
    ordered = [rows[i] for i in sorted(rows)]
    write_csv(a.out / "verified-datasets.csv", ordered)
    write_csv(a.out / "file-hashes.csv", [r for i in sorted(files) for _, r in sorted(files[i].items())])
    write_csv(a.out / "observations.csv", obs, ("Index", "kind", "count", "detail"))
    totals = {k: sum(int(r[k]) for r in ordered) for k in COUNTS}
    summary = dict(mode="pilot-only dry check" if a.pilot_only else "full cohort", datasets=len(rows),
                   original_pilots=len(pilot_ids), new_datasets=len(rows)-len(pilot_ids), overlapping_ids=0,
                   verified_record_files=sum(len(v) for v in files.values()), counts=totals,
                   max_abs_component_sum_error=max(float(r["max_abs_component_sum_error"]) for r in ordered),
                   max_abs_recomputed_p_difference=max(float(r["max_recomputed_p_diff"]) for r in ordered),
                   full_transferred_file_recheck=bool(a.records_root and not a.pilot_only))
    outputs = ("verified-datasets.csv", "file-hashes.csv", "observations.csv")
    summary["output_sha256"] = {name: fingerprint(a.out / name) for name in outputs}
    (a.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
