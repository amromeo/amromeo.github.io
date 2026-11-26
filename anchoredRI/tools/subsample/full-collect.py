"""Read checked release assets into a directory for full-cohort scoring."""
import argparse
import csv
import hashlib
import io
from pathlib import Path
import shutil
import tarfile

STUDY = Path(__file__).resolve().parent
PILOT = STUDY / "results/pilot-45854822638"
PILOT_FILES = {"COMPLETE", "decisions.csv", "draws.csv", "fits.csv", "input-id.rds",
               "medians.csv", "models.rds", "positions.rds", "timing.csv"}


def require(test, message):
    if not test:
        raise ValueError(message)


def sha(path):
    with path.open("rb") as source:
        return hashlib.file_digest(source, "sha256").hexdigest()


def rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def copy_checked(archive, member, target, digest=None):
    require(member.isfile(), f"Expected a regular file: {member.name}")
    target.parent.mkdir(parents=True, exist_ok=True)
    tmp = target.with_suffix(target.suffix + ".tmp")
    with archive.extractfile(member) as stream, tmp.open("wb") as output:
        shutil.copyfileobj(stream, output)
    if digest:
        require(sha(tmp) == digest, f"Input digest changed: {member.name}")
    if target.exists():
        require(sha(target) == sha(tmp), f"Existing output differs: {target}")
        tmp.unlink()
    else:
        tmp.replace(target)


def verify_group(asset, tag, ids, out):
    stated, filename = asset.with_suffix(".tar.sha256").read_text().strip().split()
    require(filename == asset.name and sha(asset) == stated, f"Bad outer checksum: {asset}")
    with tarfile.open(asset, "r:") as archive:
        names = archive.getnames()
        require(len(names) == len(set(names)), f"Repeated archive names in {asset}")
        with archive.extractfile("SHA256SUMS") as member:
            hashes = {name: value for value, name in (row.split("  ", 1) for row in member.read().decode().splitlines())}
        require(set(names) == set(hashes) | {"SHA256SUMS"}, f"Unverified entries in {asset}")
        present = set()
        for name, digest in hashes.items():
            parts = Path(name).parts
            if tag == "pilot":
                require(len(parts) == 3 and parts[:2] == ("pilot", "data") and parts[2].endswith(".tar.gz"),
                        f"Unexpected pilot file {name}")
                idx = int(parts[2].removesuffix(".tar.gz"))
                require(idx in ids, f"Unexpected pilot index {idx}")
                target = out / name
                copy_checked(archive, archive.getmember(name), target, digest)
                with tarfile.open(target) as pilot:
                    require(len(pilot.getnames()) == len(PILOT_FILES) and set(pilot.getnames()) == PILOT_FILES,
                            f"Unexpected files in pilot dataset {idx}")
                    for f in PILOT_FILES:
                        copy_checked(pilot, pilot.getmember(f), out / "records" / str(idx) / f)
            else:
                require(len(parts) == 3 and parts[0] in {"records", "verification"}
                        and parts[1].isdigit() and parts[2] not in (".", ".."), f"Unexpected release path {name}")
                idx = int(parts[1])
                require(idx in ids, f"Out-of-group dataset {idx}")
                copy_checked(archive, archive.getmember(name), out / name, digest)
            present.add(idx)
        require(present == ids, f"Incomplete release bundle {tag}")
    return asset.name, len(ids), asset.stat().st_size, stated


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("assets", type=Path)
    p.add_argument("output", type=Path)
    opt = p.parse_args()
    pilots = {int(r["Index"]) for r in rows(STUDY / "pilot-roster.csv")}
    groups = {"pilot": pilots, **{f"{n:02}": set() for n in range(1, 9)}}
    for r in rows(STUDY / "full-roster.csv"):
        groups[f"{1 + (int(r['block']) - 1) // 23:02}"].add(int(r["Index"]))
    seen = set()
    ledger = []
    for tag, ids in groups.items():
        require(not ids & seen, "Duplicate datasets across bundles")
        archive = opt.assets / f"task-h-records-{tag}.tar"
        ledger.append(verify_group(archive, tag, ids, opt.output))
        seen.update(ids)
    require(len(seen) == 1440 and {int(p.name) for p in (opt.output / "records").iterdir()} == seen,
            "Expected exactly 1,440 datasets")
    with (opt.output / "release-assets.csv").open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["release_asset", "datasets", "bytes", "sha256"])
        writer.writerows(ledger)


if __name__ == "__main__":
    main()
