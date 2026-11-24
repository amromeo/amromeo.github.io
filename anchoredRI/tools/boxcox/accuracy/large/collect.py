"""Choose the latest complete checkpoint for each shard, including after an Actions retry."""
import csv
import os
import re
import shutil
import sys
from pathlib import Path

root = Path(sys.argv[1])
folders = {}
for path in (root / "downloads").glob("task-e-50000-fits-*-*"):
    match = re.fullmatch(r"task-e-50000-fits-(\d+)-(\d+)", path.name)
    if match and len(list(path.glob("r-*.rds"))) == 60:
        shard, attempt = map(int, match.groups())
        if shard not in folders or attempt > folders[shard][0]:
            folders[shard] = (attempt, path)
assert set(folders) == set(range(1, 25)), "Expected 24 complete fitting batches"
fits = root / "fits"
fits.mkdir(exist_ok=True)
provenance = root / "tables" / "provenance"
provenance.mkdir(parents=True, exist_ok=True)
with (provenance / "batches.csv").open("w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["run", "shard", "attempt", "datasets", "image"])
    for shard, (attempt, folder) in sorted(folders.items()):
        for src in folder.glob("r-*.rds"):
            dest = fits / src.name
            assert not dest.exists(), "Duplicate checkpoint " + src.name
            os.link(src, dest)
        image = (folder / "runtime.txt").read_text().splitlines()[0]
        writer.writerow([os.environ["GITHUB_RUN_ID"], shard, attempt, 60, image])
        for name in ("runtime.txt", "native-source-hashes.csv", "fitter-source-hashes.csv"):
            dest = provenance / (image.replace(" ", "-") + "-" + name)
            if dest.exists():
                assert dest.read_bytes() == (folder / name).read_bytes()
            else:
                shutil.copyfile(folder / name, dest)
assert len(list(fits.glob("r-*.rds"))) == 1440
