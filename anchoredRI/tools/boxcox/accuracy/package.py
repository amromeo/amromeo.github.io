"""Publish benchmark tables as readable CSVs. Split large supplements at record boundaries."""
import csv
import io
import json
import shutil
import sys
from pathlib import Path

root = Path(sys.argv[1])
tables = root/"tables"
here = Path(__file__).resolve().parent
dest = here/"results"
assert json.loads((tables/"completion.json").read_text())["fitted"] == 3780
dest.mkdir(exist_ok=True)
shutil.copyfile(tables/"anchoredRI-comparator.csv",here.parent/"anchorRefine-comparator/anchoredRI-comparator.csv")
for p in tables.iterdir():
    if p.name in ("stage-comparator.csv","anchoredRI-comparator.csv") or not p.is_file():
        continue
    if p.suffix != ".csv" or p.stat().st_size < 510000:
        shutil.copyfile(p,dest/p.name)
        continue
    folder=dest/p.stem
    folder.mkdir(exist_ok=True)
    with p.open(newline="") as f:
        source=csv.reader(f); header=next(source); target=None
        try:
            for row in source:
                if target is None or target.tell()>475000:
                    if target is not None: target.close()
                    part=1 if target is None else part+1
                    target=(folder/f"part-{part:03d}.csv").open("w",newline="")
                    writer=csv.writer(target,lineterminator="\n");writer.writerow(header)
                writer.writerow(row)
        finally:
            if target is not None:target.close()
