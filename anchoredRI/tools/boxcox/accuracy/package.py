"""Publish benchmark tables as readable CSVs. Split large supplements at record boundaries."""
import csv
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

def replace_opening(path, opening):
    original=path.read_text()
    start=original.index("\n\n")+2
    stop=original.index("\n\n",start)
    path.write_text(original[:start]+opening+original[stop:])


opening=("The small-N RIbench and simulation-grid accuracy comparison is complete for the improved baseline "
         "and Arm B. Arm A stopped during reliability testing. N = 50,000 awaits the separate runtime review. "
         "[Accuracy results](accuracy/results/REPORT.md) · [Protocol](PROTOCOL.md) · "
         "[Baseline comparator for anchorRefine](anchorRefine-comparator/README.md)")
for file in ("README.md","STATUS.md"):
    replace_opening(here.parent/file,opening)
readme=here.parent/"anchorRefine-comparator/README.md"
original=readme.read_text()
a=original.index("| Coverage");b=original.index("## Pairing datasets",a)
readme.write_text(original[:a]+"""| Coverage at the small-N and grid review | Available | Planned |
| --- | ---: | ---: |
| RIbench N ≤ 5,000 | 2,880 | 2,880 |
| RIbench N = 50,000 | 144 previously fitted Hb datasets | 1,440 |
| Simulation grid, replicates 1 to 20 | 900 | 900 |

The file contains 3,924 unique datasets. All 396 previously exported small-N rows retained the same fits, limits and source hashes. The 144 previously exported large-N rows remain available. The full large-N comparison awaits approval. Grid rows use keys of the form `s-CC-RR` (case and replicate), the original routine seed and the family-specific screened seed. The design deliberately reuses each screened sample across scenarios in a family and replicate.

"""+original[b:])
s=readme.read_text().replace("The case and replicate fields are reserved for the simulation grid.","The case and replicate fields identify simulation grid cells and replicates.").replace("Pathological fraction is the RIbench design fraction.","Pathological fraction is the design fraction in each benchmark.").replace("[Export script](../post-round3/export-saved-comparator.R)","[Stage export script](../accuracy/extract.R) · [Accuracy report](../accuracy/results/REPORT.md)")
readme.write_text(s)
