"""Join Actions job start/finish times to in-process CPU measurements."""
import csv
from datetime import datetime
import json
from pathlib import Path
import sys

api, received, dest = map(Path, sys.argv[1:4])
timers = {}
for path in received.rglob("*.tsv"):
    if path.parent.name != "account":
        continue
    if path.stem in timers:
        raise ValueError(f"Duplicate CPU receipt {path.stem}")
    with path.open() as stream:
        timers[path.stem] = sum(float(user) + float(system) for _, _, user, system in csv.reader(stream, delimiter="\t"))

with dest.open("w", newline="") as output:
    writer = csv.writer(output)
    writer.writerow(["role", "job", "indices", "wall_seconds", "cpu_seconds"])
    jobs = json.loads(api.read_text())["jobs"]
    for job in jobs:
        if job["name"] == "report-single" and not job["completed_at"]:
            continue
        if job["conclusion"] != "success":
            raise ValueError(f'Incomplete job {job["name"]}, {job["conclusion"]}')
        delta = datetime.fromisoformat(job["completed_at"]) - datetime.fromisoformat(job["started_at"])
        name = job["name"]
        role, key = name.split("-", 1)
        if name not in timers:
            raise ValueError(f'Missing CPU receipt {name}')
        writer.writerow([role, job["id"], key if role == "pilot" else "", delta.total_seconds(), timers[name]])
