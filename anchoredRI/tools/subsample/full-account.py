"""Summarize observed study process times and GitHub runner occupancy."""
import argparse
from collections import defaultdict
import csv
from datetime import datetime
import json
from pathlib import Path
import statistics

STUDY = Path(__file__).resolve().parent


def read(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def write(path, rows):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("records", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("job_pages", nargs="+", type=Path, help="One paginated GitHub jobs JSON per attempt")
    cfg = ap.parse_args()
    cfg.out.mkdir(parents=True, exist_ok=True)
    roster = {r["Index"] for r in read(STUDY / "full-roster.csv")}
    pilot = {r["Index"] for r in read(STUDY / "pilot-roster.csv")}
    if len(roster | pilot) != 1440 or roster & pilot:
        raise ValueError("Cohort mismatch")
    phases = defaultdict(list)
    for idx in sorted(roster | pilot, key=int):
        data = read(cfg.records / idx / "timing.csv")
        totals = [r for r in data if r["level"] == "total" and r["phase"] == "dataset"]
        if len(totals) != 1 or {r["Index"] for r in data} != {idx}:
            raise ValueError(f"Missing dataset timing at {idx}")
        subtotals = defaultdict(lambda: [0., 0.])
        for r in data:
            for j, field in enumerate(("wall_seconds", "cpu_seconds")):
                subtotals[r["level"], r["phase"]][j] += float(r[field])
        for key, sums in subtotals.items():
            partition = "original pilot" if idx in pilot else "full-run jobs"
            phases[partition, *key].append(sums)
            phases["complete cohort", *key].append(sums)
    table = []
    for (group, level, phase), costs in sorted(phases.items()):
        data = dict(cohort=group, level=level, phase=phase, datasets_observed=len(costs))
        for field, values in zip(("wall", "cpu"), zip(*costs)):
            data.update({f"{field}_hours": sum(values) / 3600,
                         f"{field}_median_seconds": statistics.median(values),
                         f"{field}_min_seconds":min(values), f"{field}_max_seconds":max(values)})
        table.append(data)
    write(cfg.out / "phase-times.csv", table)
    jobs = {}
    for path in cfg.job_pages:
        payload = json.loads(path.read_text())
        for page in ([payload] if isinstance(payload, dict) else payload):
            for job in page["jobs"]:
                if job["id"] in jobs:
                    raise ValueError(f"Duplicate job {job['id']}")
                start = datetime.fromisoformat(job["started_at"])
                end = datetime.fromisoformat(job["completed_at"])
                jobs[job["id"]] = dict(id=job["id"], run=job["run_id"], attempt=job["run_attempt"],
                    job=job["name"], conclusion=job["conclusion"], start=start.isoformat(), end=end.isoformat(),
                    occupied_seconds=max(0, (end-start).total_seconds()) if job["conclusion"] != "skipped" else 0)
    write(cfg.out / "actions-jobs.csv", list(jobs.values()))
    observed = [j for j in jobs.values() if j["conclusion"] != "skipped"]
    elapsed = (max(datetime.fromisoformat(j["end"]) for j in observed) -
               min(datetime.fromisoformat(j["start"]) for j in observed)).total_seconds()
    successful = {j["job"] for j in jobs.values() if j["conclusion"] == "success"}
    if len([n for n in successful if n.startswith("Fit block ")]) != 178:
        raise ValueError("The completed workflow must contain all 178 fit blocks")
    summary = dict(full_run_occupied_runner_hours=sum(j["occupied_seconds"] for j in observed)/3600,
                   first_job_start_to_last_job_end_hours=elapsed/3600,
                   full_run_dataset_process_cpu_hours=sum(r[1] for r in phases["full-run jobs","total","dataset"])/3600,
                   full_run_dataset_wall_hours=sum(r[0] for r in phases["full-run jobs","total","dataset"])/3600,
                   complete_cohort_dataset_process_cpu_hours=sum(r[1] for r in phases["complete cohort","total","dataset"])/3600,
                   workflow_attempts=len({j["attempt"] for j in observed}))
    (cfg.out / "actual-compute.json").write_text(json.dumps(summary, indent=2)+"\n")


if __name__ == "__main__":
    main()
