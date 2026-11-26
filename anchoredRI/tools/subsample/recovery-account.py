"""Reconcile successful Task H pilot measurements and sunk attempts across runs."""
import argparse
import csv
from datetime import datetime
import json
import math
from pathlib import Path


def read_jobs(path):
    with path.open() as stream:
        pages = json.load(stream)
    if isinstance(pages, dict):
        pages = [pages]
    jobs = [job for page in pages for job in page["jobs"]]
    if not jobs or len({job["id"] for job in jobs}) != len(jobs):
        raise ValueError(f"Empty or duplicate jobs in {path}")
    if any(page["total_count"] != len(jobs) for page in pages if "total_count" in page):
        raise ValueError(f"Incomplete job pagination for {path}")
    return jobs


def read_receipts(folder):
    receipts = {}
    for path in sorted(folder.rglob("account/*.tsv")):
        with path.open(newline="") as stream:
            rows = list(csv.reader(stream, delimiter="\t"))
        if not rows or (path.stem in receipts and receipts[path.stem] != rows):
            raise ValueError(f"Empty or conflicting CPU receipts for {path.stem} in {folder}")
        if any(len(row) != 4 for row in rows) or len({row[0] for row in rows}) != len(rows):
            raise ValueError(f"Malformed timing phases in {path}")
        for phase, wall, user, system in rows:
            if not phase or any(not math.isfinite(float(t)) or float(t) < 0 for t in (wall, user, system)):
                raise ValueError(f"Invalid phase accounting in {path}")
        receipts[path.stem] = rows
    return receipts


def write(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def summarize(runs, roster):
    pilot_ids = set(roster)
    selected = set()
    all_ids = set()
    billable, ledger = [], []
    fixed_builds = 0
    archive_jobs = set()
    for run_id, api, artifacts in runs:
        receipts = read_receipts(artifacts)
        jobs = read_jobs(api)
        if set(receipts) - {job["name"] for job in jobs}:
            raise ValueError(f"Accounting receipts from another workflow in {artifacts}")
        for job in jobs:
            name = job["name"]
            if str(job["run_id"]) != run_id or job["id"] in all_ids:
                raise ValueError(f"Run or job ID mismatch for {name}")
            all_ids.add(job["id"])
            role, _, key = name.partition("-")
            result = job["conclusion"]
            skipped = result == "skipped"
            active = result is None and name == "report-single" and run_id == runs[-1][0]
            if result is None and not active:
                raise ValueError(f"Unfinished dependency {run_id}/{name}")
            row = dict(run_id=run_id, attempt=job.get("run_attempt", ""), job_id=job["id"], name=name,
                       conclusion=result or "in progress", started_at=job["started_at"] or "",
                       completed_at=job["completed_at"] or "", wall_seconds="", measured_cpu_seconds="",
                       cpu_coverage="", phases="", projection_role="", url=job["html_url"])
            if skipped or active:
                row["cpu_coverage"] = "skipped" if skipped else "current reporting job still running"
                ledger.append(row)
                continue
            if not row["started_at"] or not row["completed_at"]:
                raise ValueError(f"Missing Actions runner timestamps for {run_id}/{name}")
            wall = (datetime.fromisoformat(row["completed_at"]) - datetime.fromisoformat(row["started_at"])).total_seconds()
            if wall < 0:
                raise ValueError(f"Negative runner time for {name}")
            phases = receipts.get(name, [])
            if result == "success" and not phases:
                raise ValueError(f"Missing successful CPU receipt for {run_id}/{name}")
            cpu = sum(float(user) + float(system) for _, _, user, system in phases)
            phase_names = [r[0] for r in phases]
            projection_role = role
            selected_idx = ""
            if result == "success":
                coverage = "completed timed processes"
                if role == "pilot":
                    if key not in pilot_ids or key in selected or "pilot" not in phase_names:
                        raise ValueError(f"Duplicate or unexpected successful pilot {name}")
                    matches = list(artifacts.glob(f"**/records/{key}/COMPLETE"))
                    if not matches:
                        raise ValueError(f"No complete pilot record for {run_id}/{name}")
                    selected.add(key)
                    selected_idx = key
                if role == "build":
                    fixed_builds += 1
                if role == "archive":
                    if key in archive_jobs:
                        raise ValueError(f"Repeated archive fitting {name}")
                    archive_jobs.add(key)
                if name == "recovery-preflight":
                    projection_role = "recovery preflight"
            else:
                coverage = "partial CPU lower bound" if phases else "CPU unmeasured"
                projection_role = f"failed {role} (CPU lower bound)"
            row.update(wall_seconds=wall, measured_cpu_seconds=round(cpu, 6), phases=" | ".join(phase_names),
                       cpu_coverage=coverage, projection_role=projection_role)
            ledger.append(row)
            billable.append(dict(role=projection_role, job=job["id"], indices=selected_idx,
                                 wall_seconds=wall, cpu_seconds=round(cpu, 6)))
    if selected != pilot_ids or len(selected) != 20:
        raise ValueError(f"Pilot cohort incomplete. Missing {sorted(pilot_ids-selected)}")
    if fixed_builds != 1 or archive_jobs != {str(i) for i in range(1, 49)}:
        raise ValueError("Expected one successful original build and precisely 48 original archive jobs")
    return billable, ledger


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run", action="append", required=True, nargs=3, metavar=("ID", "JOBS_API_JSON", "ARTIFACTS"))
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--roster", type=Path, default=Path(__file__).with_name("pilot-roster.csv"))
    args = p.parse_args()
    runs = [(r[0], Path(r[1]), Path(r[2])) for r in args.run]
    if len(runs) != 3 or len({r[0] for r in runs}) != len(runs):
        p.error("Supply the aborted build, original pilot and recovery runs in order")
    with args.roster.open(newline="") as stream:
        roster = [row["Index"] for row in csv.DictReader(stream)]
    billable, ledger = summarize(runs, roster)
    args.out.mkdir(parents=True, exist_ok=True)
    write(args.out / "jobs.csv", billable)
    write(args.out / "job-ledger.csv", ledger)
    failed = [r for r in ledger if r["projection_role"].startswith("failed")]
    entries = ["# Compute used by interrupted jobs", "",
               "Runner time comes from GitHub’s job start and completion timestamps. CPU is the sum of the "
               "user and system seconds recorded by completed timed commands. A failed command can exit before "
               "its timing receipt is written. Its unrecorded CPU time cannot be recovered from the runner duration.",
               "", "| Run | Job | Occupied runner seconds | Recorded CPU seconds, lower bound | Completed timing receipts |",
               "|---|---|---:|---:|---|"]
    for r in failed:
        entries.append(f"| {r['run_id']} | {r['name']} | {r['wall_seconds']:.1f} | "
                       f"{r['measured_cpu_seconds']:.3f} | {r['phases'] or 'None'} |")
    entries += ["", f"The {len(failed)} interrupted jobs occupied runners for "
                f"{sum(r['wall_seconds'] for r in failed)/3600:.4f} hours. They have "
                f"{sum(r['measured_cpu_seconds'] for r in failed)/3600:.4f} recorded CPU-hours plus unmeasured CPU. "
                "Total study CPU is therefore a lower bound. Remaining costs use the 20 completed pilot records. "
                "Failed jobs enter spent cost once and do not affect per-stratum fitting times or pilot setup averages.",
                "", "An in-progress reporting job will enter spent cost after it completes. "
                "[Job ledger](job-ledger.csv) · [Projection input](jobs.csv)", ""]
    (args.out / "ACCOUNTING.md").write_text("\n".join(entries))


if __name__ == "__main__":
    main()
