"""Finalize a completed Task H report artifact without running R or fitting.

Usage: python3 finalize-report.py PILOT_REPORT_ARTIFACT COMPLETED_JOBS_JSON OUTPUT_DIR
The artifact directory is the one containing jobs.csv, account/, and pilot-report/.
"""
import argparse
import csv
import hashlib
import json
import math
import re
import shutil
from datetime import datetime
from pathlib import Path


def read(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write(path, rows):
    with path.open("w", newline="") as stream:
        w = csv.DictWriter(stream, list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def api_jobs(path):
    pages = json.loads(path.read_text())
    if isinstance(pages, dict):
        return pages["jobs"]
    return [j for page in pages for j in page["jobs"]]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("artifact", type=Path)
    p.add_argument("jobs_json", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument("--design", type=Path, default=Path("inputs/testsets.csv"))
    p.add_argument("--roster", type=Path, default=Path("anchoredRI/tools/subsample/pilot-roster.csv"))
    args = p.parse_args()
    source = args.artifact / "pilot-report"
    if args.out.exists() and any(args.out.iterdir()):
        p.error("Output directory must be empty.")
    jobs = api_jobs(args.jobs_json)
    report_jobs = [j for j in jobs if j["name"] == "report-single"]
    if len(report_jobs) != 1:
        raise ValueError("Expected the completed successful attempt's single report job")
    job = report_jobs[0]
    if job["conclusion"] != "success" or not job["completed_at"] or not job["started_at"]:
        raise ValueError("The report job has not completed successfully")
    wall = (datetime.fromisoformat(job["completed_at"]) - datetime.fromisoformat(job["started_at"])).total_seconds()
    with (args.artifact / "account" / "report-single.tsv").open() as stream:
        receipt = list(csv.reader(stream, delimiter="\t"))
    names = [row[0] for row in receipt]
    if not {"pilot-report", "score", "prepare"}.issubset(names) or len(names) != len(set(names)):
        raise ValueError(f"Incomplete or repeated report CPU receipts: {names}")
    cpu = sum(float(user) + float(sys) for _, _, user, sys in receipt)
    if wall <= 0 or cpu <= 0:
        raise ValueError("Report timings must be positive")

    recovery = (args.artifact / "accounting/jobs.csv").is_file()
    earlier = read(args.artifact / ("accounting/jobs.csv" if recovery else "jobs.csv"))
    if any(r["role"] == "report" for r in earlier):
        raise ValueError("The pilot artifact already includes the complete report accounting")
    previous = {r["job"] for r in earlier}
    if recovery:
        ledger = read(args.artifact / "accounting/job-ledger.csv")
        previous &= {r["job_id"] for r in ledger if r["run_id"] == str(job["run_id"])}
    if {str(j["id"]) for j in jobs if j is not job and j["conclusion"] != "skipped"} != previous:
        raise ValueError("The completed-at data must come from the same workflow run attempt as the artifact")

    # Count and preserve the original grouped RIbench 1.0.2 scores byte-for-byte.
    scored = read(source / "pilot-scores.csv")
    estimates = read(source / "pilot-endpoints.csv")
    key = lambda r: (r["Index"], r["method"])
    if len(scored) != 240 or len({key(r) for r in scored}) != 240 or {key(r) for r in estimates} != {key(r) for r in scored}:
        raise ValueError("Expected one score per estimator for 12 methods on 20 original IDs")
    native_sha = hashlib.sha256((source / "pilot-scores.csv").read_bytes()).hexdigest()
    rows = read(source / "compute-projection.csv")
    future = [r for r in rows if r["item"] == "Projected final reporting"]
    if len(future) != 1:
        raise ValueError("Missing preliminary final-report allowance")
    allowance = future[0].copy()
    # In the prospective main run the report is dispatched only after every fit
    # finishes. A temporary recovery job may occupy a runner waiting on the
    # original workflow; retain that cost in spent time only.
    wait_seconds = sum((datetime.fromisoformat(step["completed_at"]) -
                        datetime.fromisoformat(step["started_at"])).total_seconds()
                       for step in job.get("steps", [])
                       if step["name"] == "Wait for the other original pilots to finish" and step["conclusion"] == "success")
    report_allowance = wall - wait_seconds
    if report_allowance <= 0:
        raise ValueError("Invalid report wait-step timing")
    future[0].update(wall_seconds=report_allowance, cpu_seconds=cpu)
    rows.insert(next(i for i, r in enumerate(rows) if r["item"].startswith("Projected")),
                dict(item="Measured report", jobs=1, wall_seconds=wall, cpu_seconds=cpu))
    for r in rows:
        for field in ("wall_seconds", "cpu_seconds"):
            r[field] = float(r[field])
            if not math.isfinite(r[field]) or r[field] < 0:
                raise ValueError("Invalid projection entry")

    report = (source / "REPORT.md").read_text()
    section = re.search(r"(## Compute projection\n\n)(.*?)(\n## Pilot output checks)", report, re.S)
    if not section:
        raise ValueError("Missing expected preliminary report sections")
    body = section.group(2)
    tailpos = body.index("\nEstimated remaining elapsed time")
    tail = body[tailpos:]
    # Reproduce the frozen Index-order assignment and include the measured report cost.
    # The report job has completed; none of the 1,420 future datasets are fitted here.
    roster = read(args.roster)
    clock = read(source / "job-overhead.csv")
    index_walls = {r["indices"]: float(r["index_wall"]) for r in clock}
    strata = {(r["Analyte"], r["pathology_stratum"]): r["Index"] for r in roster}
    if len(strata) != 20 or set(strata.values()) != set(index_walls) or len(clock) != 20:
        raise ValueError("Pilot timing and strata do not reconcile")
    mean_overhead = sum(float(r["overhead_wall"]) for r in clock) / len(clock)
    slots = [0.0] * 20
    remaining_count = 0
    for d in sorted(read(args.design), key=lambda r: int(r["Index"])):
        if int(d["N"]) != 50000 or d["Index"] in index_walls:
            continue
        s = d["Analyte"], "<=20%" if float(d["fractionPathol"]) <= .2 else ">20%"
        first = min(range(len(slots)), key=slots.__getitem__)
        slots[first] += mean_overhead + index_walls[strata[s]]
        remaining_count += 1
    if remaining_count != 1420:
        raise ValueError("Incomplete future design")
    time_text = re.search(r"Estimated remaining elapsed time is \*\*([0-9.]+) hours\*\*", tail)
    if not time_text:
        raise ValueError("Preliminary elapsed estimate unavailable")
    prior = (max(slots) + float(allowance["wall_seconds"])) / 3600
    if f"{prior:.2f}" != time_text.group(1):
        raise ValueError("Independent readback disagrees with preliminary elapsed duration")
    adjusted = (max(slots) + report_allowance) / 3600
    tail = tail[:time_text.start()] + f"Estimated remaining elapsed time is approximately **{adjusted:.2f} hours**" + tail[time_text.end():]
    lines = ["| Work | Jobs | Occupied runner-hours | CPU-hours |", "|---|---:|---:|---:|"]
    for r in rows:
        lines.append(f"| {r['item']} | {r['jobs']} | {r['wall_seconds']/3600:.3f} | {r['cpu_seconds']/3600:.3f} |")
    for title, subset in [("Remaining", [r for r in rows if r["item"].startswith("Projected")]), ("Total", rows)]:
        lines.append(f"| **{title}** | | **{sum(r['wall_seconds'] for r in subset)/3600:.3f}** | "
                     f"**{sum(r['cpu_seconds'] for r in subset)/3600:.3f}** |")
    lines += ["", f"The finished pilot-report job occupied its runner for {wall:.1f} seconds and used "
              f"{cpu:.2f} measured CPU seconds. It is included once in the spent balance. The future report "
              f"allowance uses {report_allowance:.1f} runner seconds and {cpu:.2f} measured CPU seconds. "
              f"The {wait_seconds:.1f} seconds spent waiting for the original pilot run appear only in spent time. "
              f"This replaces the preliminary allowance "
              f"of {float(allowance['wall_seconds']):g} runner seconds and {float(allowance['cpu_seconds']):g} CPU seconds. "
              "The full cohort will require larger transfers and scoring tables. The remaining report allowance "
              "assumes equivalent runtime and is an estimate.", tail.strip(), ""]
    refreshed = report[:section.start(2)] + "\n".join(lines) + report[section.end(2):]
    shutil.copytree(source, args.out, dirs_exist_ok=True)
    if recovery:
        shutil.copy2(args.artifact / "accounting/ACCOUNTING.md", args.out / "ACCOUNTING.md")
        old_report = [r for r in ledger if r["job_id"] == str(job["id"])]
        if len(old_report) != 1:
            raise ValueError("Missing in-progress recovery report from the job ledger")
        old_report[0].update(conclusion="success", started_at=job["started_at"], completed_at=job["completed_at"],
                             wall_seconds=wall, measured_cpu_seconds=cpu, cpu_coverage="completed timed processes",
                             phases=" | ".join(names), projection_role="report")
        write(args.out / "job-ledger.csv", ledger)
        refreshed = refreshed.replace("## Compute projection\n", "## Compute projection\n\n"
                    "The total CPU column is a lower bound. Seven interrupted jobs used additional unmeasured CPU. "
                    "The remaining projection uses the 20 completed pilot jobs. "
                    "[Interrupted job accounting](ACCOUNTING.md)\n", 1)
        accounting = (args.out / "ACCOUNTING.md").read_text().replace(
            "An in-progress reporting job will enter spent cost after it completes. ", "")
        (args.out / "ACCOUNTING.md").write_text(accounting.replace("](jobs.csv)", "](jobs-final.csv)"))
    (args.out / "REPORT.md").write_text(refreshed)
    write(args.out / "compute-projection.csv", rows)
    write(args.out / "jobs-final.csv", earlier + [dict(role="report", job=job["id"], indices="",
                                                       wall_seconds=wall, cpu_seconds=cpu)])
    shutil.copy2(args.artifact / "account" / "report-single.tsv", args.out / "report-cpu.tsv")
    (args.out / "report-timing.json").write_text(json.dumps(
        dict(actions_job=job["html_url"], run_id=job["run_id"], started_at=job["started_at"],
             completed_at=job["completed_at"], wall_seconds=wall, cpu_seconds=cpu, recovery_polling_wall_seconds=wait_seconds,
             preliminary_final_report_allowance=allowance, pilot_scores_sha256=native_sha), indent=2) + "\n")
    assert hashlib.sha256((args.out / "pilot-scores.csv").read_bytes()).hexdigest() == native_sha
    print(f"{args.out}/REPORT.md  |  added measured report {wall:.0f} runner seconds, {cpu:.2f} CPU seconds")


if __name__ == "__main__":
    main()
