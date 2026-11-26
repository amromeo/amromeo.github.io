"""Summarize verified timing-pilot records; never invokes a fitter."""
import argparse
import csv
import math
import statistics
import subprocess
from collections import Counter, defaultdict
from pathlib import Path


def read(path):
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write(path, rows, fields=None):
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        out = csv.DictWriter(handle, fields or list(rows[0]))
        out.writeheader()
        out.writerows(rows)


def number(value):
    return float(value) if value not in ("", "NA", "NaN", None) else math.nan


def fmt(value, places=3):
    return f"{value:.{places}f}" if math.isfinite(value) else "—"


def finite_mean(values):
    values = [v for v in values if math.isfinite(v)]
    return statistics.mean(values) if values else math.nan


def paired_score(records, key):
    return [number(r[key]) for r in records]


def by_key(rows, *fields):
    result = {}
    for row in rows:
        key = tuple(row[f] for f in fields)
        if key in result:
            raise ValueError(f"Duplicate {fields} = {key}")
        result[key] = row
    return result


def endpoints(idx, method, row, raw=False):
    prefix = "raw_" if raw else ""
    return dict(Index=idx, method=method, lower=row[prefix + "lower"], upper=row[prefix + "upper"])


def validate_draws(fits, medians, draws, roster):
    position_hash = {}
    for draw in draws:
        bank, b = draw["seed_set"], int(draw["draw"])
        seed = (100000000 if bank == "A" else 200000000) + 100 * int(roster["Index"]) + b
        if int(draw["sampling_seed"]) != seed or int(draw["positions_n"]) != 5000:
            raise ValueError(f"Invalid draw {bank}/{b}")
        position_hash[bank, str(b)] = draw["position_hash"]
    if len(draws) != 100 or len(position_hash) != 100 or len(set(position_hash.values())) != 100:
        raise ValueError(f"Missing or repeated sample positions for Index {roster['Index']}")
    anomalies = []
    for r in fits:
        if r["routine_hash"] != roster["routine_hash"] or r["anchor_hash"] != roster["anchor_hash"]:
            raise ValueError(f"Routine/anchor mismatch for Index {roster['Index']}")
        if r["first_status"] not in ("reported", "withheld", "declined"):
            raise ValueError(f"Invalid first-stage status for {roster['Index']}")
        if r["fit_type"] == "subsample" and r["position_hash"] != position_hash[r["seed_set"], r["draw"]]:
            raise ValueError(f"Position digest mismatch for {roster['Index']}")
        expected_n = int(roster["routine_n"]) if r["fit_type"] != "subsample" else 5000
        expected_bic_n = expected_n + 40 if r["fit_type"] != "W" else 5040
        if number(r["routine_n"]) != expected_n or not math.isclose(number(r["bic_n"]), expected_bic_n,
                                                                     rel_tol=0, abs_tol=1e-8):
            raise ValueError(f"Unexpected sample size or BIC weight in {r}")
        if r["first_status"] == "withheld":
            anomalies.append(("technical first stage", r["cap"], r["fit_type"], r.get("seed_set"),
                              r.get("draw"), r["first_reason"]))
        elif not (0 <= number(r["share"]) <= 1 and 1 <= number(r["components"]) <= int(r["cap"])):
            raise ValueError(f"Missing or invalid completed share/component count on {r}")
    for m in medians:
        key = (m["cap"], m["seed_set"])
        rows = [r for r in fits if r["fit_type"] == "subsample" and (r["cap"], r["seed_set"]) == key]
        if sorted(int(r["draw"]) for r in rows) != list(range(1, 51)):
            raise ValueError(f"Missing/duplicate draw in {key}")
        statuses = Counter(r["first_status"] for r in rows)
        if any(statuses[status] != int(m[field]) for status, field in
               (("reported", "votes_pass"), ("declined", "votes_decline"), ("withheld", "votes_technical"))):
            raise ValueError(f"Vote discrepancy in {key}")
        keep = [r for r in rows if m["definition"] == "all" or r["first_status"] == "reported"]
        crp = roster["Analyte"] == "CRP"
        usable = [r for r in keep if math.isfinite(number(r["raw_upper"])) and
                  (crp or math.isfinite(number(r["raw_lower"])) and number(r["raw_lower"]) < number(r["raw_upper"]))]
        if int(m["contributors"]) != len(usable):
            raise ValueError(f"Median count mismatch in {key}")
        for field in (("upper",) if crp else ("lower", "upper")):
            if usable and not math.isclose(number(m[field]), statistics.median(number(r["raw_" + field]) for r in usable),
                                           rel_tol=1e-12, abs_tol=1e-12):
                raise ValueError(f"Exported median disagrees with draw records in {key}")
        if m["definition"] == "all" and len(usable) < 50:
            anomalies.append(("incomplete median", *key, m["definition"], "", f"{len(usable)}/50"))
        if m["definition"] == "passing" and not usable:
            anomalies.append(("no passing fits", *key, m["definition"], "", "0/50"))
    return [dict(Index=roster["Index"], finding=a[0], cap=a[1], fit_or_set=a[2],
                 set_or_definition=a[3], draw=a[4], detail=a[5]) for a in anomalies]


def load_pilot(records, roster):
    study = {}
    estimates = []
    times = []
    anomalies = []
    for meta in roster:
        idx = meta["Index"]
        base = records / idx
        fits, medians, decisions, timing = [read(base / f"{name}.csv")
                                           for name in ("fits", "medians", "decisions", "timing")]
        draws = read(base / "draws.csv")
        for table in (fits, medians, decisions, timing, draws):
            if not table or any(row["Index"] != idx for row in table):
                raise ValueError(f"Inconsistent Index in {base}")
        if len(fits) != 204 or len(medians) != 8 or len(decisions) != 48:
            raise ValueError(f"Incomplete fit, median or decision output in {base}")
        anomalies.extend(validate_draws(fits, medians, draws, meta))
        f = by_key([r for r in fits if r["fit_type"] != "subsample"], "fit_type", "cap")
        m = by_key(medians, "cap", "seed_set", "definition")
        ds = by_key(decisions, "cap", "seed_set", "definition", "estimator", "rule")
        full = f["full", "5"]
        estimates.append(endpoints(idx, "standard_raw", full, True))
        estimates.append(endpoints(idx, "package", full))
        for cap in ("5", "3"):
            estimates.append(endpoints(idx, "W" + cap, f["W", cap], True))
            for bank in ("A", "B"):
                for definition in ("all", "passing"):
                    key = cap, bank, definition
                    estimates.append(endpoints(idx, "median_" + "_".join(key), m[key]))
        study[idx] = (meta, full, f, m, ds)
        times.extend(timing)
    return study, estimates, times, anomalies


def first_table(study, scores, estimates, out):
    rows = []
    for idx, (meta, full, f, medians, ds) in study.items():
        first = full["first_status"]
        panel = {"reported": "A", "declined": "B", "withheld": "technical"}[first]
        for cap, bank, definition in medians:
            key = cap, bank, definition
            maj = ds[key + ("median", "majority")]["status"] == "reported"
            own = f["W", cap]["first_status"] == "reported"
            package = full["status"] == "reported"
            row = dict(Index=idx, analyte=meta["Analyte"], panel=panel, cap=cap, bank=bank,
                       definition=definition, median_released=maj, W_released=own, full_released=package)
            reference = "standard_raw" if panel == "A" else "package"
            for label, method in (("ref", reference), ("median", "median_" + "_".join(key)), ("W", "W" + cap)):
                row.update({label + "_" + k: number(scores[idx, method][k]) for k in ("z_lower", "z_upper", "d")})
                row.update({label + "_" + k: number(estimates[idx, method][k]) for k in ("lower", "upper")})
            triple = panel == "A" or (panel == "B" and maj and package)
            row["triple"] = triple and all(math.isfinite(row[l + "_d"]) for l in ("ref", "median", "W"))
            row["W_pair"] = panel in ("A", "B") and own and package and all(
                math.isfinite(row[l + "_d"]) for l in ("ref", "W"))
            rows.append(row)
    write(out / "paired-index.csv", rows)

    summaries, units = [], []
    for bank in ("A", "B"):
        for definition in ("all", "passing"):
            for cap in ("5", "3"):
                for panel in ("A", "B"):
                    group = [r for r in rows if (r["bank"], r["definition"], r["cap"], r["panel"])
                             == (bank, definition, cap, panel)]
                    triples = [r for r in group if r["triple"]]
                    wpairs = [r for r in group if r["W_pair"]]
                    s = dict(bank=bank, definition=definition, cap=cap, panel=panel, eligible=len(group),
                             majority_release=sum(r["median_released"] for r in group),
                             W_release=sum(r["W_released"] for r in group), common=len(triples),
                             released_W_on_triples=sum(r["W_released"] for r in triples),
                             majority_only=sum(r["median_released"] and not r["full_released"] for r in group),
                             package_only=sum(r["full_released"] and not r["median_released"] for r in group),
                             W_pairs=len(wpairs))
                    for method, gate in (("median", "median_released"), ("W", "W_released"), ("ref", "full_released")):
                        s[method + "_mean_d"] = finite_mean(paired_score(triples, method + "_d"))
                        for side in ("lower", "upper"):
                            s[method + "_extreme_" + side] = sum(abs(r[f"{method}_z_{side}"]) > 5 for r in triples)
                        def extreme(r):
                            return any(abs(r[f"{method}_z_{side}"]) > 5 for side in ("lower", "upper"))
                        s[method + "_extreme_any"] = sum(extreme(r) for r in triples)
                        s[method + "_all_eligible"] = finite_mean(
                            min(r[method + "_d"], 5) if r[gate] and math.isfinite(r[method + "_d"]) else 5
                            for r in group)
                        if method != "median":
                            s[method + "_extreme_W_pairs"] = sum(extreme(r) for r in wpairs)
                    for method in ("median", "W"):
                        s[method + "_paired_delta_d"] = finite_mean(r[method + "_d"] - r["ref_d"] for r in triples)
                    s["W_delta_d_own"] = finite_mean(r["W_d"] - r["ref_d"] for r in wpairs)
                    summaries.append(s)
                    for analyte in sorted({r["analyte"] for r in triples}):
                        a = [r for r in triples if r["analyte"] == analyte]
                        for side in ("lower", "upper"):
                            if side == "lower" and analyte == "CRP":
                                continue
                            units.append(dict(bank=bank, definition=definition, cap=cap, panel=panel, analyte=analyte,
                                              side=side, pairs=len(a), **{m + "_minus_ref": finite_mean(
                                                  r[m + "_" + side] - r["ref_" + side] for r in a)
                                                  for m in ("median", "W")}))
    write(out / "paired-summary.csv", summaries)
    write(out / "paired-limits-by-analyte.csv", units, ["bank", "definition", "cap", "panel",
                                                       "analyte", "side", "pairs", "median_minus_ref", "W_minus_ref"])
    return summaries


def projection(args, times, roster, design, out):
    pilot = {r["Index"]: r for r in roster}
    records = defaultdict(dict)
    phases = defaultdict(lambda: defaultdict(lambda: [0.0, 0.0]))
    for row in times:
        idx = row["Index"]
        values = [number(row[k]) for k in ("wall_seconds", "cpu_seconds")]
        if not all(math.isfinite(v) and v >= 0 for v in values):
            raise ValueError(f"Invalid phase timing {row}")
        if row["level"] == "total":
            if idx in records:
                raise ValueError(f"Duplicate total for {idx}")
            records[idx] = values
        else:
            key = row["level"], row["phase"]
            by_id = phases[key][idx]
            for i, v in enumerate(values):
                by_id[i] += v
    if set(records) != set(pilot):
        raise ValueError("Missing per-Index timing totals")
    phases["total", "All timed per-Index work"] = records
    phase_summary = []
    for (level, phase), obs in phases.items():
        # Fallback work on an Index without a fallback contributes zero.
        wall, cpu = zip(*(obs.get(idx, [0, 0]) for idx in pilot))
        phase_summary.append(dict(level=level, phase=phase, indices_recorded=len(obs),
                                  wall_total=sum(wall), wall_min=min(wall), wall_max=max(wall),
                                  cpu_total=sum(cpu), cpu_min=min(cpu), cpu_max=max(cpu)))
    phase_summary.sort(key=lambda r: ({"total": 0, "phase": 1, "detail": 2}.get(r["level"], 3), r["phase"]))
    write(out / "phase-times.csv", phase_summary)

    jobs = read(args.jobs)
    if not jobs or len({r["job"] for r in jobs}) != len(jobs):
        raise ValueError("Empty or duplicated job manifest")
    for j in jobs:
        for field in ("wall_seconds", "cpu_seconds"):
            j[field] = number(j[field])
            if not math.isfinite(j[field]) or j[field] < 0:
                raise ValueError(f"Missing measured job {field} for {j['job']}")
    pilot_jobs = [j for j in jobs if j["role"] == "pilot"]
    seen = [idx for job in pilot_jobs for idx in job["indices"].split()]
    if Counter(seen) != Counter(pilot.keys()):
        raise ValueError("Timed pilot jobs must cover each of the 20 IDs exactly once")
    overhead = []
    for j in pilot_jobs:
        ids = j["indices"].split()
        observed = [sum(records[idx][i] for idx in ids) for i in (0, 1)]
        residual = [j[f] - t for f, t in zip(("wall_seconds", "cpu_seconds"), observed)]
        # Actions timestamps have 1-second precision and R rounds process CPU time.
        if min(residual) < -2:
            raise ValueError(f"Job {j['job']} time is shorter than its aggregate per-Index time. Check parallelism and timing scope")
        overhead.append(dict(job=j["job"], indices=j["indices"], job_wall=j["wall_seconds"], job_cpu=j["cpu_seconds"],
                             index_wall=observed[0], index_cpu=observed[1], overhead_wall=max(0, residual[0]),
                             overhead_cpu=max(0, residual[1])))
    write(out / "job-overhead.csv", overhead)
    overhead_avg = [statistics.mean(r[k] for r in overhead) for k in ("overhead_wall", "overhead_cpu")]
    stratum = {(r["Analyte"], r["pathology_stratum"]): r["Index"] for r in roster}
    if len(stratum) != 20:
        raise ValueError("Expected 20 distinct pilot strata")
    remaining = []
    census = Counter()
    for r in sorted(design, key=lambda row: int(row["Index"])):
        if int(r["N"]) != 50000:
            continue
        s = (r["Analyte"], "<=20%" if number(r["fractionPathol"]) <= .2 else ">20%")
        sampled = stratum[s]
        census[sampled] += 1
        if r["Index"] not in pilot:
            remaining.append(records[sampled])
    if any(census[r["Index"]] != int(r["eligible_datasets"]) for r in roster) or len(remaining) != 1420:
        raise ValueError("Pilot stratum populations disagree with full design")
    jobs_remaining = [remaining[j:j + args.indices_per_job] for j in range(0, len(remaining), args.indices_per_job)]
    # First-available runner and index-order shards. Queue delays are excluded.
    slots = [0.] * args.concurrency
    for group in jobs_remaining:
        slot = min(range(len(slots)), key=slots.__getitem__)
        slots[slot] += sum(t[0] for t in group) + overhead_avg[0]
    fixed = []
    for role in sorted({j["role"] for j in jobs}):
        q = [j for j in jobs if j["role"] == role]
        fixed.append(dict(item="Measured " + role, jobs=len(q),
                          wall_seconds=sum(j["wall_seconds"] for j in q),
                          cpu_seconds=sum(j["cpu_seconds"] for j in q)))
    report_jobs = [j for j in jobs if j["role"] == "report"]
    future_report = [args.report_wall_seconds, args.report_cpu_seconds]
    if any(t is None for t in future_report):
        if not report_jobs:
            raise ValueError("Specify projected reporting wall and CPU time or provide a measured report job")
        future_report = [statistics.mean(j[key] for j in report_jobs) if t is None else t
                         for t, key in zip(future_report, ("wall_seconds", "cpu_seconds"))]
    if any(t < 0 for t in future_report):
        raise ValueError("Negative reporting duration")
    fixed.append(dict(item="Projected remaining per-Index work", jobs=len(jobs_remaining),
                      wall_seconds=sum(t[0] for t in remaining), cpu_seconds=sum(t[1] for t in remaining)))
    fixed.append(dict(item="Projected per-job overhead", jobs=len(jobs_remaining),
                      wall_seconds=overhead_avg[0] * len(jobs_remaining),
                      cpu_seconds=overhead_avg[1] * len(jobs_remaining)))
    fixed.append(dict(item="Projected final reporting", jobs=1, wall_seconds=future_report[0], cpu_seconds=future_report[1]))
    write(out / "compute-projection.csv", fixed)
    return phase_summary, fixed, overhead_avg, len(jobs_remaining), (max(slots) + future_report[0]) / 3600


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("records", type=Path)
    p.add_argument("inputs", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument("--roster", type=Path, default=Path(__file__).with_name("pilot-roster.csv"))
    p.add_argument("--jobs", type=Path, required=True)
    p.add_argument("--concurrency", type=int, required=True)
    p.add_argument("--indices-per-job", type=int, required=True)
    p.add_argument("--report-wall-seconds", type=float)
    p.add_argument("--report-cpu-seconds", type=float)
    args = p.parse_args()
    if min(args.concurrency, args.indices_per_job) <= 0:
        p.error("job concurrency and indices per job must be positive")
    roster = read(args.roster)
    if len(roster) != 20 or len(by_key(roster, "Index")) != 20:
        raise ValueError("Incomplete pilot roster")
    args.out.mkdir(parents=True, exist_ok=True)
    study, ep, timing, anomalies = load_pilot(args.records, roster)
    write(args.out / "pilot-endpoints.csv", ep)
    scores_path = args.out / "pilot-scores.csv"
    design_path = args.inputs / "testsets.csv"
    subprocess.run(["Rscript", str(Path(__file__).with_name("pilot-score.R")), str(args.out / "pilot-endpoints.csv"),
                    str(design_path), str(scores_path)], check=True)
    scores = by_key(read(scores_path), "Index", "method")
    if set(scores) != set(by_key(ep, "Index", "method")):
        raise ValueError("Missing scores")
    pairs = first_table(study, scores, by_key(ep, "Index", "method"), args.out)
    phases, fixed, overhead, jobs_remaining, elapsed = projection(args, timing, roster, read(design_path), args.out)
    write(args.out / "pilot-diagnostics.csv", anomalies,
          ["Index", "finding", "cap", "fit_or_set", "set_or_definition", "draw", "detail"])
    primary = [r for r in pairs if r["bank"] == "A" and r["definition"] == "all"]
    lines = ["# Task H timing pilot", "",
             "The 20 prespecified datasets span ten analytes and both pathology strata. Each used both banks of 50 draws. "
             "These paired scores use the pilot cohort for RIbench’s analyte-wide handling of zero limits. "
             "The approved full-cohort scorer requires all 1,440 datasets.", "",
             "| Full first stage | Maximum components | Eligible | Majority releases | W releases | Common triples | Mean paired error change, median | Mean paired error change, W | Extreme intervals, reference / median / W |",
             "|---|---:|---:|---:|---:|---:|---:|---:|---|"]
    for r in primary:
        extremes = " / ".join(str(r[m + "_extreme_any"]) for m in ("ref", "median", "W"))
        lines.append(f"| {'Released' if r['panel'] == 'A' else 'Declined'} | {r['cap']} | {r['eligible']} | "
                     f"{r['majority_release']} | {r['W_release']} | {r['common']} | "
                     f"{fmt(r['median_paired_delta_d'])} | {fmt(r['W_paired_delta_d'])} | {extremes} |")
    tech = sum(full["first_status"] == "withheld" for _, full, *_ in study.values())
    fall_calls = sum(r["phase"] == "fallback" and r["level"] == "detail" for r in timing)
    lines += ["", "Seed set A uses every usable draw. Panel A uses the standard full first-stage limits. "
              "Panel B requires both a released median and an actual released full-package fallback. "
              f"{tech} full first stages were technically withheld. W’s own check determines the W release count. "
              "Both paired differences use the same scorable triples. An extreme interval has at least one requested endpoint with |z| > 5.",
              "", "The paired tables include set B, passing-draw medians, all-eligible scores with nonreturns scored 5, "
              "W’s standalone releases, endpoint-specific extreme counts and changes in analyte units. "
              "[Per-Index pairs](paired-index.csv) · [Paired summaries](paired-summary.csv) · [Changes in analyte units](paired-limits-by-analyte.csv)",
              "", "## Measured phase times", "",
              "| Phase | Indices timed | Runner seconds, total | Runner seconds per Index, range | CPU seconds, total | CPU seconds per Index, range |",
              "|---|---:|---:|---:|---:|---:|"]
    for r in phases:
        name = r["phase"] if r["level"] != "detail" else r["phase"] + " (included above)"
        lines.append(f"| {name} | {r['indices_recorded']} | {fmt(r['wall_total'],1)} | "
                     f"{fmt(r['wall_min'],1)}–{fmt(r['wall_max'],1)} | {fmt(r['cpu_total'],1)} | "
                     f"{fmt(r['cpu_min'],1)}–{fmt(r['cpu_max'],1)} |")
    lines += ["", f"The standard full-data package called the fallback on {fall_calls} of the 20 datasets. "
              "Each stratum contributes a single observed time, so infrequent fallback calls, slow optimizations and differences in runner load contribute uncertainty.",
              "", "## Compute projection", "", "| Work | Jobs | Occupied runner-hours | CPU-hours |",
              "|---|---:|---:|---:|"]
    for r in fixed:
        lines.append(f"| {r['item']} | {r['jobs']} | {r['wall_seconds']/3600:.3f} | {r['cpu_seconds']/3600:.3f} |")
    for label, subset in (("Remaining", [r for r in fixed if r["item"].startswith("Projected")]), ("Total", fixed)):
        lines.append(f"| **{label}** | | **{sum(r['wall_seconds'] for r in subset)/3600:.3f}** | "
                     f"**{sum(r['cpu_seconds'] for r in subset)/3600:.3f}** |")
    lines += ["", f"Estimated remaining elapsed time is **{elapsed:.2f} hours** with {args.concurrency} concurrent runners, "
              f"{args.indices_per_job} datasets per job and {jobs_remaining} jobs in Index order, plus one final reporting job. "
              "This assumes immediate dispatch whenever a runner frees up. Queue delays and future failed-job retries would increase elapsed time.",
              "", f"Each remaining dataset receives the time of its pilot stratum (71 per stratum). "
              f"Measured mean setup, export and other job overhead adds {overhead[0]:.1f} runner seconds and "
              f"{overhead[1]:.1f} measured CPU seconds per job. "
              "The 20 completed pilot records are reused. The projection retains the timed standalone cap-3 fit "
              "in future per-Index work and counts prior archive verification and archived-seed reruns only in the measured total. "
              "Shared candidate fits for the two caps are timed once in each bank. CPU covers timed study processes. "
              "Runner time includes setup and artifact transfer. [Job overhead](job-overhead.csv) · [Projection inputs](compute-projection.csv)",
              "", "## Pilot output checks", "",
              "Readback verified draw counts, input hashes, the 50 votes per bank and cap, selected component bounds, "
              "and the coordinatewise raw medians against the exported draws."]
    issues = Counter(r["finding"] for r in anomalies)
    if issues:
        lines += ["", *[f"- {k} occurred {v} times." for k, v in issues.items()], "",
                  "Index, cap, bank and reason are recorded for each occurrence. [Pilot diagnostics](pilot-diagnostics.csv)"]
    else:
        lines += ["", "All 80 primary medians had 50 usable draws, all passing-only medians had a contributor, "
                  "and every first-stage fit had a statistical decision."]
    (args.out / "REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
