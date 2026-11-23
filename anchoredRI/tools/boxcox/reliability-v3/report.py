"""Compare the completed round-3 checkpoints. This script does not fit or score an interval."""
import argparse
import csv
import math
import shutil
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
IDS = set(range(1, 157))
HOSTS = ("local", "Actions")
TOL = 0.005


def read(path):
    with path.open(newline="") as source:
        return list(csv.DictReader(source))


def write(path, rows):
    if not rows:
        raise ValueError(f"No rows for {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as destination:
        out = csv.DictWriter(destination, fieldnames=list(rows[0]), lineterminator="\n")
        out.writeheader()
        out.writerows(rows)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def number(x):
    return math.nan if x in ("NA", "", None) else float(x)


def integer(row):
    return int(row["Index"])


def yes(value):
    return value in (True, 1, "1", "TRUE")


def as_key(row):
    return row["transform"], int(row["k"])


def one_per_index(rows):
    result = {integer(row): row for row in rows}
    require(len(result) == len(rows) and result.keys() == IDS, "Expected exactly Hb 1 through 156.")
    return result


def selected(rows):
    return one_per_index([r for r in rows if yes(r["selected"])])


def compare(left, right):
    lower_l, upper_l = number(left["lower"]), number(left["upper"])
    lower_a, upper_a = number(right["lower"]), number(right["upper"])
    width = upper_a - lower_a
    finite = all(math.isfinite(x) for x in (lower_l, upper_l, lower_a, upper_a)) and width > 0
    delta = max(abs(upper_l - upper_a), abs(lower_l - lower_a)) if finite else math.nan
    return dict(local_transform=left["transform"], local_k=left["k"], actions_transform=right["transform"],
                actions_k=right["k"], model_match=as_key(left) == as_key(right),
                local_lower=lower_l, local_upper=upper_l, actions_lower=lower_a, actions_upper=upper_a,
                actions_width=width, finite_pair=finite, max_endpoint_difference=delta,
                fraction_of_width=delta / width if finite else math.nan, limit_match=finite and delta <= TOL*width)


def best(rows):
    # Same family and component order as the fitted model library.
    return min(rows, key=lambda r: (number(r["bic"]), r["transform"] == "boxcox_shifted", int(r["k"])))


def matched_manifest(local, actions):
    a = {r["file"]: r["sha256"] for r in read(local)}
    b = {r["file"]: r["sha256"] for r in read(actions)}
    require(a == b and len(a) > 0, f"Source disagreement in {local}")
    return a


def by_fit(rows):
    return {(integer(r), *as_key(r)): r for r in rows}


def baseline_checks(local, actions, out):
    summary, all_audits = [], []
    for host, path in zip(HOSTS, (local, actions)):
        candidates = read(path/"base_tab.csv")
        zero = read(path/"zero.csv")
        audit = read(path/"nesting.csv")
        require(len(by_fit(candidates)) == len(candidates) == len(IDS)*9, "Incomplete baseline candidate table.")
        require(len(zero) == 468 and len(audit) == 1404, "Incomplete baseline audit.")
        counts = Counter(r["status"] for r in audit)
        grouped = defaultdict(list)
        for row in audit:
            grouped[integer(row)].append(row)
        for i, rows in sorted(grouped.items()):
            checks = [r for r in rows if r["status"] == "checked"]
            all_audits.append(dict(Index=i, runtime=host, checked=len(checks),
                             excluded=sum(r["status"] == "structurally excluded" for r in rows),
                             unavailable=sum(r["status"] == "unavailable" for r in rows),
                             violations=sum(yes(r["violation"]) for r in rows),
                             largest_checked_loss=max(number(r["loss"]) for r in checks)))
        summary.append(dict(runtime=host, datasets=len(IDS), candidates=len(candidates),
                            checked=counts["checked"], excluded=counts["structurally excluded"],
                            unavailable=counts["unavailable"], violations=sum(yes(r["violation"]) for r in audit),
                            zero_vectors=len(zero), zero_vectors_in_fitted_range=sum(yes(r["zero_eligible"]) for r in zero),
                            max_zero_loglik_difference=max(number(r["absolute_difference"]) for r in zero),
                            max_zero_floor_difference=max(number(r["physical_floor_difference"]) for r in zero),
                            min_sd_to_floor=min(number(r["minimum_sd_to_floor"]) for r in candidates),
                            min_sd_in_measurement_units=min(number(r["minimum_original_sd"])/number(r["rounding_unit"])
                                                           for r in candidates)))
        for name in ("zero", "base_tab", "nesting", "arm_tab", "continuation-blocks"):
            shutil.copyfile(path/(name+".csv"), out/(host+"-"+name+".csv"))
    write(out/"baseline-audit-by-dataset.csv", all_audits)
    write(out/"baseline-audit-summary.csv", summary)
    return summary


def armb_checks(armb, out):
    datasets = read(armb/"datasets.csv")
    require(len({integer(d) for d in datasets}) == len(datasets) == 540, "Incomplete Arm B dataset set.")
    stages = [("N <= 5,000", [r for r in datasets if r["stage"] == "small"]),
              ("N = 50,000", [r for r in datasets if r["stage"] == "50000"]), ("All", datasets)]
    candidates = read(armb/"candidates.csv")
    audit = read(armb/"nesting.csv")
    zero = read(armb/"zero.csv")
    details = []
    for label, rows in stages:
        ids = {integer(r) for r in rows}
        ze = [r for r in zero if integer(r) in ids]
        cand = [r for r in candidates if integer(r) in ids]
        au = [r for r in audit if integer(r) in ids]
        checked = sum(r["status"] == "checked" for r in au)
        violations = sum(yes(r["violation"]) for r in au)
        require(checked == sum(int(r["checked"]) for r in rows)
                and violations == sum(int(r["violations"]) for r in rows), "Arm B counts do not reconcile.")
        details.append(dict(group=label, datasets=len(rows), candidates=len(cand),
                            old_flagged=sum(int(r["old_violations"]) > 0 for r in rows),
                            new_flagged=sum(int(r["violations"]) > 0 for r in rows), comparisons=checked,
                            violations=violations, unavailable=sum(int(r["unavailable"]) for r in rows),
                            excluded=sum(int(r["excluded"]) for r in rows), changed=sum(yes(r["changed"]) for r in rows),
                            zero_vectors=len(ze), max_zero_loglik_difference=max(number(r["absolute_difference"]) for r in ze),
                            max_zero_floor_difference=max(number(r["physical_floor_difference"]) for r in ze),
                            min_sd_to_floor=min(number(r["min_sd_to_floor"]) for r in rows),
                            min_sd_in_measurement_units=min(number(r["minimum_original_sd"])/number(r["rounding_unit"])
                                                           for r in cand)))
    require([r["datasets"] for r in details] == [396, 144, 540], "Unexpected Arm B study groups.")
    write(out/"arm-b-summary.csv", details)
    write(out/"arm-b-datasets.csv", datasets)
    write(out/"arm-b-zero.csv", zero)
    for label, rows in [("candidates", candidates), ("nesting", audit)]:
        for k in range(1, 6):
            chunk = [r for r in rows if int(r["k"]) == k]
            write(out/"arm-b"/f"{label}-K{k}.csv", chunk)
    return details


def compare_round(work, out):
    phase = {"local": work/"local", "Actions": work/"actions"/"phase1"}
    cross = {"local": work/"cross-local", "Actions": work/"actions"/"cross"}
    for p in [*phase.values(), *cross.values(), work/"arm-b"]:
        require((p/"DONE.txt").is_file(), f"Missing completion marker at {p}")
    require("C2 image match TRUE" in (phase["Actions"]/"runtime.txt").read_text(), "Actions runtime unverified.")
    manifests = {}
    for a, b, filename in [(phase["local"], phase["Actions"], "sources.csv"),
                           (phase["local"], phase["Actions"], "native-sources.csv"),
                           (cross["local"], cross["Actions"], "sources.csv")]:
        for key, digest in matched_manifest(a/filename, b/filename).items():
            require(manifests.get(key, digest) == digest, f"Inconsistent source hash {key}")
            manifests[key] = digest
    write(out/"source-hashes.csv", [dict(file=k, sha256=v) for k,v in sorted(manifests.items())])
    for h in HOSTS:
        lines = (phase[h]/"runtime.txt").read_text().splitlines()
        (out/(h+"-runtime.txt")).write_text("\n".join(line.rstrip() for line in lines)+"\n")
    bsum = baseline_checks(phase["local"], phase["Actions"], out)
    arm_b = armb_checks(work/"arm-b", out)
    baseline, heads, initial, scores = {}, {}, {}, {}
    for h,p in phase.items():
        baseline[h] = selected(read(p/"base_tab.csv"))
        all_exchange = [r for f in sorted((p/"exchange").glob("*.csv")) for r in read(f)]
        require(len(by_fit(all_exchange)) == len(all_exchange) == 936, "Incomplete exchange.")
        heads[h] = selected(all_exchange)
        initial[h] = one_per_index(read(p/"initial.csv"))
        scores[h] = by_fit(read(p/"arm_tab.csv"))
        require(all(int(r["code"]) == 0 for r in heads[h].values()), "Unfinished Arm A winner.")
    old = one_per_index([r for r in read(HERE.parent/"reliability-v2/results/paired-replay.csv") if r["method"]=="armA"])
    baseline_pair, cont, old42 = [], [], set()
    for i in sorted(IDS):
        l, a = heads["local"][i], heads["Actions"][i]
        for z in ("routine_hash","anchor_hash"):
            require(l[z] == a[z] == old[i][z], f"{z} differs for Hb {i}")
        hash_pair = dict(routine_hash=l["routine_hash"], anchor_hash=l["anchor_hash"])
        baseline_pair.append(dict(Index=i, **compare(baseline["local"][i],baseline["Actions"][i]), **hash_pair))
        initial_pair = compare(initial["local"][i], initial["Actions"][i])
        o = old[i]
        for col0, col1 in [("model_match","model_match"),("limit_match","within_0_5pct")]:
            require(initial_pair[col0] == yes(o[col1]), f"Original comparison disagrees for Hb {i}.")
        orig = yes(o["model_match"]) and not yes(o["within_0_5pct"])
        if orig:
            old42.add(i)
        cont.append(dict(Index=i, original42=orig, old_model_match=yes(o["model_match"]),
                         old_limit_match=yes(o["within_0_5pct"]), old_fraction_of_width=o["fraction_of_width"],
                         old_code_local=initial["local"][i]["code"], old_code_actions=initial["Actions"][i]["code"],
                         final_code_local=l["code"], final_code_actions=a["code"],
                         local_loglik=l["loglik"], actions_loglik=a["loglik"], **compare(l,a), **hash_pair))
    require(len(old42)==42, "Frozen set of 42 could not be reproduced.")
    write(out/"baseline-paired.csv", baseline_pair)
    write(out/"arm-a-continued-pairs.csv", cont)
    expected = {int(r["Index"]) for r in cont if not (r["model_match"] and r["limit_match"])}
    outcomes, traces = [], []
    for h,p in cross.items():
        membership = one_per_index(read(p/"membership.csv"))
        require(expected=={i for i,r in membership.items() if yes(r["cross_tested"])}, "Cross-start membership disagrees.")
        if expected:
            rs=read(p/"outcomes.csv")
            require(all(r["receiver"]==h and int(r["code"])==0 for r in rs), "Reciprocal fit unresolved.")
            outcomes.extend(rs)
            traces.extend(read(p/"blocks.csv"))
    if outcomes:
        write(out/"arm-a-cross-solutions.csv", outcomes)
        write(out/"arm-a-cross-blocks.csv", traces)
    crossed=defaultdict(list)
    for r in outcomes: crossed[integer(r)].append(r)
    require(set(crossed) == expected, "Missing cross-start data.")
    classifications, retained = [], []
    for i in sorted(IDS):
        l, a = heads["local"][i], heads["Actions"][i]
        if i not in expected:
            retained.append(dict(Index=i, cross_tested=False, **compare(l,a)))
            continue
        rows=crossed[i]
        specs={as_key(l),as_key(a)}
        require({(r["receiver"],r["donor"],as_key(r)) for r in rows}
                == {(h,d,spec) for h in HOSTS for d in HOSTS for spec in specs}
                and len(rows)==len(specs)*4, f"Missing starts at Hb {i}.")
        chosen=as_key(best(rows))
        winners={(r["receiver"],r["donor"]):r for r in rows if as_key(r)==chosen}
        reciprocal_l=winners["local","Actions"]
        reciprocal_a=winners["Actions","local"]
        c=compare(reciprocal_l,reciprocal_a)
        require(c["finite_pair"], "Missing cross-start limits.")
        dll=abs(number(reciprocal_l["ll"])-number(reciprocal_a["ll"]))
        cls="missed better fit" if dll>1 else ("weak identification" if not c["limit_match"] else "concordant")
        gains={}
        for h,tag in (("local","local"),("Actions","actions")):
            other=next(x for x in HOSTS if x!=h)
            own_stored=number(scores[h][(i,*chosen)]["ll"])
            gains[tag+"_foreign_minus_pre_exchange"]=number(winners[h,other]["ll"])-own_stored
            gains[tag+"_foreign_minus_own_restart"]=number(winners[h,other]["ll"])-number(winners[h,h]["ll"])
        classifications.append(dict(Index=i, original42=i in old42, competing_models=len(specs), common_transform=chosen[0],
                                    common_k=chosen[1], reciprocal_local_ll=reciprocal_l["ll"],
                                    reciprocal_actions_ll=reciprocal_a["ll"], abs_loglik_difference=dll,
                                    classification=cls, **c, **gains))
        best_by_host={}
        for h in HOSTS:
            # Retain the best result for every family/K, including saved candidates that were never nominated.
            pool=[r for (j,*_),r in scores[h].items() if j==i]+[r for r in rows if r["receiver"]==h]
            best_by_host[h]=best(pool)
        retained.append(dict(Index=i,cross_tested=True,**compare(best_by_host["local"],best_by_host["Actions"])))
    if classifications:
        write(out/"arm-a-cross-decisions.csv", classifications)
    write(out/"arm-a-after-retaining-best.csv", retained)
    return bsum,arm_b,baseline_pair,cont,classifications,retained,outcomes,manifests


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work",type=Path,help="Directory containing local, actions, cross-local and arm-b")
    parser.add_argument("output",type=Path)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    data=compare_round(args.work,args.output)
    from report_text import render
    render(args.work,args.output,*data)


if __name__=="__main__":
    main()
