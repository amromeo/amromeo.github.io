# SPDX-License-Identifier: GPL-3.0-only
"""Rebuild the four study routes from each public fit and its 100 before-mask draw records."""

import csv
from hashlib import file_digest
import math
from pathlib import Path
import re
from statistics import median

BANKS = ("A", "B")
DESIGN_SHA = "5f4a78f982c3a882c8cadb742ada0f9b191e7a4ba658651fa5a69dddc48652c7"
GENERATOR_SHA = "b6b754278afb630ff707b154eff48e8517dd60b52f36a410db558b0401215964"
SHA = re.compile(r"^[0-9a-f]{64}$")
STATUSES = {"reported", "declined", "withheld"}
NAN = float("nan")


def ensure(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    with path.open("rb") as f:
        return file_digest(f, "sha256").hexdigest()


def read(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def number(v):
    return NAN if v in ("", "NA", None) else float(v)


def close(x, y):
    return x == y or (math.isnan(x) and math.isnan(y)) or math.isclose(x, y, rel_tol=2e-15, abs_tol=2e-15)


def get_ends(row, prefix=""):
    return number(row[prefix + "lower"]), number(row[prefix + "upper"])


def valid(limits):
    lower, upper = limits
    return math.isfinite(lower) and math.isfinite(upper) and lower < upper


def missing(limits):
    return all(math.isnan(v) for v in limits)


def load_design(design, hashes):
    ensure(sha(design) == DESIGN_SHA, "Stage 2 design digest differs")
    design_rows = [r for r in read(design) if int(r["N"]) == 50000]
    roster = {int(r["j"]): r for r in design_rows}
    reference = read(hashes)
    by_j = {int(r["j"]): r for r in reference}
    ensure(len(design_rows) == len(roster) == len(reference) == len(by_j) == 1440
           and set(roster) == set(by_j), "Original 1,440 design rows and hash IDs must match")
    for j, row in roster.items():
        for name in ("routine_hash", "anchor_hash"):
            ensure(SHA.fullmatch(by_j[j][name]), f"Malformed reference {name} for {j}")
            row[name] = by_j[j][name]
    return roster


def read_receipts(folder, j, design):
    rows = read(folder / "receipt.csv")
    ensure(len(rows) == 1, f"Case receipt count for {j}")
    receipt = rows[0]
    ensure(int(receipt["j"]) == j and int(receipt["routine_n"]) == 50000
           and int(receipt["anchor_n"]) == 40 and close(float(receipt["rounding_unit"]), float(design["h"]))
           and all(receipt[key] == design[key] for key in ("routine_hash", "anchor_hash")),
           f"Generation or archive identity mismatch on {j}")
    return receipt


def check_draws(j, bank, rows):
    ordered = sorted(rows, key=lambda r: int(r["draw"]))
    ensure([int(r["draw"]) for r in ordered] == list(range(1,51)), f"Incomplete draws at {j}/{bank}")
    base = 100_000_000 if bank == "A" else 200_000_000
    for i, row in enumerate(ordered, 1):
        ensure(int(row["sample_seed"]) == base + 100*j + i and SHA.fullmatch(row["positions_hash"]),
               f"Sampler identity at {j}/{bank}/{i}")
        ensure(row["first_status"] == row["status"] and row["status"] in STATUSES,
               f"Invalid subsample decision at {j}/{bank}/{i}")
        if row["status"] == "reported":
            ensure(valid(get_ends(row)) and get_ends(row) == get_ends(row,"raw_"),
                   f"Subsample public/raw limits differ at {j}/{bank}/{i}")
        if row["status"] == "declined":
            ensure(row["lrt_flag"].upper() == "TRUE" or row["lower_flag"].upper() == "TRUE",
                   f"Subsample decline without either flag at {j}/{bank}/{i}")
    pairs = [p for p in (get_ends(r,"raw_") for r in ordered) if valid(p)]
    median_limits = tuple(median(x) for x in zip(*pairs)) if pairs else (NAN,NAN)
    decline = sum(r["first_status"] == "declined" for r in ordered)
    passed = sum(r["first_status"] == "reported" for r in ordered)
    technical = len(ordered) - passed - decline
    majority = "declined" if decline >= 26 else ("reported" if decline + technical <= 25 else "withheld")
    return median_limits, dict(contributors=len(pairs), votes_pass=passed, votes_decline=decline,
                               votes_technical=technical, majority=majority)


def read_case(folder, design):
    j = int(design["j"])
    receipt = read_receipts(folder, j, design)
    fits = read(folder / "fits.csv")
    routes = read(folder / "routes.csv")
    ensure(len(fits) == 101 and len(routes) == 2 and all(int(r["j"]) == j for r in fits + routes)
           and {r["bank"] for r in routes} == set(BANKS), f"Records/route count differ for {j}")
    full, = [r for r in fits if r["fit_type"] == "full"]
    ensure(full["first_status"] in STATUSES and full["status"] in STATUSES, f"Unexpected public status for {j}")
    f_ends = get_ends(full)
    f_return = full["status"] == "reported"
    ensure(f_return or missing(f_ends), f"Package nonreturn retained limits for {j}")
    outcomes = []
    for bank in BANKS:
        rr, = [r for r in routes if r["bank"] == bank]
        draws = [r for r in fits if r["fit_type"] == "subsample" and r["bank"] == bank]
        m, votes = check_draws(j, bank, draws)
        ensure(all(int(rr[k]) == v for k,v in votes.items() if k != "majority")
               and rr["majority"] == votes["majority"], f"Vote discrepancy at {j}/{bank}")
        ensure(all(close(x,y) for x,y in zip(m,get_ends(rr,"median_"))), f"Median discrepancy for {j}/{bank}")
        if full["first_status"] == "reported":
            step, limits = "a", m
        elif f_return and valid(f_ends):
            step, limits = "b", f_ends
        elif votes["majority"] == "reported" and valid(m):
            step, limits = "c", m
        else:
            step, limits = "d", (NAN,NAN)
        ensure(rr["step"] == step and all(close(x,y) for x,y in zip(limits, get_ends(rr)))
               and (rr["status"] == "reported") == valid(limits),
               f"Policy differs from independently assembled route at {j}/{bank}")
        for key, source in [("first_status","first_status"), ("package_status","status"),
                            ("package_method","method"), ("package_reason","reason"),
                            ("first_reason","first_reason"), ("fallback_status","fallback_status"),
                            ("fallback_reason","fallback_reason")]:
            ensure(rr[key] == full[source], f"{key} differs for {j}/{bank}")
        ensure(all(close(x,y) for x,y in zip(f_ends,get_ends(rr,"package_"))), f"Package endpoints differ for {j}/{bank}")
        frac = float(design["fraction"])
        item = dict(j=j, bank=bank, cell_id=(j-1)//8+1, replicate=int(design["replicate"]), family=design["family"],
                    disease=design["disease"], fraction=frac, pathology="<=20%" if frac <= .2 else ">20%",
                    rounding_fraction=float(design["rounding_fraction"]), rounding_unit=float(design["h"]),
                    routine_hash=receipt["routine_hash"], anchor_hash=receipt["anchor_hash"],
                    first_status=full["first_status"], first_reason=full["first_reason"],
                    package_status=full["status"], package_reason=full["reason"], package_method=full["method"],
                    fallback_status=full["fallback_status"], fallback_reason=full["fallback_reason"],
                    first_lrt_flag=full["lrt_flag"], first_lower_flag=full["lower_flag"], first_p=full["p"],
                    policy_status=rr["status"], policy_reason=rr["reason"], step=step, **votes)
        for k in ("lambda","mu","sigma","truth_lower","truth_upper"):
            item[k] = float(design[k])
        for name, pair, reported in (("policy",limits,rr["status"]=="reported"),
                                     ("package",f_ends,f_return),
                                     ("through_b",limits if step in ("a","b") else (NAN,NAN),valid(limits) and step in ("a","b"))):
            item.update({name+"_lower":pair[0],name+"_upper":pair[1],name+"_returned":reported,
                         name+"_invalid":reported and not valid(pair)})
        outcomes.append(item)
    return outcomes


def collect(records, roster):
    receipts = list(records.rglob("receipt.csv"))
    cases = {int(p.parent.name): p.parent for p in receipts if p.parent.name.isdecimal()}
    ensure(len(cases) == len([p for p in receipts if p.parent.name.isdecimal()]) == len(roster)
           and cases.keys() == roster.keys(), "Unexpected, incomplete, or duplicated case directories")
    return [r for j in sorted(roster) for r in read_case(cases[j],roster[j])]
