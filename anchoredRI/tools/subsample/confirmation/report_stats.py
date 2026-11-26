# SPDX-License-Identifier: GPL-3.0-only
"""Paired metrics and whole-cell uncertainty for the Task H confirmation."""

from collections import Counter
import math
import numpy as np
from scipy.stats import beta
from report_records import ensure

N_BOOT = 2000
SEED = 37000000
NAN = float("nan")


def avg(a):
    return float(np.mean(a)) if len(a) else NAN


def middle(a):
    return percentile(a, .5) if len(a) else NAN


def percentile(a, q):
    """Linear percentile of defined values, including valid unbounded errors."""
    a = np.sort(np.asarray(a, dtype=float))
    a = a[~np.isnan(a)]
    if not len(a):
        return NAN
    h = (len(a)-1)*q
    left, right = a[int(math.floor(h))], a[int(math.ceil(h))]
    if left == right or h.is_integer():
        return float(left)
    if np.isneginf(left):
        return NAN if np.isposinf(right) else -math.inf
    if np.isposinf(right):
        return math.inf
    return float(left + (h % 1)*(right-left))


def bootstrap(cell, pairs, cand, ref, loss_diff):
    ids, inverse = np.unique(cell, return_inverse=True)
    k = len(ids)
    # One stream per reported population; both banks, both estimators and every
    # contrast within it use identical selections from the same ordered cells.
    samples = np.random.Generator(np.random.PCG64(SEED)).integers(k, size=(N_BOOT,k))
    def cell_sum(vals):
        return np.bincount(inverse, weights=vals.astype(float), minlength=k)[samples].sum(axis=1)
    denominator = cell_sum(pairs)
    assigned = cell_sum(np.ones_like(pairs))
    def common_mean(x):
        y = cell_sum(np.where(pairs & np.isfinite(x), x, 0.))
        inf_n = cell_sum(pairs & np.isposinf(x))
        with np.errstate(divide="ignore",invalid="ignore"):
            averages = y / denominator
        averages[inf_n > 0] = math.inf
        return averages
    with np.errstate(invalid="ignore",divide="ignore"):
        box = {"uncapped":common_mean(cand)-common_mean(ref),
               "capped":common_mean(np.minimum(cand,5))-common_mean(np.minimum(ref,5)),
               "all_assignment":cell_sum(loss_diff)/assigned}
    metrics = {}
    informative = len(np.unique(cell[pairs]))
    for key, draws in box.items():
        defined = ~np.isnan(draws)
        valid = int(defined.sum())
        enough = (k if key == "all_assignment" else informative) >= 2 and valid > 0
        metrics.update({key+"_bootstrap_valid":valid,key+"_bootstrap_undefined":N_BOOT-valid,
                        key+"_ci_lower":percentile(draws,.025) if enough else NAN,
                        key+"_ci_upper":percentile(draws,.975) if enough else NAN})
    return metrics


def groups(data):
    yield "all","all",data
    for axis in ("pathology","fraction","family","disease","rounding_fraction"):
        for value in sorted({r[axis] for r in data}):
            yield axis, str(value), [r for r in data if r[axis]==value]


def metric_row(rows, estimator, bank, axis, group):
    def ar(name, dtype=float):
        return np.asarray([r[name] for r in rows], dtype=dtype)
    p, c = ar("package_returned",bool), ar(estimator+"_returned",bool)
    pairs = p & c
    pd, cd = ar("package_d"), ar(estimator+"_d")
    ensure(not np.any(np.isnan(pd[p])) and not np.any(np.isnan(cd[c])), "Missing score on returned interval")
    ref, cand = pd[pairs], cd[pairs]
    finite = np.isfinite(ref) & np.isfinite(cand)
    with np.errstate(invalid="ignore"):
        d = cand-ref
        difference = avg(cand)-avg(ref)
    distinct = ~np.isnan(d)
    p_loss=np.where(p,np.minimum(pd,5),5)
    c_loss=np.where(c,np.minimum(cd,5),5)
    p_ext=ar("package_extreme",bool); c_ext=ar(estimator+"_extreme",bool)
    row = dict(bank=bank, estimator=estimator, axis=axis, group=group, assigned=len(rows),
               cells=len(set(ar("cell_id",int))), package_returns=int(p.sum()), policy_returns=int(c.sum()),
               lost_package_returns=int((p & ~c).sum()), added_returns=int((~p & c).sum()),
               package_invalid_returns=int(ar("package_invalid",bool).sum()), policy_invalid_returns=int(ar(estimator+"_invalid",bool).sum()),
               package_extremes=int(p_ext.sum()), policy_extremes=int(c_ext.sum()),
               package_mean_d_released=avg(pd[p]), policy_mean_d_released=avg(cd[c]), common_returns=int(pairs.sum()),
               package_common_mean_d=avg(ref), policy_common_mean_d=avg(cand), common_uncapped_difference=difference,
               package_common_capped=avg(np.minimum(ref,5)),policy_common_capped=avg(np.minimum(cand,5)),
               common_capped_difference=avg(np.minimum(cand,5)-np.minimum(ref,5)), common_finite_pairs=int(finite.sum()),
               common_finite_pair_uncapped_difference=avg(d[finite]), common_package_infinite_d=int(np.isinf(ref).sum()),
               common_policy_infinite_d=int(np.isinf(cand).sum()), common_both_infinite_d=int((~np.isfinite(cand)&~np.isfinite(ref)).sum()),
               median_paired_change=middle(d[distinct]),median_absolute_paired_change=middle(np.abs(d[distinct])),
               paired_change_defined=int(distinct.sum()), wins=int((cand<ref).sum()),ties=int((cand==ref).sum()),
               losses=int((cand>ref).sum()), policy_all_assignment=avg(c_loss),package_all_assignment=avg(p_loss),
               all_assignment_difference=avg(c_loss-p_loss), informative_pair_cells=len(set(ar("cell_id",int)[pairs])))
    for side in ("lower","upper"):
        row["package_extreme_"+side]=int(np.sum(p & (np.abs(ar("package_z_"+side))>5)))
        row["policy_extreme_"+side]=int(np.sum(c & (np.abs(ar(estimator+"_z_"+side))>5)))
    row.update(bootstrap(ar("cell_id",int),pairs,cd,pd,c_loss-p_loss))
    return row


def rescue_metrics(rows):
    eligible = [r for r in rows if not r["package_returned"]]
    rescues = [r for r in eligible if r["step"]=="c"]
    n=len(rescues)
    extreme=sum(r["policy_extreme"] for r in rescues)
    errors=[r["policy_d"] for r in rescues]
    result=dict(package_nonreturns=len(eligible),rescues=n,recovered_fraction=n/len(eligible) if eligible else NAN,
                resolved_majorities=sum(r["majority"]!="withheld" for r in eligible),
                unresolved_majorities=sum(r["majority"]=="withheld" for r in eligible),
                package_statistical_nonreturns=sum(r["package_status"]=="declined" for r in eligible),
                package_technical_nonreturns=sum(r["package_status"]=="withheld" for r in eligible),
                rescue_mean_d=avg(errors),rescue_median_d=middle(errors),rescue_max_d=max(errors,default=NAN),
                rescue_endpoint_extremes=extreme,rescue_extreme_rate=extreme/n if n else NAN,
                exact95_low=(0. if extreme==0 else float(beta.ppf(.025,extreme,n-extreme+1))) if n else NAN,
                exact95_high=(1. if extreme==n else float(beta.ppf(.975,extreme+1,n-extreme))) if n else NAN)
    for key in ("votes_pass","votes_decline","votes_technical","contributors"):
        a=[r[key] for r in rescues]
        result.update({key+"_min":min(a,default=NAN),key+"_median":middle(a),key+"_max":max(a,default=NAN)})
    for side in ("lower","upper"):
        a=[r["policy_z_"+side] for r in rescues]
        result.update({side+"_signed_median":middle(a),side+"_signed_min":min(a,default=NAN),side+"_signed_max":max(a,default=NAN)})
    return result


def gate_rows(metrics,rescues,bank):
    use=[r for r in metrics if r["bank"]==bank and r["estimator"]=="policy"]
    overall,=[r for r in use if r["axis"]=="all"]
    strat=[r for r in use if r["axis"]=="pathology"]
    ensure(len(strat)==2 and {r["group"] for r in strat}=={"<=20%",">20%"}, "Missing pathology strata")
    cp=overall["common_uncapped_difference"]
    low_obs=cp<0
    improved=not math.isnan(cp) and overall["uncapped_ci_upper"]<0 and not math.isnan(overall["uncapped_ci_lower"])
    no_worse=not math.isnan(overall["common_capped_difference"]) and overall["capped_ci_upper"]<.02 and not math.isnan(overall["capped_ci_lower"])
    accuracy=low_obs and all(r["common_returns"]>0 and r["common_capped_difference"]<=.02 for r in strat)
    coverage=overall["lost_package_returns"]==0 and overall["policy_invalid_returns"]==0 and all(
        r["policy_extremes"]<=r["package_extremes"] for r in [overall,*strat])
    all_assignment=overall["all_assignment_difference"]<0
    rescue=rescues
    step_c=("insufficient evidence" if rescue["rescues"]<10 else
            ("pass" if rescue["rescue_mean_d"]<.5 and rescue["rescue_endpoint_extremes"]==0 else "fail"))
    # The accuracy-row point-estimate gate and the three inferential/descriptive
    # labels have distinct estimands. See the approved protocol.
    complete=("fail" if not (accuracy and coverage and all_assignment) or step_c=="fail"
              else ("insufficient evidence" if step_c=="insufficient evidence" else "pass"))
    yn=lambda v:"pass" if v else "fail"
    return dict(bank=bank,primary=bank=="A",lower_observed_mean_error=low_obs,
                improved_accuracy=improved,no_worse=no_worse,accuracy_gate=yn(accuracy),
                coverage_and_extremes_gate=yn(coverage),all_assignment_gate=yn(all_assignment),
                step_c_evidence=step_c,complete_four_step_gate=complete)


def tables(rows):
    metrics, rescues, reasons, routes, gates = [],[],[],[],[]
    for bank in ("A","B"):
        bank_rows=[r for r in rows if r["bank"]==bank]
        for axis, group, subset in groups(bank_rows):
            for estimator in ("policy","through_b"):
                metrics.append(metric_row(subset,estimator,bank,axis,group))
            for step in "abcd":
                cases = [r for r in subset if r["step"]==step]
                released=[r for r in cases if r["policy_returned"]]
                routes.append(dict(bank=bank,axis=axis,group=group,step=step,n=len(cases),released=len(released),
                                   technical_nonreturns=sum(r["policy_status"]=="withheld" for r in cases),
                                   statistical_nonreturns=sum(r["policy_status"]=="declined" for r in cases),
                                   mean_returned_d=avg([r["policy_d"] for r in released]),
                                   any_extreme=sum(r["policy_extreme"] for r in released)))
            rescue=dict(bank=bank,axis=axis,group=group,**rescue_metrics(subset))
            rescues.append(rescue)
            if axis=="all":
                for key,count in Counter((r["first_reason"],r["package_status"],r["package_reason"],
                                          r["fallback_status"],r["fallback_reason"],r["step"],r["majority"])
                                         for r in subset if not r["package_returned"]).items():
                    reasons.append(dict(bank=bank,first_reason=key[0],package_status=key[1],package_reason=key[2],
                                        fallback_status=key[3],fallback_reason=key[4],step=key[5],majority=key[6],n=count))
        total,=[r for r in rescues if r["bank"]==bank and r["axis"]=="all"]
        gates.append(gate_rows(metrics,total,bank))
    return metrics,rescues,reasons,routes,gates
