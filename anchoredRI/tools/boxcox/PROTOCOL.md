# Task E — Box-Cox shape and mixture capacity

## Status and frozen baseline

Approved for implementation by the user on 22 November 2026. The following protocol was committed and reviewed before coding or fitting. Create `study/boxcox` from main commit `55eee86b3f4111197b0043cc062ef0b781dc2057`. This commit defines current anchoredRI for this study. All new code, manifests, checks and reports belong under `anchoredRI/tools/boxcox/`. Record source hashes and package versions. Keep `R/` and the simulation design at the pinned commit.

Current main requires two routine settings. RIbench supplies one vector. For RIbench, use a study-local one-group adaptation matching the already executed C2 adapter at `f72ea14fc7978bf2478f54c87548cb6c127fa766`. Verify that this adaptation reproduces the archived baseline fits before comparing arms. For the simulation grid, use the pinned two-setting estimator directly. Describe this distinction in every report. A mismatch stops the affected comparison for investigation.

Reuse the independent-anchor small-N outputs from run 45096614442 and N = 50,000 outputs from run 45101564181, including the existing anchors-only and refineR results. Validate dataset keys, data hashes, anchor seeds and source provenance before joining. No shared-anchor results from the original RIbench run enter the comparison.

## Models

Arm A has one Box-Cox transform shared by all mixture components. With offset d, transform x as g(x) = ((x-d)^lambda - 1)/lambda, with the continuous log(x-d) limit at lambda = 0. Component 1 receives the screened anchors. Each routine setting has its own mixture weights. RIbench uses one setting and the simulation grid uses two.

Estimate lambda jointly with component means, SDs and weights. Compare K = 1, 2 and 3 and two offset specifications by BIC. The unshifted specification fixes d = 0 and is eligible when all routine and anchor values are positive. The shifted specification estimates d below every observed value. Follow the existing shifted fitter's pooled-data spread calculation, including its fallback, with bounds min(x,anchor) + spread * c(-5,-0.001). Count one extra parameter for a fitted offset. These are additive-offset Box-Cox models with the convention x-d.

Proposed numerical settings fixed by this protocol are lambda in [-2,2], starts at -1, -0.5, 0, 0.5, 1 and 2, and the existing shifted fitter's offset starts. Use the frozen fitters' component initializations, numerical scale safeguards and optimizer tolerances where applicable. Include frozen identity/log/shifted-log solutions as valid starts when they lie in the candidate domain. Record all starts, convergence codes, attained likelihoods and bound hits. Unit changes must be accounted for exactly in transformations and original-scale likelihoods. Test the transformation and objective gradients numerically before benchmark runs. Any numerical revision after results are inspected must be versioned and identified as exploratory.

The log likelihood sums the routine mixture log density, the anchors' component-1 log density, and (lambda-1)*log(x-d) for every routine value and every anchor. Retain multiplicities when compressing repeated values. BIC is -2*logLik + p*log(n_routine+n_anchor), with p = 2*K + G*(K-1) + 1 for an estimated lambda, plus 1 when d is estimated. Fixed lambda or offset parameters do not count as estimated parameters.

Use the same untruncated normal-on-transformed-scale likelihood convention as the current estimator. For nonzero lambda, the inverse Box-Cox has a restricted domain. Record this limitation and treat a requested endpoint outside that domain as an unavailable interval. Do not clip or substitute endpoints. Reject a candidate whose likelihood or fitted parameters are invalid. Select the lowest-BIC valid candidate, then apply the safeguard to that selected fit without selecting a replacement model.

Arm B keeps identity, log and shifted-log transformations and their original-scale BIC calculation, while extending the component search to K = 1 through 5. Retain the frozen starts and settings for K <= 3 and use the same initialization construction for K = 4 and 5. The parameter count is 2*K + G*(K-1), plus 1 for shifted-log. Keep component 1 anchored. Report component counts and cap hits for both arms and the baseline.

## Equivalence checks before outcome analysis

For Arm A with lambda fixed at 1 and d = 0, map the additive constant in g(x) = x-1 into all component means, initializations and bounds. This special equivalence test accepts the identity fitter's full real-valued input domain. With lambda fixed at 0 and d = 0, reproduce the current log fitter on positive inputs. Match component counts, fitted limits and original-scale log likelihoods using the same starts and optimizer settings. Fixed-lambda BIC excludes a fitted shape parameter.

Check each component-count candidate as well as the selected result. Limits must agree within 1e-8 absolute plus 1e-10 relative tolerance, and log likelihoods within 1e-8 per observation. Record maximum discrepancies. Stop on a mismatch before interpreting arm performance. Apply these checks to every study dataset where the corresponding frozen transform is eligible. Verify Arm B restricted to K <= 3 against the frozen auto estimator with the same tolerances.

## Datasets and paired inputs

Run both arms on all 2,880 RIbench N <= 5,000 datasets and all 1,440 N = 50,000 datasets. Draw 40 anchors with seed 17000000 + Index using the original C2 inverse Box-Cox construction, shift, removal of negative/nonfinite draws, and rounding. Do not replenish removed draws. Check against archived anchor hashes and pass identical vectors to all anchor-based methods. Use the archived routine-data hashes and RIbench 1.0.2 scoring settings. CRP uses the one-sided 95th percentile. Other analytes use the central 95% interval. Report lambda = 0.432 separately.

Also run all 45 simulation-grid cases on replicates 1–20 using the frozen stress_data() generator and its original anchor seeds. Preserve its intentional anchor reuse across scenarios within a family/replicate. Use paired comparisons and Monte Carlo SE based on the 20 replicate differences within each cell. These are previously explored datasets, so label this as development and regression checking. This protocol does not convert the previously examined datasets into an independent confirmation sample.

## Study-local safeguard

Implement the two checks under tools/boxcox/, independently of Task F. On the selected transform and offset, compare the screened anchors under component 1's fitted mean and SD with their own normal MLE mean and SD. Compute LR = max(0, 2*(logLik_own-logLik_component1)) and the upper-tail chi-squared p-value with 2 df. Hold lambda and offset fixed in both terms. For Arm A this calculation uses the fitted Box-Cox scale. Transformation Jacobians cancel. The p-value is a nominal diagnostic after model fitting and selection.

Decline if p < 0.01 OR the estimated lower limit is strictly below the minimum routine value. Skip the lower-limit check for a one-sided interval. An unavailable LRT or interval is a technical withholding, recorded separately. When declined, reported limits are missing, fitted limits remain diagnostic, and both flags and the reason are retained. Keep the anchored estimate and anchors-only estimate as diagnostics. Apply the same study-local checks to the current anchoredRI baseline. Anchors-only and refineR remain comparators with their original reporting behavior.

Report performance before the safeguard over all attempts and after it among reported intervals. Show the reporting denominator, decline rate, technical failures, and error-group flag counts. Retain rejected estimates for the unfiltered analysis. Do not replace declined fits with anchors-only intervals.

Document that heavily mixed populations containing many mildly abnormal patients, including ICU-only populations, can produce inaccurate intervals that pass both checks.

## Metrics and outputs

Preserve RIbench's native standardized endpoint deviations and capped mean absolute z deviation. CRP contributes its upper endpoint only. Report uncapped and capped results, failures and the number scored. The native implausible-result criterion used for the success threshold is per-dataset mean absolute z deviation > 5. Also export every dataset with either endpoint |z| > 5, making the distinction explicit.

For each arm versus current anchoredRI, anchors-only and refineR, report overall and distribution type by pathological fraction <=20% versus >20% results. Include paired mean differences, wins, ties and losses, numbers of complete pairs, and decline rates. Give before-safeguard comparisons on the common dataset set and after-safeguard comparisons on the same retained pairs, with denominators. Report each arm separately.

For the simulation grid, retain stress_metrics() and stress_row() metrics, including error as a percentage of the true interval width. For every family by overlap by affected-fraction cell, report paired changes from the frozen baseline, wins and MCSE = SD(differences)/sqrt(20). Flag cells with a paired mean error increase greater than 2 MCSE as evidence of regression, and show all cells without treating unflagged changes as proof of equivalence.

Save per-dataset selected transform, lambda, offset, component count, component parameters and weights, BIC, convergence and bound hits, anchor seed/hash, limits, routine minimum, LRT statistic/p, lower-limit flag, decline status/reason, z deviations and runtime. Save per-candidate fits needed to audit selection. Separate generation, fitting and total elapsed runtime. Preserve source IDs for every reused comparator result.

## Fixed success criteria

Evaluate each arm separately. Primary success uses the unfiltered, before-safeguard native capped metric over the specified dataset set, so selective declines cannot create a primary success.

1. At N = 50,000, capped mean absolute z deviation must be no greater than existing refineR's value on the same 1,440 datasets, and at most 1% of attempts may be implausible (at most 14 of 1,440).
2. At N <= 5,000, capped mean absolute z deviation must be no greater than current anchoredRI's value on the same 2,880 datasets.

Use full-precision archived comparators, not rounded report values. Report any missing/failed fits explicitly. If an arm has missing fits, give the native score and an additional all-attempt score assigning missing fits the cutoff of 5; require the all-attempt score to meet the threshold and count technical fit failures with implausible results for the 1% check. Report safeguard-retained performance and coverage separately from the primary verdict. Show simulation regressions alongside the verdict.

The fixed thresholds, numerical choices and interpretation of the primary criterion were reviewed before implementation. The user authorized proceeding on 22 November 2026. They remain unchanged.
