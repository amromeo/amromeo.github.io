# Task B2 protocol for the routine-data separation check

**Status.** Awaiting review before implementation or study runs. The development choice is open.

Can the fitted mixture identify an inaccurate reference interval before it is reported? The primary result comes from fresh replicates 21–40. The rule succeeds if it declines at least two-thirds of estimates with error over 10% and no more than 5% of estimates with error under 5%.

## Datasets and frozen sources

| Phase | Replicates per case | Datasets | Use |
| --- | --- | ---: | --- |
| Development | 1–20 | 900 | Select one of 14 prespecified candidates |
| Confirmation | 21–40 | 900 | Evaluate the committed choice once |

Each phase covers the same 45 cases. The design crosses three healthy distributions (normal, lognormal and shifted lognormal), three degrees of overlap (separated, moderate and heavy), and five affected fractions (10%, 30%, 50%, 70% and 90%). Each dataset has 2,500 routine results per service and 40 independently screened healthy patients. The two service-specific affected fractions are five percentage points below and above the overall fraction. Within a family and replicate, the design reuses the same 40 anchors across cases. Each new replicate draws fresh anchors and routine patients.

Use frozen anchoredRI 0.2.0 and `tools/stress-design.R` from commit `f463a2cb84439ef9a9bb71e44f76661ceb0df00f`. Use Step B's `lrt()` in `tools/disagreement/run.R` and its `ao_normal()` helper in `tools/anchors-only/estimators.R` from commit `55eee86b3f4111197b0043cc062ef0b781dc2057`. All new study code belongs in `tools/disagreement-v2/`. Use `R/` and `tools/stress-design.R` as read-only inputs.

## Fit and equivalence check

Generate each dataset with `stress_data(case, replicate)` at the design defaults. For development, match case, replicate and seed to the existing 900 `anchored40` records in `assets/comparison-study/summary/comparison-results.zip`. Fit each dataset twice using the same routine values, service labels and 40 anchors.

1. Run frozen `fit_anchored_ri(..., transform = "auto")` with default shift bounds.
2. Run the same selection using study-local copies of `.fit_anchored_fixed` and `.fit_anchored_shifted`. The copies also return every selected component's mean and SD on its fitted internal scale and its weights in the two services. Keep the frozen fitting and BIC selection code, including the component numbering. Component 1 is the healthy component that receives the anchors.
3. Compare the two fits on **every dataset before scoring it**. The selected transform and component count must match exactly. Each of the two limits and the selected offset on the original measurement scale must be finite and agree within an absolute `1e-10`, using unrounded values. If a fit cannot be verified or any comparison fails, stop and identify the case, replicate, field and values from both fits.

## Primary separation statistic

Let $\mu_j$ and $\sigma_j$ be the fitted mean and SD of component $j$ on the selected model's internal scale. The fitted healthy interval on that scale has endpoints $L=\mu_1+\Phi^{-1}(0.025)\sigma_1$ and $U=\mu_1+\Phi^{-1}(0.975)\sigma_1$. For each disease component $j=2,\ldots,K$, compute its probability of falling inside that interval and its pooled routine weight.

$$
p_j=\Phi\!\left(\frac{U-\mu_j}{\sigma_j}\right)-\Phi\!\left(\frac{L-\mu_j}{\sigma_j}\right),\qquad
w_j=\frac{n_a w_{ja}+n_b w_{jb}}{n_a+n_b}.
$$

Here $n_a$ and $n_b$ count routine patients only. All components share the same internal scale, so these calculations require no back-transformation. The primary statistic is the fitted fraction of the disease mixture inside the healthy interval.

$$
O=\frac{\sum_{j=2}^{K} w_j p_j}{\sum_{j=2}^{K} w_j}.
$$

Compare fitted $O$ with `stress_design$disease_inside_true_ri`. The design values are approximately 5% for separated, 49% for moderate and 82% for heavy overlap. When BIC selects one component, record $O$ as undefined, keep the dataset in the analysis and use the LRT alone under either candidate rule.

## Secondary check and candidate rules

Use Step B's likelihood-ratio test (LRT) unchanged. On the selected transform and offset, it compares the normal log likelihood of the 40 anchors under their own maximum-likelihood mean and SD with their log likelihood under the parameters implied by the frozen healthy limits. Step B takes twice the log-likelihood difference, bounds it below at zero and uses the upper tail of a chi-squared distribution with 2 degrees of freedom. A finite $p<0.01$ raises the LRT flag. The bootstrap CI check is dropped.

| Candidate | Two or three components | One component |
| --- | --- | --- |
| a | Decline if $O>t$ | LRT alone |
| b | Decline if $O>t$ or LRT $p<0.01$ | LRT alone |

Evaluate each rule at $t$ equal to 0.50, 0.55, 0.60, 0.65, 0.70, 0.75 and 0.80. Use strict inequalities on unrounded values. As in Step B, an unavailable LRT withholds an interval when the decision requires the LRT. Record it as a technical withholding, separate from a positive LRT. For rule a with a defined $O$, the LRT is diagnostic only. The standalone comparator always uses the Step B LRT decision.

## Development and rule selection

Score each **original frozen estimate**, including estimates subsequently withheld, using the stress design's error measure. For true limits $L_0,U_0$ and estimated limits $\hat L,\hat U$ on the measurement scale, the error is

$$
E=100\,\frac{\max(|\hat L-L_0|,|\hat U-U_0|)}{U_0-L_0}.
$$

For each candidate, pool the 900 development datasets. A candidate qualifies if it declines at most 5% of all estimates with $E<5$. Among qualifying candidates, select the one that declines the largest **number** of estimates with $E>10$. Break ties first by higher $t$, then by rule a. Use exact counts for the cap and the ranking. Estimates at a cutoff stay in all applicable summaries but do not enter a strict inequality at that cutoff. Truth enters only development selection and retrospective scoring. The fitter and the decision rules use only observed inputs. Simulated routine disease labels, true limits, true overlap, affected fractions and error are available only for study scoring.

**Development choice to commit before confirmation**

| Rule | Threshold | Declined / all with $E>10$ | Declined / all with $E<5$ |
| --- | --- | --- | --- |
| Pending development | Pending | Pending | Pending |

Record all 14 development results. Enter the winner and its counts above and commit the protocol and analysis code before generating, fitting or inspecting any replicate from 21–40. If no candidate meets the 5% cap, record that result and stop before confirmation.

## Confirmation and interval reporting

After the selection commit, generate all 45 cases for replicates 21–40. Run both fits and their equivalence check on all 900 datasets. Apply the committed rule and the LRT-alone comparator to each dataset. Save both decisions before joining simulation truth for scoring. The rule, threshold, fitting code, analysis and success criterion are fixed from that commit onward.

When a rule reports an interval, return the frozen lower and upper limits. When it declines, leave **both reported limits missing** and retain the fitted limits in diagnostic fields only. Use these messages verbatim.

| Reason for decline | Message |
| --- | --- |
| Separation flag | routine results are too mixed to separate healthy from affected patients |
| LRT flag | routine data disagree with the screened sample |

If both flags trigger rule b, display the separation message, retain both flags and count one decline. If an unavailable LRT is the sole reason for withholding, record that the LRT could not be evaluated. For every withheld interval, recommend establishing an interval with at least 120 screened patients or verifying a suitable published interval with 20 screened patients.

## Report

Lead with confirmation and show development separately. For the selected rule and the LRT alone, report the following in each phase.

- For each family, give a 3 × 5 table with overlap in rows and affected fraction in columns. For each cell, show mean fitted $O$ beside the true overlap, the number with defined $O$ out of 20, and the decline count and rate for each rule. The fits and overlap values are common to both rules. Report one-component counts per cell and overall.
- Give the number declined divided by the total in each error group, $E>10$, $E>25$ and $E<5$, with percentages. The denominator includes every scored estimate in that group, including withheld estimates. Count technical withholdings in decline rates and identify their number separately.
- Among reported intervals, give the number reported, mean, median, 90th percentile (R quantile type 7) and maximum $E$, and the number and share with $E>10$. Show counts against all 900 attempts and identify any fit failures or unavailable statistics. Mark a zero denominator as undefined.
- For confirmation, give a single verdict for the committed rule using unrounded counts. A pass requires at least two-thirds declined among $E>10$ and at most 5% declined among $E<5$. Show both fractions. An empty required error group makes the verdict unevaluable. Report the outcome without adjusting the rule.
