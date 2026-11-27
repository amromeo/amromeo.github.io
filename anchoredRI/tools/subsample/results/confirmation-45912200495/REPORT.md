# Task H independent confirmation at N = 50,000

**The subsampling policy failed independent confirmation.** The policy and the package returned intervals on exactly the same 1,437 of 1,440 datasets. Mean error on those returns rose from 0.0983 with the package to 0.1222 with the primary subsample seeds. The cell-resampled 95% interval for the increase was +0.0148 to +0.0338. The upper bound exceeded the prespecified +0.02 margin. The interval includes +0.02, so the size of the increase relative to that margin remains uncertain. Seed set B on the same datasets and anchors gave the same decision. The RIbench accuracy gain did not reproduce on these independent synthetic cases.

| Method | Returns | Mean error on returns | All-assignment score | Intervals with an extreme endpoint |
| --- | --- | --- | --- | --- |
| Package | 1437 | 0.0983 | 0.1085 | 0 |
| A through step c | 1437 | 0.1222 | 0.1324 | 0 |
| A through step b | 1437 | 0.1222 | 0.1324 | 0 |
| B through step c | 1437 | 0.1221 | 0.1323 | 0 |
| B through step b | 1437 | 0.1221 | 0.1323 | 0 |

Errors use the two true healthy endpoints for each dataset. The all-assignment score caps the mean absolute endpoint error at 5 and assigns 5 to every nonreturn. An extreme interval has at least one endpoint with absolute healthy z error above 5. [All methods, breakdowns and bootstrap denominators](metrics.csv)

## Paired accuracy and decisions

| Seed set | Common returns | Policy mean d | Package mean d | Uncapped difference [95% interval] | Capped difference [95% interval] | Better / tied / worse |
| --- | --- | --- | --- | --- | --- | --- |
| A | 1437 | 0.1222 | 0.0983 | 0.0239 [0.0148, 0.0338] | 0.0239 [0.0148, 0.0338] | 509 / 33 / 895 |
| B | 1437 | 0.1221 | 0.0983 | 0.0238 [0.0146, 0.0338] | 0.0238 [0.0146, 0.0338] | 508 / 33 / 896 |

| Seed set | Lower observed mean error | Improved accuracy | No worse | Accuracy gate | Coverage and extremes | All-assignment gate | Step c | Four-step judgment |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A | no | no | no | fail | pass | fail | insufficient evidence | fail |
| B | no | no | no | fail | pass | fail | insufficient evidence | fail |

The paired differences subtract the package error. Improved accuracy requires the entire two-sided cell-resampled 95% interval for the uncapped difference below zero. No worse requires the upper 95% bound for the capped difference below +0.02. The accuracy gate also checks the capped point estimates in both pathology strata. The bootstrap resampled factorial cells 2,000 times and kept the eight replicates in each sampled cell together. [Decision table](gates.csv) · [Protocol](https://github.com/amromeo/amromeo.github.io/blob/140bd43/anchoredRI/tools/subsample/confirmation/PROTOCOL.md)

| Set A stratum | Common n | Capped change | Policy extremes | Package extremes | Lost package returns |
| --- | --- | --- | --- | --- | --- |
| <=20% | 576 | 0.0075 | 0 | 0 | 0 |
| >20% | 861 | 0.0349 | 0 | 0 | 0 |

Across all assignments the score difference was 0.0239 [0.0147, 0.0337] for A, and 0.0238 [0.0146, 0.0338] for B. All bounds and gate decisions use full precision. Signed infinite differences enter the bootstrap, and the tables count undefined replicates. [Common-return, finite-pair and all-assignment statistics](metrics.csv)

## How the policy routed datasets

| Step | Output | A assignments | B assignments | A returns | B returns |
| --- | --- | --- | --- | --- | --- |
| a | Median after full-data pass | 1404 | 1404 | 1404 | 1404 |
| b | Package fallback | 33 | 33 | 33 | 33 |
| c | Median after package nonreturn | 0 | 0 | 0 | 0 |
| d | No interval | 3 | 3 | 0 | 0 |

| Set | Package nonreturns | Step c returns | Mean d | Median d | Maximum d | Extreme intervals | Extreme rate, exact 95% interval | Judgment |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A | 3 | 0 | undefined | undefined | undefined | 0 | [undefined, undefined] | insufficient evidence |
| B | 3 | 0 | undefined | undefined | undefined | 0 | [undefined, undefined] | insufficient evidence |

With fewer than 10 recoveries the step-c evidence is insufficient. At 10 or more, success requires mean uncapped error below 0.5 and zero extreme endpoints. The exact binomial interval describes uncertainty among the observed recoveries in this factorial design. Each package nonreturn’s flags, votes, contributing draws and reasons are included. [Each package nonreturn](package-nonreturns.csv.gz) · [Reasons and majority outcomes](nonreturn-reasons.csv) · [Extreme endpoints in all methods](extreme-endpoints.csv) · [Step-c summaries](recoveries.csv) · [Routes by factor](routes.csv)

The only three package nonreturns occurred with shoulder disease. For all three the first-stage anchor check and fallback check declined, and 40–47 of the 50 subsamples declined across the two banks. Since neither bank recovered a case, the rescue-error mean and the exact binomial interval for the extreme rate among rescues have no denominator. All paired accuracy changes came from the 1,404 step-a replacements. The 33 step-b intervals copied the package exactly. [Recovery status by dataset](seed-comparison.csv)

## Healthy family

| Set A family | Assigned | Common returns | Uncapped change [95% interval] | All-assignment change | Additional returns |
| --- | --- | --- | --- | --- | --- |
| linear | 480 | 480 | 0.0340 [0.0195, 0.0497] | 0.0340 | 0 |
| lognormal | 480 | 479 | 0.0007 [-0.0085, 0.0100] | 0.0007 | 0 |
| power_0.5 | 480 | 478 | 0.0371 [0.0168, 0.0568] | 0.0369 | 0 |

The power-0.5 generating family has no exact match among anchoredRI’s fitting transforms. It showed higher error under subsampling, as did the linear family. The lognormal family showed little change. Seed set B gave the same pattern. These family comparisons are descriptive. [Exact fraction, disease pattern, rounding and family tables](metrics.csv) · [Per-case decisions, endpoints and scores](cases.csv.gz) · [Readback and independent score verification](verification.json)

## Separate RIbench component analysis

The archived Task E Arm B analysis used its own 941 full-data fits that selected five components. Of these, 659 (70.0%) had at least one extra component centered inside the true reference limits, and 436 (46.3%) had one centered in the true healthy interquartile interval. When present, these components carried 36.9% and 35.2% of the routine weight on average, respectively. With a center in the healthy interquartile interval, component 1’s mean shortfall was 37.80 percentage points and mean endpoint error was 1.036. Without one, the corresponding values were 3.36 points and 0.162. Even with zero designed pathology, 29 of 46 five-component fits had a central extra component. This supports splitting of the healthy distribution as one possible explanation. Disease can also overlap the center. The synthetic confirmation tested the release policy independently of these descriptive fits. [Task E component analysis, weights and per-dataset evidence](../../review/component-centers/REPORT.md)

## Verification and data

The approved protocol preceded the source freeze and launch. All 192 GitHub Actions jobs succeeded on their first attempt. Both fixture hosts produced identical decisions and sample positions, with zero endpoint and p-value differences. Independent readback covered all 145,440 fit records, 720 million subsample positions and 2,880 policy routes. The cohort had zero technical first-stage fit failures. A second implementation rebuilt the routes and scores from the original fit rows and compared 164 results, including all three accuracy statements, both banks’ rescue decisions and the whole-cell bootstrap bounds. All matched. [Protocol](../../confirmation/PROTOCOL.md) · [Source hashes](../../confirmation/SOURCE.sha256) · [Launch receipt](../../confirmation/launch.json) · [Workflow](https://github.com/amromeo/amromeo.github.io/actions/runs/45912200495) · [Verification counts](verification/summary.json) · [Independent calculation](verification/independent-gates.py) · [Independent check results](verification/independent-gates.json)

The permanent release contains eight archives of per-case records and one archive with the original input references and fixture results. Each archive has its own SHA-256 sidecar. [Data release](https://github.com/amromeo/amromeo.github.io/releases/tag/task-h-confirmation-50000) · [Archive membership and digests](release-assets.csv) · [Release checksums](RELEASE-SHA256SUMS) · [Report and table checksums](SHA256SUMS)
