# Task E at 50,000 routine results

Improved baseline missed the score criterion and the implausible-fit criterion. Arm B missed the score criterion.

| Arm | Native score before | refineR threshold | Score criterion | Implausible or failed / 1,440 | Count criterion (≤14) | Overall |
| --- | --- | --- | --- | --- | --- | --- |
| Improved baseline | 0.36176078 | 0.22818185 | Miss | 37 | Miss | Miss |
| Arm B | 0.37436430 | 0.22818185 | Miss | 1 | Pass | Miss |

Each arm had 1,440 attempts. RIbench omits a dataset from its native cutoff mean when its mean absolute endpoint z deviation exceeds 5. The primary comparison uses raw fitted limits before either safeguard. Technical failures and implausible fits count once toward the limit of 14. The full-precision threshold and the prespecified additional score for missing fits are in [the criteria table](large-success-criteria.csv).

## Accuracy and coverage

| Method | Native before | Native after | In cutoff before | In cutoff after | Uncapped before | Uncapped after | Implausible before | Implausible after | Reported / attempted | Declined | Technical |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Improved baseline | 0.3618 | 0.1706 | 1403 | 1066 | 0.9254 | 0.1706 | 37 | 0 | 1066/1440 | 374 (26.0%) | 0 |
| Arm B | 0.3744 | 0.1307 | 1439 | 1135 | 0.4026 | 0.1307 | 1 | 0 | 1135/1440 | 305 (21.2%) | 0 |
| Original frozen fit, same Actions runner | 0.3689 | 0.1742 | 1400 | 1041 | 0.9143 | 0.1742 | 40 | 0 | 1041/1440 | 399 (27.7%) | 0 |
| Original archived C2 | 0.3656 | 0.1735 | 1400 | 1038 | 0.9572 | 0.1735 | 40 | 0 | 1038/1440 | 402 (27.9%) | 0 |
| Anchors only | 0.2961 | 0.2961 | 1440 | 1440 | 0.2961 | 0.2961 | 0 | 0 | 1440/1440 | 0 (0.0%) | 0 |
| refineR | 0.2282 | 0.2282 | 1438 | 1438 | 0.3099 | 0.3099 | 2 | 2 | 1440/1440 | 0 (0.0%) | 0 |

After-safeguard scores describe the intervals each method released. RIbench uses two endpoints for nine analytes and the 95th percentile alone for CRP. The LRT declines at p < 0.01. The lower-limit check declines two-sided intervals with a lower limit below every routine observation. Heavily mixed populations can still yield inaccurate intervals that pass both checks.

### Skewed FT4 datasets with lambda = 0.432

| Method | Native before | Native after | In cutoff before | In cutoff after | Uncapped before | Uncapped after | Implausible before | Implausible after | Reported / attempted | Declined | Technical |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Improved baseline | 0.2528 | 0.1708 | 144 | 127 | 0.2528 | 0.1708 | 0 | 0 | 127/144 | 17 (11.8%) | 0 |
| Arm B | 0.1613 | 0.1596 | 144 | 141 | 0.1613 | 0.1596 | 0 | 0 | 141/144 | 3 (2.1%) | 0 |
| Original frozen fit, same Actions runner | 0.2567 | 0.1706 | 142 | 125 | 0.5149 | 0.1706 | 2 | 0 | 125/144 | 19 (13.2%) | 0 |
| Original archived C2 | 0.2575 | 0.1713 | 142 | 125 | 0.5154 | 0.1713 | 2 | 0 | 125/144 | 19 (13.2%) | 0 |
| Anchors only | 0.2456 | 0.2456 | 144 | 144 | 0.2456 | 0.2456 | 0 | 0 | 144/144 | 0 (0.0%) | 0 |
| refineR | 0.1792 | 0.1792 | 144 | 144 | 0.1792 | 0.1792 | 0 | 0 | 144/144 | 0 (0.0%) | 0 |

## Component selection in Arm B

| Components | Share of all attempts | Native before | Uncapped before | Implausible before | Declines among selected | Native after | Reported |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 to 3 | 140/1440 (9.7%) | 0.0422 | 0.0422 | 0 | 3/140 (2.1%) | 0.0429 | 137 |
| 4 | 359/1440 (24.9%) | 0.1125 | 0.1125 | 0 | 13/359 (3.6%) | 0.0939 | 346 |
| 5 | 941/1440 (65.3%) | 0.5238 | 0.5668 | 1 | 289/941 (30.7%) | 0.1687 | 652 |
| 4 or 5 | 1300/1440 (90.3%) | 0.4102 | 0.4414 | 1 | 302/1300 (23.2%) | 0.1427 | 998 |

The improved baseline allows up to three components. The pooled four-or-five row includes both individual component groups. Shares use all 1,440 attempts. Decline rates use the number selecting that component count. [Individual component counts and errors for both methods](component-counts-and-errors.csv)

## Paired comparisons

| Arm | Comparator | Filter | Complete pairs | Mean change in uncapped absolute z | Wins / ties / losses |
| --- | --- | --- | --- | --- | --- |
| Improved baseline | Arm B | before | 1440 | 0.5228 | 534 / 140 / 766 |
| Improved baseline | Arm B | after | 946 | 0.0458 | 290 / 137 / 519 |
| Improved baseline | Original frozen fit, same Actions runner | before | 1440 | 0.0111 | 766 / 0 / 674 |
| Improved baseline | Original frozen fit, same Actions runner | after | 999 | -0.0012 | 527 / 0 / 472 |
| Improved baseline | Original archived C2 | before | 1440 | -0.0319 | 764 / 0 / 676 |
| Improved baseline | Original archived C2 | after | 998 | -0.0016 | 531 / 0 / 467 |
| Improved baseline | Anchors only | before | 1440 | 0.6292 | 867 / 0 / 573 |
| Improved baseline | Anchors only | after | 1066 | -0.0826 | 732 / 0 / 334 |
| Improved baseline | refineR | before | 1440 | 0.6154 | 643 / 0 / 797 |
| Improved baseline | refineR | after | 1066 | -0.1512 | 595 / 0 / 471 |
| Arm B | Improved baseline | before | 1440 | -0.5228 | 766 / 140 / 534 |
| Arm B | Improved baseline | after | 946 | -0.0458 | 519 / 137 / 290 |
| Arm B | Original frozen fit, same Actions runner | before | 1440 | -0.5117 | 846 / 0 / 594 |
| Arm B | Original frozen fit, same Actions runner | after | 932 | -0.0492 | 574 / 0 / 358 |
| Arm B | Original archived C2 | before | 1440 | -0.5547 | 850 / 0 / 590 |
| Arm B | Original archived C2 | after | 929 | -0.0491 | 580 / 0 / 349 |
| Arm B | Anchors only | before | 1440 | 0.1064 | 953 / 0 / 487 |
| Arm B | Anchors only | after | 1135 | -0.1170 | 866 / 0 / 269 |
| Arm B | refineR | before | 1440 | 0.0926 | 825 / 0 / 615 |
| Arm B | refineR | after | 1135 | -0.0928 | 759 / 0 / 376 |

Negative differences favor the arm. After safeguarding, each pair includes only datasets for which both methods reported an interval. The paired tables also give each distribution by pathological fraction at or below 20% versus above 20%. [Paired results](ribench-pairs.csv) · [Subgroup scores and denominators](ribench-summary.csv) · [Safeguard flags by error](ribench-safeguard-by-error.csv)

## Fit diagnostics and runtime

| Method | Finite runner-up | Within 2 BIC | Median endpoint movement in near ties (% of width) | Iteration-limit wins | Line-search-error wins |
| --- | --- | --- | --- | --- | --- |
| Improved baseline | 1440 | 14/1440 (1.0%) | 0.7% | 8 | 0 |
| Arm B | 1440 | 48/1440 (3.3%) | 0.8% | 454 | 0 |

For one-sided CRP, the reference width for a near tie runs from zero to the fitted upper limit. The recorded exits describe the optimizer for the selected solution. [Candidate gaps](runner-up-gaps.csv) · [Winning optimizer exits](winning-optimizer-exits.csv)

GitHub Actions fitted all 1,440 datasets in 24 batches, with four workers per batch. Input preparation took 2.8 minutes. The fitting window was 32.7 minutes and the two arms used 2.04 CPU hours. Extraction, scoring and report generation took 0.4 minutes. Runner image revisions and counts were 20261116.133.1 (1080), 20261116.376.1 (360). [Actions run](https://github.com/amromeo/amromeo.github.io/actions/runs/45298562883) · [Timing data](runtime.csv)

RIbench supplied a single pooled routine sample and 40 independent screened anchors per dataset. Each routine vector and anchor vector matched its archived hash. The unchanged round 3 fitter considered identity, log when eligible, and shifted log. Baseline and Arm B used the same first three component candidates. The original frozen estimator and its one-group study adapter were also fitted on each Actions runner. The archived C2, anchors-only and refineR results came from [the independent-anchor C2 run](https://github.com/amromeo/amromeo.github.io/actions/runs/45101564181). The previously accepted two-setting simulation grid was evaluated separately.

The comparator contains 5,220 unique datasets, including all 1,440 large-N baseline fits. The Actions estimates replace 144 previously exported large-N rows after both input hashes match. The 2,880 small-N rows and 900 grid rows are preserved byte for byte. [Comparator CSV](anchoredRI-comparator.csv)
