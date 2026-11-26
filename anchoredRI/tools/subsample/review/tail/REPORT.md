# Where the five-component full fit underestimated healthy share

Across the 934 full-data fits that selected five components, the true healthy group made up 69.25% of the routine data on average. Component 1 accounted for 50.61%. Of the 18.64 percentage-point net deficit, **18.30 points (98.16%) fell inside the true reference limits** and **0.34 points (1.84%) outside**. Even allowing for the different amounts of healthy mass in the two regions, the deficit was larger inside. Component 1 fell short of the true mass by 27.81% inside and 9.89% outside.

| Region of measurement | True healthy mass, % of routine observations | Component 1 mass, % of routine observations | Deficit, percentage points |
| --- | ---: | ---: | ---: |
| Inside the true limits | 65.790 | 47.494 | 18.296 |
| Outside the true limits | 3.462 | 3.120 | 0.343 |
| Below the lower limit | 1.447 | 1.597 | −0.150 |
| Above the upper limit | 2.015 | 1.523 | 0.492 |
| All values | 69.252 | 50.613 | 18.639 |

The lower and upper rows partition the outside row. Component 1 had excess mass below the lower limit on average. The deficit above the upper limit exceeded that excess. These comparisons localize the *net* shortfall in the model’s distribution. They cannot identify which individual healthy observations contributed to other fitted components. A model could miss some healthy observations while assigning other observations in the same region to component 1.

Among the 701 fits where component 1 underestimated the healthy share overall, 25.34 of the 26.17 missing percentage points (96.84%) fell inside the true limits. The other 233 fits overestimated the healthy share by 4.01 points on average. Among the 633 five-component fits that passed the safeguards, the deficit averaged 8.77 points inside and 0.24 outside. Among the 301 that declined, it averaged 38.32 inside and 0.56 outside. [Summaries by analyte, first-stage decision and sign of the deficit](summary.csv) · [Each dataset and rounding sensitivity](per-dataset.csv)

## Calculation

Each dataset contributes equally, before the safeguards. For a region of values A, its deficit in percentage points of the **routine population** is

$$
100\{h\,P_H(X\in A)-w_1\,P_1(X\in A)\}.
$$

The designed healthy share is $h=1-\mathrm{fractionPathol}-\mathrm{bg\_fraction}$, with 0.1% uniform background. The weight $w_1$ and component-1 Gaussian location, spread and transform come from the unweighted, cap-5 full fit. $P_1$ uses the archived identity, log or shifted-log transform with its saved support. $P_H$ uses the RIbench 1.0.2 physiological distribution, a Gaussian on the generating Box–Cox scale with the specified shift, retaining nonnegative values as in the generator. The calculation uses the numeric ground-truth limits from each design record. For one-sided CRP, all values at or below the true upper limit are inside and the lower tail has zero mass. [RIbench truth table](https://github.com/amromeo/anchorRefine/blob/638e2d3500faa01d4e2816895a95d7a6685567a0/stage1/large/testsets.csv) · [RIbench 1.0.2 generator and transforms, in the pinned source archive](https://github.com/amromeo/anchorRefine/blob/638e2d3500faa01d4e2816895a95d7a6685567a0/vendor/RIbench_1.0.2.tar.gz) · [Fitted transforms and quantile calculation](https://github.com/amromeo/amromeo.github.io/blob/fe8105ac6961f1b45bbfbe034648938f227c4726/anchoredRI/R/safeguard.R)

The primary analysis integrates continuous distributions up to the published limits. The rounding sensitivity counts a value as inside when its rounded measurement is within those limits, inclusive. It uses the outer bin edges for the smallest and largest qualifying grid values. Under this definition, the deficit was 18.304 points inside and 0.334 outside.

All regional deficits retain their sign and sum to the dataset’s total healthy-share deficit. Reported percentages of the overall deficit divide summed regional deficits by the summed overall deficit. The corresponding percentages of true healthy mass divide by summed healthy mass in that region. [Calculation](decompose.py)

## Reproduction

Use Python 3.10 or later. Give the script the verified Task H release directory with one `<Index>/fits.csv` per dataset, and the pinned `testsets.csv`. It validates the truth file’s SHA-256, the original release manifest and every per-dataset fit table. On every selected five-component fit it checks that the saved component-1 parameters reproduce the raw requested endpoint probabilities, then checks that the regional masses add up. [Task H data release](https://github.com/amromeo/amromeo.github.io/releases/tag/task-h-subsample-50000) · [Per-file verification manifest](../../results/full-45860591866/verification/file-hashes.csv)

```sh
python decompose.py /path/to/verified/records /path/to/testsets.csv --out /path/to/output
```
