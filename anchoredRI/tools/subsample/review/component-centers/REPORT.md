# Where the extra components are centered in the archived Task E Arm B fits

Among the 941 full-data Task E Arm B fits that selected five components at N = 50,000, **436 (46.3%) placed at least one extra component in the central 50% of the true healthy distribution**. On those 436 datasets, the qualifying components accounted for 35.25% of the routine mixture on average. Component 1 fell short of the true healthy share by 37.80 percentage points on average, compared with 3.36 on the 505 datasets without a centrally located extra component. The corresponding mean endpoint errors were 1.036 and 0.162. All errors and shares come from these same Task E fits, including their first-stage declines. [Per-dataset results](datasets.csv)

| Location of extra component center | Fits with at least one | Qualifying components out of 3,764 | Mean combined routine weight, all 941 fits | Mean combined routine weight, fits with at least one |
|---|---:|---:|---:|---:|
| Within the true reference limits | 659 (70.0%) | 1,186 | 25.83% | 36.88% |
| Within the true healthy interquartile interval | 436 (46.3%) | 537 | 16.33% | 35.25% |

The two rows overlap. Within each row, weights sum over all qualifying components in a fit. A fit with no qualifying center contributes zero to the all-fit mean. The weights describe the entire routine mixture mass of those components, which also extends beyond the specified region. [Centers, weights and healthy percentiles for all four extra components](components.csv)

| Center region | Mean component-1 shortfall, any / none (percentage points) | Mean uncapped endpoint error d, any / none | Spearman ρ, qualifying weight vs. shortfall | Spearman ρ, qualifying weight vs. d |
|---|---:|---:|---:|---:|
| True reference interval | 27.78 / −0.47 | 0.767 / 0.099 | 0.853 | 0.773 |
| True central 50% | 37.80 / 3.36 | 1.036 / 0.162 | 0.820 | 0.665 |

Of the 436 fits with a center in the true central 50%, 337 had one, 97 had two and two had three. Greater combined weight in these components tracked larger component-1 shortfalls and larger endpoint errors within each of the eight pathological-fraction groups. The associations varied by analyte. For example, in LDH the central weight correlated positively with shortfall (ρ = 0.63) and negatively with endpoint error (ρ = −0.24). [Summary by pathological fraction, analyte and first-stage decision](summary.csv) · [Counts and errors by number of qualifying centers](by-count.csv)

Even at 0% designed pathology, 29 of 46 five-component fits placed an extra center in the true central 50%. Those centers carried a combined 49.65% mean routine weight on those 29 datasets. The generating design comprised 99.9% healthy observations and 0.1% uniform background. This supports the hypothesis that the mixture can describe a single healthy distribution with several components. With disease present, an overlapping pathological component could also have a center in this region. Component weights sum to one, which itself links a large extra weight with a smaller weight for component 1. These comparisons cannot assign individual healthy observations to a particular fitted component or determine the cause of a specific error.

## Cohort and calculation

The archived Task E Arm B run selected five components on 941 of its 1,440 datasets. The later Task H run selected five components on 934 datasets. This analysis uses Task E's fitted parameters and matches its own saved endpoints and errors by original RIbench Index. Component 1 is the component anchored by the screened patients. The other four are unscreened mixture components. For each, the center is the Gaussian location in its selected transform, mapped back to measurement units. For identity fits the Task E export already puts its locations in measurement units. For log and shifted-log fits, the inverse-transformed location is the component’s original-scale median, given by the shift plus the exponential of the fitted location. Every component-1 inverse transform reproduced its original raw requested limits. [Task E parameter exporter](https://github.com/amromeo/amromeo.github.io/blob/8bcd49c3c7977c50982524b9dbc48d16bdfa5914/anchoredRI/tools/boxcox/accuracy/large/extract.R#L60-L71) · [Selected parameters](https://github.com/amromeo/amromeo.github.io/blob/8bcd49c3c7977c50982524b9dbc48d16bdfa5914/anchoredRI/tools/boxcox/accuracy/large/results/selected-component-parameters.csv) · [Task E share report](../../archived-shares/REPORT.md)

The reference interval uses the published RIbench limits, inclusive. For CRP it includes values at or below its upper reference limit. Each healthy interquartile interval includes both bounds, even for CRP. We classify each center by its percentile under the continuous generating healthy distribution using the Box–Cox parameters, shift and support restriction of RIbench 1.0.2. The component-1 shortfall is the true designed healthy share (one minus the pathological fraction and the 0.1% background) minus the fitted routine share, expressed in percentage points. A positive value indicates underestimation. The error d is the uncapped mean absolute pre-safeguard Task E z deviation of the requested limits, using only the upper limit for CRP. Each fit has equal weight. Correlations use Spearman ranks with averaged ties, including fits with zero qualifying weight. [RIbench truth](https://github.com/amromeo/anchorRefine/blob/638e2d3500faa01d4e2816895a95d7a6685567a0/stage1/large/testsets.csv) · [Task E fits](https://github.com/amromeo/amromeo.github.io/blob/8bcd49c3c7977c50982524b9dbc48d16bdfa5914/anchoredRI/tools/boxcox/accuracy/large/results/armB-fits.csv) · [Task E error records](https://github.com/amromeo/amromeo.github.io/blob/8bcd49c3c7977c50982524b9dbc48d16bdfa5914/anchoredRI/tools/boxcox/accuracy/large/results/ribench-outcomes.csv)

## Reproduction and checks

Download the four pinned CSV files linked above into the same directory. From this report’s directory run the commands below with Python 3.10 or later and the Task H R environment (RIbench 1.0.2). The Python script checks the original input SHA-256 hashes and writes the four output tables. The independent R calculation maps the true quartiles back to the original scale with RIbench’s inverse Box–Cox function. It reproduced both classifications for all 3,764 centers, both qualifying weight sums on all 941 datasets, the shortfalls and raw endpoint errors, and the four pooled Spearman correlations. [Analysis](analyze.py) · [Independent verification](verify.R) · [Original file hashes](../../archived-shares/sources.json) · [Report checksums](SHA256SUMS)

```sh
python3 analyze.py /path/to/archived-inputs
Rscript verify.R /path/to/archived-inputs .
```
