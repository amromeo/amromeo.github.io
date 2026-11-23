# Nesting audit of the saved fits

The frozen baseline had at least one likelihood inversion in 299 of 5,220 datasets (5.73%). A likelihood inversion means that the saved fit for a larger model had a lower original-scale log likelihood than a smaller model it should contain. BIC can legitimately favor the smaller model because it penalizes additional parameters. This audit tests the likelihood before that penalty.

| Saved fits | Datasets | With an inversion | Winner changed using known nested likelihoods |
|---|---:|---:|---:|
| Baseline, RIbench N ≤ 5,000 | 2,880 | 154 (5.35%) | 0 |
| Baseline, RIbench N = 50,000 | 1,440 | 73 (5.07%) | 0 |
| Baseline, simulation grid | 900 | 72 (8.00%) | 0 |
| Arm A, first 156 Hb datasets, local | 156 | 31 (19.87%) | 0 |
| Arm A, same 156 Hb datasets, Actions | 156 | 31 (19.87%) | 0 |
| Arm A, other saved small-N fits | 240 | 66 (27.50%) | 4 |
| Arm A, other saved N = 50,000 fits | 144 | 76 (52.78%) | 11 |
| Arm B, all saved fits before the hold | 540 | 333 (61.67%) | 0 |

For the last column, each larger candidate receives at least the best known nested likelihood, including bounds propagated through adjacent component counts, and keeps its own parameter penalty. The winner changed in 15 of the 384 additional local Arm A fits. These are lower-bound calculations. The actual BIC winner after a correct fit cannot be recovered from the saved candidates. The improved fitter provides an observed new choice on the 156 paired datasets.

The comparisons use statistical nesting, the allowed shift interval, the counts of retained observations and a floating-point allowance. Some historical numerical bounds on SDs and weight ratios differ between candidate optimizations. A flagged inversion identifies a problem with the attained fit or those bounds. The original RDS files do not contain a feasible parameter witness for every losing candidate. A shift-versus-zero comparison is excluded when zero falls outside its search interval. Identity can be embedded at lambda one by translating component means.

The Actions Arm A checkpoints retained the fixed-model starting parameters. We evaluated the likelihood at those saved parameters on the same verified routine and anchor vectors. The Actions fixed-model seed inventory is incomplete for 29 of 156 datasets. We count only comparisons with both likelihoods available when estimating inversion rates. The original C2 archive kept the selected model. A separate archived-image preflight saved all candidates for Hb 100 and 148.

[Dataset counts](violation-summary.csv) · [Counts by nesting relationship](by-comparison.csv) · [Dataset flags and BIC selections](datasets-audited.csv.gz)

The [pair-level tables](pairs) give both likelihoods, the signed difference, the floating-point allowance and the status of every comparison.
