# Component 1 routine shares

Each fit is paired with the unweighted full-data fit for the same dataset and component cap. The table uses the designed physiological fraction, which excludes the 0.1% uniform background. Positive gap means the model assigned too much routine probability to component 1. A negative change in absolute gap or endpoint error favors the candidate.

| Search | Estimate | Seed set | Definition | Datasets assigned | Mean gap (pp) | Mean change in absolute gap (pp) | Mean change in endpoint error |
| --- | --- | --- | --- | ---: | ---: | ---: | ---: |
| 3 | W | — | all | 1,440 | -1.2304 | -3.3132 | -0.3380 |
| 5 | W | — | all | 1,440 | -9.3619 | -2.5985 | -0.2123 |
| 3 | full | — | all | 1,440 | -4.5718 | 0.0000 | 0.0000 |
| 5 | full | — | all | 1,440 | -13.1369 | 0.0000 | 0.0000 |
| 3 | median | A | all | 1,440 | -0.8293 | -4.9622 | -0.4302 |
| 3 | median | A | passing | 1,440 | -1.2080 | -6.3594 | -0.7356 |
| 3 | median | B | all | 1,440 | -0.7899 | -4.9209 | -0.4569 |
| 3 | median | B | passing | 1,440 | -1.3734 | -6.1303 | -0.6893 |
| 5 | median | A | all | 1,440 | -8.1187 | -3.8298 | -0.2366 |
| 5 | median | A | passing | 1,440 | -8.7815 | -3.1438 | -0.2397 |
| 5 | median | B | all | 1,440 | -8.2065 | -3.7650 | -0.2373 |
| 5 | median | B | passing | 1,440 | -8.8801 | -3.0902 | -0.2405 |

At cap 5, the mean absolute share gap was 14.95 percentage points for the full fit, 12.35 for W and 11.12 for the primary median. At cap 3, it was 13.50, 10.19 and 8.54, respectively. The cap-3 median came closer to the designed share, while its pre-check endpoint error was higher than the cap-5 median’s (0.421 versus 0.174). The table counts assigned datasets. For the passing-only medians, the comparable counts for seed sets A and B were 1,292 and 1,295 at cap 3, and 1,430 and 1,433 at cap 5.

The signed and absolute Spearman correlations use the share and endpoint errors of the same usable fits. For a median, the share and limits used for these correlations are medians over the same contributing draws. A blank correlation indicates fewer than two distinct values on at least one axis. Each draw contributes once within its seed set. Each median contributes once per dataset, cap, bank and definition. Individual models are grouped by their selected component count and medians by their full fit’s count. The released count for a median uses its majority rule. [Component groups](by-components.csv) · [Groups by pathological fraction](by-components-and-fraction.csv) · [Absolute gap bands](gap-bands.csv)

0 datasets had a discrepancy between the reconstructed component counts and observed routine length. The realized-fraction sensitivity uses the full routine vector and is reported separately from the primary design-fraction analyses. [Count discrepancies](count-discrepancies.csv)

The per-fit export includes requested endpoints, safeguard flags, contribution counts, component counts across draws, the parent full fit’s gap and all paired changes. Passing-only individual-model summaries select the draws that passed the first-stage checks. [Paired per-fit and per-median records](paired-fits.csv.gz)

Share versus error figures use the component count of the individual fit, or the full fit for a median. Columns show the two seed sets separately. [Cap 3, all](share-gap-cap3-all.png) · [Cap 3, passing](share-gap-cap3-passing.png) · [Cap 5, all](share-gap-cap5-all.png) · [Cap 5, passing](share-gap-cap5-passing.png)

[Prior Task E archived component-share report](https://github.com/amromeo/amromeo.github.io/blob/e5a4577c1fe73869c6a24ce64c8759f190de3197/anchoredRI/tools/subsample/archived-shares/REPORT.md)
