# Task E validation and reproducibility

All 5,220 same-runtime comparisons matched the frozen C2 source to an absolute tolerance of 1e-10, with the same selected transform and component count.

The archived comparison included 4,320 RIbench datasets. Sixty-two datasets (1.44%) selected a different transform or component count from the archived run. This passed the prespecified 5% review gate. Limits differed by more than 0.1% of the archived interval width in 447 datasets. The largest absolute difference was 4,635.05 for IgE Index 4450. The largest relative difference was 688.37 archived widths for the one-sided CRP Index 4926, whose archived zero-to-upper span was 0.01008.

The local baseline had a runner-up within two BIC points in 269 of 5,220 fits (5.15%). This included 208 of 2,880 small-N RIbench fits, 15 of 1,440 N = 50,000 fits and 46 of 900 simulation fits. Across those 269 fits, the median change in the more affected endpoint was 6.18% of the selected interval width.

The nine [candidate tables](baseline-candidates) contain all 46,980 baseline candidate rows in their original order, including each model’s BIC and limits. Each dataset appears in one file. The [BIC-gap table](anchoredRI-bic-gaps.csv.gz) gives the selected model and runner-up for every dataset. The [Hb 100 and Hb 148 comparison](../../validation/HB-BIC.md) gives their local and replay BICs. The original C2 archive retained only the selected model for each dataset.

The Arm A replay used the verified C2 GitHub Actions runtime and exactly RIbench Index 1 through 156. All source, routine-data and anchor hashes matched. Twenty-five datasets selected a different component count or shift specification across runtimes. Ninety-six limit pairs differed by more than 0.1% of the Actions interval width. The largest absolute endpoint difference was 1.75194 Hb units at Index 69. The largest difference as a fraction of the Actions interval width was 27.99% at Index 91.

The first replay attempts for five shards ran on a newer Actions image and were rejected by the runtime guard. Their retries passed. The original collection job selected the five small failed-attempt artifacts. A collection-only recovery selected the complete artifact for every shard and verified all 156 unique indices. No fit was repeated during collection.

The validation timing file covers the last resumed validation segment. Earlier segments also contributed checkpoints.

[C2-image replay run](https://github.com/amromeo/amromeo.github.io/actions/runs/45131636241) · [Completed replay artifact](https://github.com/amromeo/amromeo.github.io/actions/runs/45135677194)
