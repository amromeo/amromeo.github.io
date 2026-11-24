# Improved anchoredRI comparator for anchorRefine

[Download the per-dataset CSV](anchoredRI-comparator.csv). Each row contains the improved anchoredRI baseline fit from GitHub Actions. The fitter selected the lowest-BIC eligible identity, log or shifted-log model with one to three components.

| Benchmark | Datasets |
| --- | ---: |
| RIbench N ≤ 5,000 | 2,880 |
| RIbench N = 50,000 | 1,440 |
| Simulation grid, replicates 1 to 20 | 900 |
| Total unique datasets | 5,220 |

The large-N run refitted all 1,440 datasets, including 144 previously exported Hb datasets with matching routine and anchor hashes. The small-N and grid rows are byte-identical to the accepted export. The two-setting grid uses keys of the form `s-CC-RR` (case and replicate) and reuses the screened sample across scenarios within each family and replicate. [Large-N results](../accuracy/large/results/REPORT.md) · [Small-N and grid results](../accuracy/results/REPORT.md)

## Pairing datasets

Use the unique `dataset_id`, which is `r-<Index>` for RIbench. Verify both input hashes. The 40 independent anchors were drawn with seed `17000000 + Index`. The RIbench routine seed repeats across datasets. The `routine_hash` and `anchor_hash` are SHA-256 digests of the serialized R input vectors using `digest::digest(x, algo="sha256")` under R 4.5.0. The original RIbench `HashCode` is also present as `ribench_hashcode`.

```r
ref <- read.csv("anchoredRI-comparator.csv", na.strings = "", stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(ref$dataset_id), !anyDuplicated(anchorRefine_results$dataset_id))
keys <- c("dataset_id", "routine_hash", "anchor_hash")
overlap <- intersect(anchorRefine_results$dataset_id, ref$dataset_id)
paired <- merge(anchorRefine_results, ref, by = keys, suffixes = c("_anchorRefine", "_anchoredRI"))
stopifnot(setequal(paired$dataset_id, overlap))
```

| Columns | Meaning |
| --- | --- |
| `schema_version`, `method`, `fitter`, `fitter_source_commit`, `max_components` | Comparator schema 1, method anchoredRI, round 3 multistart fitter, exact source revision and a three-component cap. |
| `benchmark`, `stage`, `Index`, `case`, `replicate` | Benchmark and dataset identifiers. The case and replicate fields identify simulation grid cells and replicates. |
| `Analyte`, `N`, `n_routine`, `n_anchor`, `distribution`, `pathological_fraction`, `scenario`, `rounding_step` | Study definition, nominal sample size, actual input counts and recording precision. Pathological fraction is the design fraction in each benchmark. |
| `transform`, `lambda`, `shift`, `components`, `loglik`, `bic`, `convergence_code`, `at_bound` | Selected candidate and fit diagnostics. |
| `lower_probability`, `upper_probability`, `lower`, `upper` | Requested percentiles and raw fitted limits. Use these for comparisons before the safeguard. |
| `routine_min`, `lrt_statistic`, `lrt_p`, `lower_flag`, `declined`, `technical`, `report_status`, `reason` | The original Task E safeguard and its result for this fit. |
| `reported_lower`, `reported_upper` | Limits released by the safeguard. Withheld limits are blank. A one-sided interval also has a blank lower endpoint. |

All finite numeric fields preserve double precision. The exporter rereads the CSV and validates limits, likelihoods and LRT values against the saved fits. [Small-N export](../accuracy/extract.R) · [Large-N export](../accuracy/large/extract.R)
