# Box-Cox shape and mixture capacity study

The accuracy comparison on 2,880 small RIbench and 900 grid datasets is complete for the improved baseline and Arm B. Arm A stopped during reliability testing. N = 50,000 awaits the separate runtime review. [Accuracy results](accuracy/results/REPORT.md) · [Protocol](PROTOCOL.md) · [Baseline comparator for anchorRefine](anchorRefine-comparator/README.md)

Implementation of the approved Task E protocol. The package code in `R/` is unchanged. The original frozen baseline is `55eee86b3f4111197b0043cc062ef0b781dc2057` and the one-group adapter is copied from `f72ea14fc7978bf2478f54c87548cb6c127fa766`.

Run from the package directory with RIbench 1.0.2 and the original study's R dependencies. `frozen/` contains verbatim sources. `capacity-fixed.R` and `capacity-shifted.R` add diagnostic returns and extend the existing initialization patterns to four and five components. Restricting them to three components must reproduce the frozen estimator.

`models.R` implements the shared Box-Cox likelihood, analytic gradients, original-scale Jacobian and BIC, fixed-shape checks, and study-local safeguards. Every candidate retains its starts, optimizer results and bound hits. For numerical stability each shape optimization uses fixed centering and scaling from its starting transform. SD bounds are .005 to 5 on that numerical scale, as in the shifted fitter. Mean bounds cover the observed range at the endpoints of the shape and offset search regions. The parameter count excludes these fixed numerical coordinates.

Fixed lambda 1 uses the identity optimizer with the required translation of means. Fixed lambda 0 uses the log optimizer. An independently written Box-Cox likelihood verifies every candidate's original-scale log likelihood. The tests also check the analytic gradient against finite differences, multiplicities and inverse-transform domains.

The configuration RDS contains `baseline`, `output`, `workers`, `small_archive`, `large_archive`, `small_data`, and `large_data`. The data directories contain `Data/<Analyte>/*.csv`. Archives are the completed C2 small-N and N=50,000 artifacts. The simulation comparators come from `assets/anchors-only/anchors-only-results.zip` and `assets/comparison-study/summary/comparison-results.zip` at the pinned commit. Each contains 900 comparator fits on the same case, replicate and routine seed.

```
Rscript tools/boxcox/test-models.R
Rscript tools/boxcox/run-study.R config.rds validate
Rscript tools/boxcox/run-study.R config.rds fit
Rscript tools/boxcox/score-study.R config.rds
```

Validation checks all 4,320 independent anchor hashes before fitting. Every routine vector must match RIbench's original hash. All 5,220 datasets must pass candidate and selected-fit equivalence before experimental fitting begins. The same-machine check has absolute tolerance 1e-10 and requires identical selected transform and component count. The full cross-machine audit reports limit differences above 0.1% of archived width and all model changes. The study pauses after all 4,320 RIbench comparisons if more than 5% select a different model. A same-machine failure stops validation immediately. The baseline in all comparisons is the local frozen refit. See the final protocol amendment for the denominator and one-sided convention.

A mismatch writes `HALT.txt`. Do not remove this file or loosen tolerances to continue. Investigate and record the cause. A completed validation writes `VALIDATED.txt`. Every dataset is checkpointed, so interrupted work can resume with identical source hashes. Run only one controller for a given output directory.

The scoring script copies the existing anchors-only and refineR results and uses RIbench's unchanged native scoring functions. The primary success criteria use estimates before the safeguard. Reports retain failure denominators and list endpoint extremes separately from the native mean-absolute-z definition. No N=500,000 data are included in this study.

## Candidate models and transitions between stages

Each local checkpoint includes the BIC, fitted limits and rank of every candidate and the gap to the runner-up. Baseline candidate tables, the share within two BIC points and all corresponding limit changes are exported when validation finishes. Both experimental arms have the same reports after full fitting. The two Hb examples are in `validation/hb-100-148-local-candidates.csv`.

The GitHub Actions replay runs on `study/boxcox-c2-runtime`. Its workflow is copied from `replay-workflow.yml` under this study directory. A preflight replays Hb 100 and 148, compares with the original estimates and records the image and numerical libraries. The original C2 candidate BICs were not saved. Once local validation passes, run `Rscript tools/boxcox/replay-release.R <local-output> <release-file>` and commit that generated file as `anchoredRI/tools/boxcox/replay-release.csv` **on the CI branch** to release the 156 Arm A repeats. Download the `armA-replay-156` artifact and run `Rscript tools/boxcox/compare-runtimes.R config.rds <extracted-artifact>` when the local Arm A fits for all 156 are ready. The comparison verifies both data and model-source hashes and saves every paired result, exception and largest difference. Only these 156 cases enter the cross-runtime Arm A comparison.
