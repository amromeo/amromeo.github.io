# Box-Cox shape and mixture capacity study

Implementation of the approved Task E protocol. The package code in `R/` is unchanged. The baseline remains `55eee86b3f4111197b0043cc062ef0b781dc2057` and the one-group adapter is copied from `f72ea14fc7978bf2478f54c87548cb6c127fa766`.

Run from the package directory with RIbench 1.0.2 and the original study's R dependencies. `frozen/` contains verbatim sources. `capacity-fixed.R` and `capacity-shifted.R` add diagnostic returns and extend the existing initialization patterns to four and five components. Restricting them to three components must reproduce the frozen estimator.

`models.R` implements the shared Box-Cox likelihood, analytic gradients, original-scale Jacobian and BIC, fixed-shape checks, and study-local safeguards. Every candidate retains its starts, optimizer results and bound hits. For numerical stability each shape optimization uses fixed centering and scaling from its starting transform. SD bounds are .005 to 5 on that numerical scale, as in the shifted fitter. Mean bounds cover the observed range at the endpoints of the shape and offset search regions. The parameter count excludes these fixed numerical coordinates.

Fixed lambda 1 uses the identity optimizer with the required translation of means. Fixed lambda 0 uses the log optimizer. An independently written Box-Cox likelihood verifies every candidate's original-scale log likelihood. The tests also check the analytic gradient against finite differences, multiplicities and inverse-transform domains.

The configuration RDS contains `baseline`, `output`, `workers`, `small_archive`, `large_archive`, `small_data`, `large_data`, and `simulation_archive`. The data directories contain `Data/<Analyte>/*.csv`. Archives are the completed C2 small-N and N=50,000 artifacts. The simulation archive is the original `repeated-fits.csv`, with all 45 cases and replicates 1–20.

```
Rscript tools/boxcox/test-models.R
Rscript tools/boxcox/run-study.R config.rds validate
Rscript tools/boxcox/run-study.R config.rds fit
Rscript tools/boxcox/score-study.R config.rds
```

Validation checks all 4,320 independent anchor hashes before fitting. Every routine vector must match RIbench's original hash. All 5,220 datasets must pass candidate and selected-fit equivalence before experimental fitting begins. The same-machine check has absolute tolerance 1e-10 and requires identical selected transform and component count. Cross-machine comparisons use the approved 0.1%-of-archived-width tolerance and report every exception. Any selected model change or more than 1% outside that tolerance pauses the study for review. Frozen C2 is refitted locally and supplies the comparison baseline. Archived anchors-only and refineR remain unchanged. See the dated protocol amendment.

A mismatch writes `HALT.txt`. Do not remove this file or loosen tolerances to continue. Investigate and record the cause. A completed validation writes `VALIDATED.txt`. Every dataset is checkpointed, so interrupted work can resume with identical source hashes. Run only one controller for a given output directory.

The scoring script copies the existing anchors-only and refineR results and uses RIbench's unchanged native scoring functions. The primary success criteria use estimates before the safeguard. Reports retain failure denominators and list endpoint extremes separately from the native mean-absolute-z definition. No N=500,000 data are included in this study.

The first numerical tests passed. A synthetic fitted-shape smoke check also passed. Full validation then stopped on an archived-baseline mismatch, documented in STATUS.md. That does not establish benchmark performance. Results will be available only after validation, fitting and scoring finish.
