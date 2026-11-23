# Task F safeguard validation

The package now applies the Step B likelihood-ratio test and observed-minimum check by default. The implementation is on `study/safeguard`, based on main `55eee86b3f4111197b0043cc062ef0b781dc2057`. The optimization algorithm is unchanged. One-group support uses the already tested C2 adaptation.

## Archived decisions

All 5,220 archived decisions reproduced exactly using the package's check and reporting functions. These checks replay the saved fits rather than refitting them on a different numerical runtime.

| Dataset | Fits | LRT flags | Lower-limit flags | Either flag | Both flags |
|---|---:|---:|---:|---:|---:|
| B2 confirmation, replicates 21–40 | 900 | 91 | 0 | 91 | 0 |
| Independent-anchor RIbench, N <= 5,000 | 2,880 | 533 | 9 | 533 | 9 |
| Independent-anchor RIbench, N = 50,000 | 1,440 | 402 | 21 | 402 | 21 |

B2's committed rule a at threshold 0.80 declined 122/900. Its LRT-alone comparator declined 91/900, including 86/312 datasets with error >10% and 5/453 with error <5%. Task F implements the LRT plus observed-minimum rule, not B2's separation statistic. The lower-limit check adds no B2 declines in this confirmation set.

For small-N RIbench, the combined rule flags all 23 fits whose mean absolute z deviation exceeds 5 and 224/1,952 whose mean absolute z deviation is below 0.25. For the alternative endpoint definition, it flags 32/33 fits with either endpoint |z| >5. Those definitions are deliberately kept separate.

CRP is one-sided. Its archive contains only the reported upper endpoint, so the full reference-component parameters cannot be recovered. For CRP the package calculates the p-value from the archived LRT statistic and applies the one-sided reporting rule. For the remaining RIbench fits and every B2 fit, the package recomputes the LRT from the saved central limits and regenerated original anchors. All flags agree with the archives. Public-function tests independently check LRT computation and the one-sided skip on complete fits.

## Package behavior

Tests cover each decline reason, simultaneous flags, both diagnostic estimates, no anchors-only substitution, technical withholding, the strict lower-limit inequality, one-sided requests, the research override, and the refineR recommendation beginning at 50,000 routine results. Bootstrap tests cover every draw's checks and diagnostics, accepted-draw confidence intervals, count accounting, random-state restoration and agreement between sequential and parallel execution.

`successful` counts accepted bootstrap draws. `failed` counts numerical failures. `declined_draws` and `withheld_draws` are separate. An accepted point fit needs at least 80% accepted bootstrap draws for confidence intervals. A declined point fit reports neither the reference interval nor its confidence intervals. These intervals are conditional on passing the safeguard.

The safeguard can miss wrong intervals in heavily mixed populations with many mildly abnormal patients, including ICU-only populations. A passing result does not establish clinical validity.

## Reproduction

R 4.5.0 completed `R CMD check --no-manual` with all package tests passing. The three two-setting parity tests (identity, log and shifted-log) had zero difference in their limits. Run the package tests after installation. `test-parity.R` compares research-mode fits against frozen main source on the same runtime. `replay-archives.R` accepts an RDS configuration with paths named `b2`, `small`, `large`, `small_diagnostics`, `large_diagnostics`, and `output`. It expects the original C2 small-N artifact, the completed 50,000 artifact, the B2 confirmation output, and Task D's lower-limit diagnostics. RIbench 1.0.2 is required only for this archive replay, not for package use.

The branch includes the aggregate decline counts. Running the replay also writes a local CSV with one decision per dataset.
