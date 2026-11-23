# Task E fitting reliability

Neither improved fitter met all three requirements on the 156 hemoglobin datasets. The baseline met the model and endpoint agreement targets. Its shifted-log model failed the nesting check. Arm A passed nesting and missed both agreement targets. The 5,220-dataset accuracy comparison is on hold.

## Improved fitters on the same 156 datasets

Both runtimes fit the same 156 routine vectors and the same 40 screened anchors per dataset. A model match required the same transform or shift specification and the same component count. Endpoint agreement required both limits to be within 0.5% of the Actions interval width. All 156 interval pairs were available for each method.

| Prespecified requirement | Improved baseline | Arm A, fitted Box-Cox |
|---|---:|---:|
| Same model in at least 152/156 | 155/156 (99.4%) | 132/156 (84.6%) |
| Both limits agree in at least 149/156 | 155/156 (99.4%) | 91/156 (58.3%) |
| Zero nesting violations locally | 3 in 1,206 comparisons | 0 in 3,036 |
| Zero nesting violations on Actions | 4 in 1,206 comparisons | 0 in 3,036 |

Every eligible likelihood comparison had both estimates. For the baseline, Hb 102 was the one dataset where the model and endpoints disagreed. The local fit selected three identity-scale components and Actions selected three log-scale components. The largest endpoint difference was 0.187 Hb units, or 4.54% of the Actions interval width.

Three baseline nesting violations recurred on both machines, at Hb 22, 82 and 118. Actions also flagged Hb 102. In every case, the shift range included zero. The unshifted log fit reached a component SD between 0.00362 and 0.00439. The shifted-log fitter imposed a minimum SD of 0.005 on that scale and therefore could not represent the unshifted solution. All seven observed inversions count toward the prespecified gate. These examples show a mismatch in the numerical bounds of the two models.

For Arm A, 42 datasets exceeded the endpoint tolerance even when both runtimes selected the same model. The largest absolute endpoint difference was 1.006 Hb units at Hb 85, where both selected a three-component shifted model. The largest relative difference was 22.24% of the Actions interval width at Hb 21. At the same 0.5% tolerance, the original Arm A replay had 96/156 agreeing limit pairs and 131/156 matching models. Multistart v2.1 gave 91 and 132, respectively.

The fitter retained the highest finite feasible likelihood attained. In 74/156 local Arm A fits, the winning optimization had a nonzero convergence code. Seventy-three hit the 1,000-iteration limit and one ended with a line-search error. On Actions, 69/156 winners hit the iteration limit. The corresponding baseline counts were one locally and two on Actions, all line-search errors. A retained likelihood is an attained value. The optimizer does not certify a global maximum.

[Individual paired limits, model choices and BIC gaps](paired-replay.csv) · [All candidate BICs and likelihoods](candidates) · [Every nesting violation and its bound](violations.csv) · [Convergence by method and runtime](convergence.csv) · [Convergence codes for every start](convergence-codes)

## What the saved fits showed before refitting

The original local baseline had a nesting violation in 154/2,880 small-N RIbench datasets (5.35%), 73/1,440 N = 50,000 datasets (5.07%) and 72/900 simulations (8.00%). Across the study, that was 299/5,220 datasets (5.73%). In the 156-dataset Arm A replay, 31/156 were flagged in each runtime. Fourteen were flagged on both.

Giving each candidate the best likelihood already demonstrated by a nested model changed none of the baseline winners or the 156 replay winners. Among the 384 additional saved local Arm A fits, it changed 15 BIC choices. This calculation cannot establish which model a correct global fit would select. The v2.1 refit actually changed 24/156 baseline choices in each runtime relative to its historical reference. For Arm A, 50/156 choices changed locally and 45/156 on Actions. The Actions baseline reference is the original archived C2 choice. The endpoint reproducibility tests above are comparisons of two independent v2.1 fits.

[Saved-fit audit by dataset group and analyte](../saved-audit/REPORT.md) · [Original and improved choice for each dataset](original-and-improved-selections.csv)

## Reproducibility records

The corrected fitter used the original fixed-model and grid starting recipes, nested solutions, a split of the widest component and 20 Index-seeded perturbations per candidate. The original pilot omitted frozen seed fits and applied the original SD start limits incorrectly. Its results and the reason for repeating it are recorded in the protocol. The corrected run used fresh outputs and the same thresholds.

The Actions job matched the archived C2 image, R version, numerical libraries and required package versions. All eight fitting-source hashes and all 156 pairs of routine-data and anchor hashes matched across runtimes. The downloaded artifact contained 156 fit records and matched the Actions SHA-256 digest `2a0f2dcf0f455c9096300a884f98b968ad374047d4a7b693c2432e01fb254acd`.

We also compared the starting vectors with their historical counterparts on each runtime. All 26,208 grid starts matched to 1e-10 in both environments. All 1,907 local frozen-model seeds matched. On Actions, 1,546 of 1,905 paired frozen-model seeds matched their older values, and one historical shifted-log seed for Hb 148 was absent from the new inventory. The same source for generating a seed can itself reach a different fitted solution. This comparison uses the 156 Hb datasets and does not establish reproducibility for other analytes.

The convergence files give the number of starts, the chosen start and every nonzero start index and code. All unlisted start indices in each row have code zero. The Actions artifact also retains every starting vector, fitted vector, optimizer message and numerical bound.

[Corrected Actions run and raw fits](https://github.com/amromeo/amromeo.github.io/actions/runs/45175838337) · [Fitter source commit](https://github.com/amromeo/amromeo.github.io/commit/fc420fe0dea22f4aad8a2c8e921a810f88a9d828) · [Original-start checks](original-starts-summary.csv) · [Original-start exceptions](original-start-exceptions.csv) · [Source hashes](source-hashes.csv) · [Protocol](../../PROTOCOL.md)

## Next reliability check for review

For the baseline, align the log and shifted-log SD bounds at zero shift, then verify the identical parameter vector gives the identical likelihood. For Arm A, exchange the best parameter vectors between runtimes on the discordant Hb datasets and fit from those common starts. Compare the likelihood, limits and convergence status. This will help distinguish a missed better solution from several parameter sets giving nearly the same likelihood and different limits. The accuracy benchmark requires a separate release after the reliability criteria pass.
