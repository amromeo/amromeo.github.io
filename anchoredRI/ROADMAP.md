# anchoredRI roadmap

Last updated: November 25, 2026 · Current version: 0.4.0

This file records where the project stands, the decisions already made, and the studies still to do. Update it when a study finishes or a plan changes. Each study keeps its own protocol and report.

## Where things stand

- **What the package does.** `fit_anchored_ri()` estimates a reference interval from routine results plus results from 40 chart-reviewed patients who meet the reference criteria. It tries anchoredRI first. If that fit cannot report an interval, it tries a fallback that adjusts a refineR estimate toward the screened patients and checks it against their results. The fallback was developed in the [anchorRefine repository](https://github.com/amromeo/anchorRefine). If neither method reports an interval, the package suggests collecting at least 120 screened patients for a direct interval or verifying a suitable published interval with 20.
- **Evidence.** The combined method passed all 16 prespecified conditions on 7,488 independent simulated datasets at 1,000, 1,500, 5,000 and 50,000 routine results. Depending on dataset size, it returned an interval on 98.2–99.5% of confirmation datasets and 98.7–99.7% of RIbench datasets. The RIbench comparison is descriptive and used 4,320 datasets. [Task G report](https://github.com/amromeo/amromeo.github.io/blob/655a007f6d5ec67d41262c45f21293b3d0b945b5/anchoredRI/tools/combined/results/REPORT.md) · [Accuracy and return rates](https://github.com/amromeo/amromeo.github.io/blob/655a007f6d5ec67d41262c45f21293b3d0b945b5/anchoredRI/tools/combined/results/accuracy.csv)
- **What remains untested.** All evaluations so far use simulated data. The combined method has been evaluated with exactly 40 screened patients. Confidence-interval coverage is unknown. The combined method has yet to be tested directly between 5,000 and 50,000 routine results or above 50,000.

## Decisions already made

1. **No fallback to the screened patients alone.** When neither method can report an interval, the package suggests a direct study with at least 120 screened patients or verification of a suitable published interval with 20.
2. **One rule at every dataset size.** Try anchoredRI first and the fallback second. The package makes this decision without a known share of patients with disease.
3. **Keep anchoredRI first.** With 20% of patients with disease or less, trying the fallback first had slightly lower loss at each RIbench size, with all three paired 95% intervals including zero. In the same subgroup of the independent confirmation data, trying anchoredRI first had lower loss at each size, with all four paired 95% intervals below zero. Loss counts a missing interval as an error of 5. Study 10 tests whether choosing the order from the data does better. [Task G order comparison](https://github.com/amromeo/amromeo.github.io/blob/655a007f6d5ec67d41262c45f21293b3d0b945b5/anchoredRI/tools/combined/results/REPORT.md#a-versus-r-co-primary-order-comparison)
4. **Include settings where refineR struggles.** The studies include cases where most patients have disease. Report results separately for a share of patients with disease of 20% or less and more than 20%.
5. **Keep anchorRefine in its own research repository.** Study changes to the fallback there. Bring package changes here through a pull request.

## How studies are run

- Write the protocol before running anything. Set the question, datasets, success rule and stopping rule. Record any departures.
- Choose settings on development data and judge them on separate confirmation data. Keep the confirmation data out of tuning.
- Use RIbench as a descriptive check alongside the confirmation data unless the protocol specifies a different role.
- Score a dataset with no interval as an error of 5 when comparing overall performance. Report the share of datasets with an interval and the error among returned intervals separately.
- Repeat a prespecified audit sample in the same pinned container used for the main run. Apply the audit’s prespecified matching requirements.
- Work on a branch and open a pull request. Am reviews and merges. Keep direct pushes off `main`, which serves the live site.

## Planned work

| # | Study | Question | Status | Priority |
|---:|---|---|---|---|
| 1 | Website revision | Can a first-time visitor see what the method does, how well it works and how to install it? | Proposal for review | Now |
| 2 | Screened sample size | How does the method perform with 20, 30, 60, 80 and 120 screened patients? | Not started | Next |
| 3 | Confidence intervals | Do the 95% confidence intervals cover the true limits at the stated rate, and can fallback intervals get confidence intervals? | Not started | Next |
| 4 | Untested dataset sizes | How accurate is it between 5,000 and 50,000 routine results and above 50,000? | Not started | Next |
| 5 | Imperfect screening | What happens when some screened patients have disease or were chosen partly by their result? | Not started | Next |
| 6 | Blend weight by size | Should the fallback's adjustment toward the screened patients depend on dataset size? | Not started | Later |
| 7 | Measurement units | Can the fallback give the same interval across measurement units? | Not started | Later |
| 8 | Same answer on every computer | Can model selection be made stable across computers? | Not started | Later |
| 9 | Real pediatric data | Does it work on real CHOP data, judged against independently screened children? | Not started | Later (start planning now) |
| 10 | Order chosen from the data | Should large datasets that look mostly free of disease try the fallback first? | Not started | Later |

### 1. Website revision

Rebuild the home page around the current method, add an accuracy table from the published Task G files, use the same menu on every page, move older studies under “Earlier studies,” and remove internal study labels from the main text. Use the published fits and simulations.

**Done when** the pull request is merged and every number on the new pages links to its source file.

### 2. Screened sample size

The combined method was evaluated with exactly 40 screened patients. The package warns at other counts and skips the fallback below 40.

**Plan.** Test 20, 30, 60, 80 and 120 screened patients, with 40 as the reference, for both methods at the same routine dataset sizes as Task G. Compare against a direct interval from 120 screened patients.

**Done when** we know how often an interval is returned and how accurate it is at each count, and can set the warning and fallback minimum from that evidence.

### 3. Confidence intervals

The bootstrap gives confidence intervals for limits supplied by anchoredRI. Coverage has yet to be evaluated. Fallback limits currently come without confidence intervals.

**Plan.** Measure whether the confidence interval for each limit covers the true limit at the stated rate across dataset sizes and shares of patients with disease. Design and test a bootstrap for the fallback.

**Done when** coverage is reported by size and share of patients with disease, and the package gives fallback confidence intervals or explains why they are unavailable.

### 4. Untested dataset sizes

Task G tested 1,000, 1,500, 5,000 and 50,000 routine results. The package prints a message above 5,000 about the untested range between 5,000 and 50,000.

**Plan.** Test sizes such as 10,000 and 20,000, and at least one size above 50,000, with the frozen 0.4.0 method.

**Done when** the size message can be removed or narrowed based on the results.

### 5. Imperfect screening

The combined-method simulations assumed that every screened patient met the reference criteria. In practice, chart review can miss relevant disease, and seeing a result can influence whether a reviewer includes the patient.

**Plan.** Mix patients with disease into the screened sample at several rates, such as 2.5%, 5% and 10%. Include a scenario where screened patients are chosen partly by their result. Measure the effect on both methods and on their decisions to report an interval.

**Done when** the help page can describe how screening errors affect performance.

### 6. Blend weight by dataset size

The fallback adjusts the refineR estimate toward the screened patients with a fixed weight of 0.25. In development, a weight of 0.5 reduced control error at smaller sizes and increased it by about 0.043 at 50,000 routine results. [Fallback study report](https://github.com/amromeo/anchorRefine/blob/638e2d3500faa01d4e2816895a95d7a6685567a0/results/stage3-confirmation-45463873162/REPORT.md)

**Plan.** Use new development data to choose a size-dependent rule and new confirmation data to evaluate it. The anchorRefine repository is the likely home for this study.

**Done when** the new rule beats the fixed weight on independent confirmation data or we record why the package keeps the fixed weight.

### 7. Measurement units

The fallback inherits a small sensitivity to measurement units from refineR. Changing units and converting back moved a limit by more than 0.1% of the interval width on 36 of 720 development datasets for the blend the package uses. This measurement covers the fallback. [Fallback study report](https://github.com/amromeo/anchorRefine/blob/638e2d3500faa01d4e2816895a95d7a6685567a0/results/stage3-confirmation-45463873162/REPORT.md)

**Plan.** Find the cause in the refineR fitting path and test a unit-invariant version. Measure anchoredRI's unit sensitivity the same way.

**Done when** both methods give the same interval within a stated tolerance across units, or the remaining sensitivity is documented.

### 8. Same answer on every computer

Model choice can occasionally differ between computers. RIbench calcium dataset 870 is an example from earlier validation. In Task G, the first confirmation audit failed when jobs ran on different system images. The repeat matched in the pinned container. [Calcium example](README.md) · [Task G audit](https://github.com/amromeo/amromeo.github.io/blob/655a007f6d5ec67d41262c45f21293b3d0b945b5/anchoredRI/tools/combined/results/REPORT.md#study-population-and-reproducibility)

**Plan.** Measure how often the reported interval changes across Windows, macOS and Linux with recent R versions. Test whether tie-breaking rules or more optimizer starts make model choice more stable.

**Done when** we can state how often results differ across computers and by how much, and have reduced differences where possible.

### 9. Real pediatric data

This study tests clinical performance and will take time to set up.

**Plan.** Choose a few analytes and age groups. Agree on eligibility rules, comparisons and a success rule in a protocol, and obtain approvals for patient data use before collecting data. Collect routine results from CHOP and chart-review two independent groups of children. Use 40 screened children as inputs to the method and use a separate group of at least 120 to estimate the interval directly for comparison.

**Done when** the method's intervals are compared with the direct intervals on real data and the findings are reported regardless of the outcome.

### 10. Order chosen from the data

The package always tries anchoredRI first. On RIbench, with 20% or fewer patients with disease, refineR's limits were closer to the true values at 5,000 and 50,000 routine results. At 50,000 the error was 0.095 for refineR alone against 0.117 for the combined method. With more than 20% disease, refineR did much worse (0.365 against 0.173). The package can't know the true share of patients with disease, but the first fit produces an estimate of it. [RIbench evidence](https://amromeo.github.io/anchoredRI/evidence.html#ribench) · [Numbers and sources](https://amromeo.github.io/anchoredRI/website-sources.html#accuracy-table-and-figure)

**Plan.** Use the first fit's estimate of the share of routine results outside the reference group, together with dataset size, to decide which method to try first. Measure how well that estimate tracks the true share, especially when disease overlaps heavily with the reference group. Choose the rule and its cut-offs on new development data. Then test it against the fixed anchoredRI-first order on new confirmation data not used in Task G.

**Done when** the rule beats the fixed order on independent confirmation data without losing accuracy when disease is common, or we keep the fixed order and record why.

## Known limits with no study planned yet

- **Heavy overlap when most patients have disease.** A wrong interval can pass both anchoredRI checks when the results of patients with disease overlap heavily with the reference distribution. An ICU-only cohort could have this pattern.
- **Reference distributions outside the model families.** The confirmation data included a Box-Cox shape (λ = 0.5) that anchoredRI's models cannot represent exactly. Other shapes remain to be studied systematically.
- **Age as a continuous variable.** The method estimates one interval for one age group. Pediatric intervals often change smoothly with age. Extending the method to continuous age would be a separate project.
- **Distribution.** The package is installed from GitHub or the downloadable tarball. Start with a tagged release for each version. CRAN or r-universe could follow real-data validation.

## Related records

- [Task G protocol](https://github.com/amromeo/amromeo.github.io/blob/655a007f6d5ec67d41262c45f21293b3d0b945b5/anchoredRI/tools/combined/PROTOCOL.md) and [report](https://github.com/amromeo/amromeo.github.io/blob/655a007f6d5ec67d41262c45f21293b3d0b945b5/anchoredRI/tools/combined/results/REPORT.md)
- [anchorRefine repository](https://github.com/amromeo/anchorRefine) for fallback development
- [Change log](NEWS.md)
