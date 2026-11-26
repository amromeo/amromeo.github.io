# Task H subsampling results

The cohort contains 1,440 RIbench datasets at nominal N = 50,000. Each dataset contributes one equally weighted comparison. Each subsample draws 5,000 routine observations and retains its 40 original anchors. Arm W fits the full routine vector once with a total routine weight of 5,000 and an anchor weight of 1 per person. Both W and the subsamples use 5,040 as the BIC sample size. The standard fit allows up to five components.

| Standard full first stage | Cap | Eligible | Majority releases | W releases | Common triples | Mean error change, median | Mean error change, W | Extreme intervals, reference / median / W |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| Release | 5 | 1126 | 1125 | 1118 | 1126 | -0.021 | -0.004 | 0 / 0 / 0 |
| Decline | 5 | 314 | 263 | 252 | 250 | 0.163 | 0.209 | 0 / 0 / 0 |
| Release | 3 | 1126 | 1019 | 989 | 1126 | 0.286 | 0.402 | 0 / 14 / 17 |
| Decline | 3 | 314 | 154 | 151 | 147 | -0.105 | -0.091 | 0 / 0 / 0 |

Seed set A uses the coordinatewise medians of all usable raw draws. Panel A compares with the five-component-cap full first-stage estimate. In panel B the median must pass its majority check and the full package must return a fallback. 0 full first stages had technical withholdings. The endpoint columns count intervals with any requested |z| > 5. W contributes its raw limits to the common triples. Its own two checks determine its release count. [All split panels and W’s standalone releases](paired-summary.csv) · [Each pair](paired-index.csv) · [Limit changes by analyte](paired-limits-by-analyte.csv)

For the 1,126 datasets whose standard full fit passed, mean endpoint error was 0.131 for the full fit, 0.110 for the subsample median and 0.127 for W. The order reversed when the full fit declined and both the majority and the package fallback released an interval. On those 250 datasets, mean error was 0.233 for the package fallback, 0.396 for the median and 0.442 for W. None of these cap-5 estimates had an endpoint with |z| above 5 in either comparison. [Matched comparisons](paired-summary.csv)

## Accuracy and release

| Estimator and gate | Returned / assigned | Declined | Technical | Pre-check mean error | Released mean error | Released cutoff mean (n) | All-assignment score after checking | Released extreme intervals |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Full, first stage | 1,126 / 1,440 | 314 | 0 | 0.411 | 0.131 | 0.131 (1,126) | 1.193 | 0 |
| Full, package | 1,423 / 1,440 | 17 | 0 | 0.148 | 0.148 | 0.148 (1,423) | 0.205 | 0 |
| W, cap 5, own checks | 1,370 / 1,440 | 70 | 0 | 0.198 | 0.179 | 0.179 (1,370) | 0.413 | 0 |
| Median, cap 5, full gate | 1,126 / 1,440 | 314 | 0 | 0.174 | 0.110 | 0.110 (1,126) | 1.176 | 0 |
| Median, cap 5, majority | 1,388 / 1,440 | 52 | 0 | 0.174 | 0.163 | 0.163 (1,388) | 0.337 | 0 |
| Majority median or package, cap 5 | 1,436 / 1,440 | 4 | 0 | 0.160 | 0.160 | 0.160 (1,436) | 0.174 | 0 |
| Full, cap 3, first stage | 1,069 / 1,440 | 371 | 0 | 0.851 | 0.172 | 0.172 (1,069) | 1.416 | 0 |
| W, cap 3, own checks | 1,140 / 1,440 | 300 | 0 | 0.513 | 0.163 | 0.163 (1,140) | 1.171 | 0 |
| Median, cap 3, full gate | 1,069 / 1,440 | 371 | 0 | 0.421 | 0.168 | 0.148 (1,068) | 1.401 | 1 |
| Median, cap 3, majority | 1,173 / 1,440 | 267 | 0 | 0.421 | 0.153 | 0.153 (1,173) | 1.051 | 0 |
| Majority median or package, cap 3 | 1,430 / 1,440 | 10 | 0 | 0.153 | 0.153 | 0.153 (1,430) | 0.187 | 0 |
| Archived Task E Arm B | 1,135 / 1,440 | 305 | 0 | 0.403 | 0.131 | 0.131 (1,135) | 1.162 | 0 |
| Archived Task E baseline | 1,066 / 1,440 | 374 | 0 | 0.925 | 0.171 | 0.171 (1,066) | 1.425 | 0 |
| refineR, native limits | 1,440 / 1,440 | 0 | 0 | 0.230 | 0.230 | 0.230 (1,440) | 0.230 | 1 |
| anchorRefine, w = 1 | 1,345 / 1,440 | 95 | 0 | 0.196 | 0.171 | 0.171 (1,345) | 0.490 | 0 |
| anchorRefine, t = 0.25 | 1,387 / 1,440 | 53 | 0 | 0.198 | 0.176 | 0.176 (1,387) | 0.354 | 0 |
| Anchors only | 1,440 / 1,440 | 0 | 0 | 0.296 | 0.296 | 0.296 (1,440) | 0.296 | 0 |

The dataset error is the mean absolute deviation on RIbench’s z scale over requested endpoints, with the upper limit alone for CRP. The native cutoff excludes errors above 5. The all-assignment score caps errors at 5 and assigns 5 to a nonreturn. For the package and combined policy, the score uses the actual returned endpoints over the entire cohort. Native refineR limits are primary. Archived anchorRefine rows use the primary Stage 1 check. The historical table reproduces Task E’s original extraction. [All strata and filters](task-e-metrics.csv.gz) · [Historical Task E results](task-e-historical-reproduction.csv) · [Historical paired differences](historical-paired-contrasts.csv.gz) · [Paired differences, wins, ties, losses and coverage](paired-contrasts.csv.gz)

The cap-5 majority-median-or-package policy kept every one of the full package’s 1,423 returns and added 13. On those same 1,423 datasets, it raised mean endpoint error by 0.012. It improved 688 intervals, tied on 48 and worsened 687. With each nonreturn scored as 5, the score across all 1,440 fell from 0.205 to 0.174. Neither policy released an interval with an extreme endpoint. [Paired comparisons](paired-contrasts.csv.gz)

## Safeguards, component caps and seed stability

- At cap 5, among the 314 standard full-data statistical declines, the majority released 13 when the package withheld and the package released 47 when the majority withheld. The combined policy returned 1,388 medians, 1 full first-stage estimate and 47 released fallbacks, with 4 nonreturns.
- At cap 3, among the 314 standard full-data statistical declines, the majority released 7 when the package withheld and the package released 150 when the majority withheld. The combined policy returned 1,173 medians, 107 full first-stage estimates and 150 released fallbacks, with 10 nonreturns. 18 primary majority decisions changed between seed sets.

Changing the subsample seeds changed the primary cap-5 majority decision on 9 of 1,440 datasets (0.63%). The 95th percentile of the absolute movement between the two raw medians was 0.067 z units for the lower limit and 0.084 for the upper. [Seed movement summaries](seed-movement-summary.csv)

Each model selection uses the same ordered routine values, anchors and positioning seeds across the two caps. The movement tables give B-minus-A changes by side in analyte units and absolute changes on the z scale. They include passing-fits-only medians, the largest movements by Index, and the subset released under both seeds. [Seed movement summaries](seed-movement-summary.csv) · [Per-Index movement and votes](seed-movement.csv.gz) · [Seed decision changes](seed-decisions.csv) · [Cap and gate decision changes](release-changes.csv) · [Changed component counts](cap-selection-changes.csv) · [Component-count pairs](cap-component-counts.csv)

The vote includes all 50 draws. Technical votes can leave the majority unresolved. A statistically declined fit can contribute valid raw limits to the primary median. A missing passing-only median takes the saved full-package route. [Vote and flag records](gate-records.csv.gz) · [Missing contributors](median-contributors-below-50.csv.gz) · [Policy routes](policy-routes.csv) · [Full versus majority](full-vs-majority.csv) · [Flags by error](flags-by-error.csv)

## Endpoint and fit records

Each endpoint and fit uses the frozen RIbench 1.0.2 scorer. For individual draws, each cap, seed set and draw number defines one raw method on the cohort. All scores precede release masking. The extreme lists retain declined raw intervals and identify each gate. [Endpoint scores](endpoint-scores.csv.gz) · [Fit scores](fit-scores.csv.gz) · [Method decisions](method-decisions.csv.gz) · [Every extreme endpoint](extreme-endpoints.csv.gz) · [Datasets with mean error above 5](datasets-above-five.csv.gz)

Before safeguards, 11 cap-5 full fits had an extreme endpoint, compared with zero cap-5 W fits or primary medians in either seed set. Individual subsamples could still pass the safeguards with an extreme endpoint. Five such cap-5 draws came from LDH dataset 5609, with lower-limit z from +5.64 to +7.31. The cap-3 primary median had 17 extreme intervals under seed set A and 15 under set B before the checks. [Extreme endpoints and decisions](extreme-endpoints.csv.gz)

The physiological-share analysis pairs each fit with the unweighted fit at the same cap. The optimizer analysis compares iteration-limit and line-search exits with converged fits. [Component shares](shares/REPORT.md) · [Paired share and error changes](shares/paired-fits.csv.gz) · [Optimizer exits](optimizer/OPTIMIZER.md)

Across both caps, 19,332 selected fits reached the iteration limit and 798 had a line-search exit. Among cap-5 full fits, 479 of the 481 iteration-limit exits selected four or five components. In cap-5 subsamples, the mean error for iteration-limit fits exceeded the converged-fit mean by 0.109 under seed set A and 0.099 under set B. Comparing draws from the same dataset reduced these differences to +0.006 and −0.002. [Exits and endpoint errors by fit type, cap and component count](optimizer/OPTIMIZER.md)

## Verification, compute and data

The frozen implementation was `558001e`. Launch commit `6442ca4` assigned the remaining 1,420 datasets to 178 consecutive Index blocks. All 187 Actions jobs succeeded on the first attempt. They used 106.96 runner-hours, including setup and bundling, and the dataset processes used 104.27 CPU-hours. From the first job starting to the last job ending, 5 hours 46 minutes 42 seconds elapsed. [Launch receipt](../../launch-full.json) · [Actions run](https://github.com/amromeo/amromeo.github.io/actions/runs/45860591866) · [Compute totals](compute/actual-compute.json) · [Phase times](compute/phase-times.csv)

The readback covered all 1,440 datasets, including the 20 original pilots. Verification checked 12,960 source output files, 144,000 draw seeds, 720 million stored positions, 293,760 cap-specific fit records and 69,120 decisions. In the earlier cap-3 audit, 1,354 of 1,440 matched the archived Task E baseline on every check. Rerunning the 86 exceptions with the archive seed resolved 37 and left 49 exceptions. Both comparisons remain available by Index. [Verification counts](verification/summary.json) · [Cap-3 audit and exceptions](../pilot-45854822638/archive/AUDIT.md)

The nine release bundles hold the complete per-dataset records. The pilot bundle preserves the original 20. The eight numbered bundles hold the remaining 1,420. The release includes a SHA-256 sidecar for each bundle. [Permanent data release](https://github.com/amromeo/amromeo.github.io/releases/tag/task-h-subsample-50000) · [Bundle sizes and dataset counts](release-assets.csv) · [Release SHA-256 checksums](RELEASE-SHA256SUMS) · [Summary file SHA-256 checksums](SHA256SUMS)

RIbench informed the subsampling proposal. Adoption would require independent confirmation and a prespecified sample-size rule for 5,000–50,000 routine results.
