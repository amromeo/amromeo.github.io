# Task G confirmation results

**16 of 16 conditions passed.** The combined method met the prespecified success rule.

The blend-first order did not have lower confirmation loss in any prespecified stratum. The study stops here for review.

Order flag on descriptive RIbench: R had lower loss at N = 1,000 (<=20%), N = 5,000 (<=20%), N = 50,000 (<=20%).

A tries safeguarded anchoredRI first. B is the checked Stage 3 post-fit blend alone. R tries B first and uses anchoredRI as fallback. A and R can return on the same datasets. Their order matters when both component methods return different intervals. C is checked refineR. All losses count a nonreturn as 5. Accuracy uses only paired common returns.

## A versus B: success conditions

Differences are A minus the comparator. Negative values favor A. Parentheses give the paired denominator. Pass/fail uses unrounded confirmation values.

| N | Pathology | Metric | Confirmation | RIbench (descriptive) |
|---:|---|---|---|---|
| 1,000 | <=20% | error | -0.0271 (563) pass | +0.0076 (685) |
| 1,000 | <=20% | loss | -0.1199 (576) pass | -0.2212 (720) |
| 1,000 | overall | error | -0.0068 (1,645) pass | -0.0970 (1,323) |
| 1,000 | overall | loss | -0.5114 (1,872) pass | -0.4580 (1,440) |
| 1,500 | <=20% | error | -0.0285 (564) pass | Unavailable |
| 1,500 | <=20% | loss | -0.1082 (576) pass | Unavailable |
| 1,500 | overall | error | -0.0080 (1,697) pass | Unavailable |
| 1,500 | overall | loss | -0.3680 (1,872) pass | Unavailable |
| 5,000 | <=20% | error | -0.0158 (573) pass | +0.0471 (709) |
| 5,000 | <=20% | loss | -0.0328 (576) pass | -0.0271 (720) |
| 5,000 | overall | error | -0.0136 (1,767) pass | -0.0422 (1,378) |
| 5,000 | overall | loss | -0.2162 (1,872) pass | -0.2127 (1,440) |
| 50,000 | <=20% | error | -0.0132 (575) pass | +0.0169 (702) |
| 50,000 | <=20% | loss | -0.0217 (576) pass | -0.0627 (720) |
| 50,000 | overall | error | -0.0400 (1,842) pass | -0.0337 (1,387) |
| 50,000 | overall | loss | -0.0904 (1,872) pass | -0.1450 (1,440) |

## A versus R: co-primary order comparison

Differences are A minus the comparator. Negative values favor A. Parentheses give the paired denominator. Brackets give paired 95% bootstrap intervals.

| N | Pathology | Metric | Confirmation | RIbench (descriptive) |
|---:|---|---|---|---|
| 1,000 | <=20% | error | -0.0266 (574) [-0.0394, -0.0132] | +0.0073 (720) [-0.0518, +0.0768] |
| 1,000 | <=20% | loss | -0.0265 (576) [-0.0393, -0.0132] | +0.0073 (720) [-0.0518, +0.0768] |
| 1,000 | overall | error | -0.0061 (1,843) [-0.0133, +0.0017] | -0.0894 (1,435) [-0.1504, -0.0254] |
| 1,000 | overall | loss | -0.0060 (1,872) [-0.0132, +0.0017] | -0.0891 (1,440) [-0.1499, -0.0253] |
| 1,500 | <=20% | error | -0.0280 (574) [-0.0399, -0.0161] | Unavailable |
| 1,500 | <=20% | loss | -0.0279 (576) [-0.0398, -0.0160] | Unavailable |
| 1,500 | overall | error | -0.0074 (1,839) [-0.0144, +0.0001] | Unavailable |
| 1,500 | overall | loss | -0.0072 (1,872) [-0.0141, +0.0001] | Unavailable |
| 5,000 | <=20% | error | -0.0157 (575) [-0.0241, -0.0073] | +0.0464 (720) [-0.0412, +0.1570] |
| 5,000 | <=20% | loss | -0.0157 (576) [-0.0241, -0.0073] | +0.0464 (720) [-0.0412, +0.1570] |
| 5,000 | overall | error | -0.0130 (1,846) [-0.0191, -0.0068] | -0.0407 (1,430) [-0.1111, +0.0423] |
| 5,000 | overall | loss | -0.0128 (1,872) [-0.0188, -0.0067] | -0.0404 (1,440) [-0.1100, +0.0421] |
| 50,000 | <=20% | error | -0.0131 (576) [-0.0195, -0.0068] | +0.0167 (714) [-0.0114, +0.0488] |
| 50,000 | <=20% | loss | -0.0131 (576) [-0.0195, -0.0068] | +0.0165 (720) [-0.0113, +0.0483] |
| 50,000 | overall | error | -0.0395 (1,862) [-0.0439, -0.0353] | -0.0329 (1,421) [-0.0657, -0.0029] |
| 50,000 | overall | loss | -0.0393 (1,872) [-0.0437, -0.0350] | -0.0325 (1,440) [-0.0648, -0.0029] |

## Study population and reproducibility

Confirmation comprised 7,488 Stage 3 synthetic datasets. The healthy families follow anchorRefine’s Box-Cox design, including lambda = 0.5, which anchoredRI’s transforms do not represent exactly. Extra replication of high-pathology shoulder and overlapping patterns contributes to the overall means. RIbench comprised 4,320 datasets and informed the original hypothesis. It has no N = 1,500 stratum. Its results are descriptive.

The first run used different Ubuntu image revisions across some paired jobs and failed its audit. All primary and audit jobs were repeated using the same immutable R 4.5.0 container and the same compiled library. No fits were selected between runs. All 749 prespecified audit repeats in the pinned run matched the primary anchoredRI fits and decisions at the protocol tolerance. Independent reconstruction checked all six routes against the shard releases. Standalone B and C were checked against the exported native scores. The paired intervals resample confirmation datasets within the original design cells and RIbench parent scenarios as clusters. [Actions run](https://github.com/amromeo/amromeo.github.io/actions/runs/45614043039) · [First-run discrepancies](first-audit-discrepancies.csv) · [Audit IDs](../audit.csv) · [Protocol](../PROTOCOL.md)

The tables give each method’s loss, returned-only error, return count and share by size, pathological fraction and disease pattern. Route counts identify which method supplied each interval. [All six methods by stratum](accuracy.csv) · [Route counts](routes.csv) · [Paired headline comparisons and denominators](headline.csv)

[Download the per-dataset scores, fit records, audit repeats and source manifests](task-g-evidence.zip).
