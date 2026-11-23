# Task E reliability round 3

The baseline passed the Hb reproducibility requirements. Stop Arm A under the prespecified weak-identification rule. The accuracy comparison remains on hold.

These checks ask whether the same routine and screened results give the same fitted reference limits on two computers. They cover 156 simulated hemoglobin datasets, each with 40 screened anchors. Arm B uses the 540 datasets with previously saved fits.

## 1. Baseline with a rounding-based SD floor

| Across the 156 Hb datasets | Result | Required |
| --- | --- | --- |
| Same transform and component count | 155/156 (99.4%) | At least 152/156 |
| Both reference limits within 0.5% of the Actions interval width | 155/156 (99.4%) | At least 149/156 |
| Eligible nesting violations | 0 local, 0 Actions | Zero in each runtime |

The Hb rounding step was 0.1. In each candidate model, a component had to span at least one rounding unit in original-scale SD. The shared bound also preserved the stronger pre-existing floor where it applied. At zero shift, identity, log and shifted log use the same original-unit floor at a given median. The likelihood identity is between log and shifted log evaluated at the same parameters with shift zero.

| Runtime | Fitted candidates | Eligible nesting checks | Excluded by shift bounds | Unavailable | Zero-shift vectors | Maximum likelihood difference | Maximum SD-floor difference |
| --- | --- | --- | --- | --- | --- | --- | --- |
| local | 1404 | 1206 | 198 | 0 | 468 | 0 | 0 |
| Actions | 1404 | 1206 | 198 | 0 | 468 | 0 | 0 |

The largest baseline endpoint difference was 0.881756 Hb units (14.02% of interval width) at Index 142. The local fit selected log with 3 components (11.92 to 17.39). Actions selected identity with 3 components (11.04 to 17.33). Of the 468 zero-shift vectors per runtime, 270 fell inside the permitted shift range. The others tested the algebraic identity only. All fitted component SDs met the shared floor to numerical precision. [Paired baseline fits](baseline-paired.csv) · [Nesting and floor summary](baseline-audit-summary.csv)

## 2. Arm A after continuing unfinished winners

| Runtime | Original iteration-limit winners | Original line-search-error winners | L-BFGS-B blocks | Bounded bridge calls | Datasets using affine retry | Final winners with code 0 |
| --- | --- | --- | --- | --- | --- | --- |
| local | 73 | 1 | 78 | 1 | 1 | 156/156 |
| Actions | 69 | 0 | 78 | 2 | 0 | 156/156 |

Of the original 42 datasets where both runtimes selected the same model but produced different limits, **29/42** still had limits outside tolerance after continuation. 13/42 came within tolerance. Counting either model or limit disagreement, 29/42 still differed.

Across all 156, 132/156 selected the same model and 97/156 had both limits within tolerance. 60 had a model or endpoint disagreement and entered the reciprocal-start experiment. Each unfinished winning fit was resumed from its last parameters with up to 10,000 iterations per block. A bounded bridge or curvature-scaled retry was used when a numerical line search stalled. Only an actual terminal code 0 counted as converged. [All 156 before-and-after records](arm-a-continued-pairs.csv) · [Local blocks](local-continuation-blocks.csv) · [Actions blocks](Actions-continuation-blocks.csv)

## 3. Arm A reciprocal starts

| Classification of reciprocal solutions | Cross-tested datasets | All 156 Hb datasets |
| --- | --- | --- |
| missed better fit | 47/60 (78.3%) | 47/156 (30.1%) |
| weak identification | 12/60 (20.0%) | 12/156 (7.7%) |
| concordant | 1/60 (1.7%) | 1/156 (0.6%) |

The stopping boundary was more than 5% of all 156 datasets, or eight datasets. At Hb Index 144, the log-likelihoods differed by just 0.105. The two lower limits were 11.35 locally and 10.95 on Actions. The largest endpoint change was 6.06% of the interval width.

The two runtimes completed 336 own-start and transferred-start solutions, all with code 0. The largest absolute difference between a transferred starting log-likelihood and its recomputation was 0. The largest reciprocal terminal log-likelihood difference was 227.996 at Index 150.

On the common selected model, the foreign start improved the best pre-exchange likelihood on at least one host by more than 1 in 48/60 datasets (largest gain 228.684). It beat a fresh restart from that host’s own parameters by more than 1 in 48/60. When each host retained its best witnessed BIC across candidates, 155/156 selected the same model and 154/156 gave limits within tolerance.

Each runtime started from both saved solutions for each distinct winning specification. If the specifications differed, all comparisons used the same shift specification and component count selected by the lowest BIC achieved on either machine. A reciprocal log-likelihood difference greater than 1 counted as a missed better fit. A difference of 1 or less with either endpoint more than 0.5% of the Actions width apart counted as weak identification. The 96 datasets already concordant after continuation contribute to the full-study denominator. [Dataset classifications](arm-a-cross-decisions.csv) · [All terminal solutions](arm-a-cross-solutions.csv) · [Best fit retained by each runtime](arm-a-after-retaining-best.csv)

## 4. Arm B nesting audit

| Saved datasets | Refit | Previously flagged | Flagged after refitting | Eligible nesting checks | Unavailable | Changed selected model |
| --- | --- | --- | --- | --- | --- | --- |
| N <= 5,000 | 396 | 251 | 0 | 5742 | 0 | 132 |
| N = 50,000 | 144 | 82 | 0 | 2384 | 0 | 69 |
| All | 540 | 333 | 0 | 8126 | 0 | 201 |

The 540 local refits covered 288 small-N Hb, 108 small-N Ca and 144 Hb datasets at N = 50,000. Across 8,095 fitted candidates, 0 of 8,126 eligible nesting comparisons failed. The shift or positivity constraints excluded 1,050 comparisons. All component SDs met the shared floor to numerical precision. The zero-shift likelihood and floor differences were zero across 2,695 available vectors. [All 540 dataset records](arm-b-datasets.csv) · [Arm B summary](arm-b-summary.csv) · [Candidate and nesting tables](arm-b/)

## Reproducibility

Routine-input and anchor hashes matched the archived hashes for all 156 pairs. The two runtimes also matched on 14 distinct source-file hashes, including the native likelihood and convergence recovery. The local fits used Debian 13 with reference BLAS. Actions matched the original C2 Ubuntu image, R 4.5.0, BLAS, LAPACK and package versions. The Actions run retains the raw fit parameters and trace files. [Source hashes](source-hashes.csv) · [Local donor parameters](../exchange/local/) · [Actions donor parameters](../exchange/actions/) · [Local environment](local-runtime.txt) · [Actions environment](Actions-runtime.txt) · [GitHub Actions run](https://github.com/amromeo/amromeo.github.io/actions/runs/45180631883)

The thresholds and the Arm A cross-start rule were recorded before this round. The protocol records the numerical recovery for Hb 104. The two-runtime findings describe the 156 Hb datasets. Arm B was audited locally on its saved subset. [Protocol](../../PROTOCOL.md)
