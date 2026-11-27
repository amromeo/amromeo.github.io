# Task H. Subsampling large routine datasets

**Closed on 27 November 2026 after failed independent confirmation.** On RIbench, the subsampling policy improved mean endpoint accuracy and returned 13 additional intervals. On the 1,440 independent synthetic datasets, both methods returned intervals for the same 1,437 datasets, and subsampling increased mean endpoint error. Both subsample seed sets gave the same decision. anchoredRI retains its existing full-data method and fallback. [Supplement summary](SUPPLEMENT.md) · [Confirmation report](results/confirmation-45912200495/REPORT.md)

| Study record | Protocol and results | Archived data |
|---|---|---|
| RIbench at N = 50,000 | [Study protocol](PROTOCOL.md) · [Full benchmark](results/full-45860591866/REPORT.md) · [Revised policy](review/policy/REPORT.md) | [RIbench data release](https://github.com/amromeo/amromeo.github.io/releases/tag/task-h-subsample-50000) |
| Independent confirmation at N = 50,000 | [Confirmation protocol](confirmation/PROTOCOL.md) · [Results and confidence intervals](results/confirmation-45912200495/REPORT.md) | [Confirmation data release](https://github.com/amromeo/amromeo.github.io/releases/tag/task-h-confirmation-50000) |

A separate analysis of archived Task E Arm B fits examined how often an extra mixture component was centered in the healthy range and how this related to first-stage error. It supports a hypothesis about how a mixture can split a healthy population. Independent confirmation of that explanation is outstanding. [Component report and safeguard table](review/component-centers/REPORT.md)
