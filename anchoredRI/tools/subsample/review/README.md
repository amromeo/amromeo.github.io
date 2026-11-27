# Task H review of the combined release rule

The revised policy uses the subsample median when the full-data first stage passes. When the first stage declines, the package fallback takes priority. The subsample majority can release a median only after the package also withholds. These results use the saved cap-5 records for all 1,440 RIbench datasets at N = 50,000.

With seed set A, the policy returned 1,436 intervals. Mean endpoint error on returns was 0.1316, and the all-assignment score was 0.1451. Seed set B returned 1,435, with mean error 0.1317 and all-assignment score 0.1486. Neither bank released an extreme endpoint. Both preserved the package’s 1,423 returns. On those same datasets, mean error fell from 0.1476 to 0.1308 for A and 0.1311 for B. Step c added 13 of 17 package nonreturns under A and 12 under B. [Routes, paired comparisons and per-dataset results](policy/REPORT.md)

The 934 full-data fits selecting five components underestimated the healthy share by an average of 18.639 percentage points. Integrating the generating healthy distribution and the fitted first component at the true reference limits placed 18.296 points of this net deficit inside and 0.343 points outside. As a fraction of true healthy mass in each region, the shortfall was 27.8% inside and 9.9% outside. The evidence places most of the shortfall within the reference interval. It leaves open whether it concentrates near the endpoints within that interval. [Tail decomposition, rounding sensitivity and per-dataset results](tail/REPORT.md)

Task H closed on 27 November 2026. The revised policy failed independent confirmation on anchorRefine’s Stage 2 synthetic datasets at N = 50,000. Arm W and the three-component variants concluded with the RIbench analysis. [Confirmation report](../results/confirmation-45912200495/REPORT.md) · [RIbench results](../results/full-45860591866/REPORT.md) · [Study summary](../SUPPLEMENT.md)

[SHA-256 checksums for the review files](SHA256SUMS)
