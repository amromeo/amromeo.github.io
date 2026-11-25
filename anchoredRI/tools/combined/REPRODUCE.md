# Reproduce Task G

Use the study/combined branch. The workflow that produced the anchoredRI fits is pinned to commit 1be215f96c416b9ffde17d008dce5fa1e503b86a. It builds anchoredRI 0.3.1 from the published tarball on Ubuntu 24.04 with R 4.5.0. Each primary or audit artifact contains the selected ID list, session information and per-dataset records. Audit artifacts also contain the original shard’s primary guidance-panel checks.

Run the following in the directory containing this file, with R, Python, NumPy and pandas on the path. The download scripts refer to the published Stage 3 commit and the Task G Actions run. The analysis uses only the primary and audit files for that run.

```sh
python3 inputs.py
python3 verify-zips.py
python3 collect.py
bash run-analysis.sh
```

The headline table places RIbench and confirmation on adjacent rows for each comparison. The report displays them in adjacent columns. In the detailed tables, a missing returned-only error denotes zero returns. The loss denominator always includes every assigned dataset.
