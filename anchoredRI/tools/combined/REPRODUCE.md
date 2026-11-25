# Reproduce Task G

Use the study/combined branch. The workflow that produced the anchoredRI fits is pinned to commit 5862f650ac49a657951fa3cbb0430a81249ed65a. All jobs use the same R 4.5.0 container and the same compiled library. See [the computational audit](RUN-NOTE.md) for the first run and reason for repeating it. Each primary or audit artifact contains the selected ID list, session information and per-dataset records. Audit artifacts also contain the original shard’s primary guidance-panel checks.

Run the following in the directory containing this file, with R, Python, NumPy and pandas on the path. The download scripts refer to the published Stage 3 commit and the Task G Actions run. The analysis uses only the primary and audit files for that run.

```sh
python3 inputs.py
python3 verify-zips.py
python3 collect.py
bash run-analysis.sh
```

The headline table places RIbench and confirmation on adjacent rows for each comparison. The report displays them in adjacent columns. In the detailed tables, a missing returned-only error denotes zero returns. The loss denominator always includes every assigned dataset.
