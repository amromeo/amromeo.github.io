# Reproduce Step B

Run from anchoredRI/ using R (base packages only) and Python with pandas/numpy:

```sh
Rscript tools/disagreement/run.R test
seq 1 60 | xargs -P 8 -I '{}' Rscript tools/disagreement/run.R '{}'
Rscript tools/disagreement/run.R
python tools/disagreement/summarize.py
```

The 60 bootstrap jobs are independent. Each writes its own RDS checkpoint. The final R command validates matching anchors and scores the existing 900 anchoredRI estimates. The Python command checks withholding behavior, writes summaries, records source hashes, and packages the results. No production R code is changed.

The vectorized bootstrap profile uses exactly the same 1,001-point grid and local refinement as the Step A estimator. The test compares full selection results with that estimator on 60 bootstrap samples and checks endpoint boundaries, unavailable checks, and the likelihood ratio at the anchors' own fitted parameters.

`datasets.csv` separates diagnostic estimates from reported intervals. A declined result has missing reported limits, with a recommendation. Secondary declines are hypothetical and do not control the primary reported interval.
