# Task D follow-up

B2 confirmation retains its committed rule a at t = 0.80 and LRT-alone comparator. Its replicates 21–40 run separately without changing that selection.

For C2 small-N fits, compare the saved lower limit with the minimum routine value. CRP has a one-sided upper reference limit, so its lower-limit flag is false. Report the strict mean absolute z > 5 group (23 fits) and < 0.25 group (1,952 fits). Also retain the either-endpoint > 5 group to reconcile the earlier 33-fit audit. Report the lower-limit flag, LRT p < 0.01, their OR and their AND. No fits are repeated for this check.

The lower-limit flag catches 9/23 implausible and 0/1952 accurate fits. LRT catches 23/23 and 224/1952, respectively. The OR catches the same counts as LRT alone. All nine lower-limit flags have LRT p < 0.01.

Run the 1,440 N = 50,000 definitions using the C2 estimators and seed 17000000 + Index. Verify all 1,440 rounded anchor samples are distinct before fitting. Generate each shard through RIbench's original generator, preserving its random sequence and hash checks. Refit all three methods because the larger-N routine datasets require new fits. Retain the C2 LRT as a diagnostic with no declines. Save native metrics, paired comparisons, extreme-error audit, per-dataset diagnostics, generation time and fit time.

The 50,000 workflow is on a separate branch to avoid rerunning C2. Report the completed stage runtime before launching the 1,440 N = 500,000 datasets. The 500,000 stage is authorized but has no automatic trigger in this workflow.
