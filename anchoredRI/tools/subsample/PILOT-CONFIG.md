# Task H timing pilot

The sampling frame contains 1,440 datasets at nominal N = 50,000. The 20 pilot IDs were chosen from design metadata before fitting. Each analyte contributes one dataset with at most 20% pathology and one above 20%. Each of the 20 strata contains 72 datasets. The roster records the observed routine length, original Task G shard and both ordered vector hashes. Four pilot datasets have 49,999 routine observations. [Pilot roster](pilot-roster.csv)

Within each stratum, sort the lowercase hexadecimal SHA-256 of the UTF-8 string `TaskH-pilot-v1|<Index>` and take the first entry, breaking ties by smaller numeric Index. The decimal Index has no padding or newline. Roster SHA-256 is `61f6d0c8aa46aaf55b8d2d19d88b7399ae68a1f923a096a060b3f84ad40d091d`.

| Execution setting | Frozen configuration |
| --- | --- |
| Machine | Standard GitHub-hosted Linux x64 runner, `ubuntu-24.04` |
| Container | `rocker/r-ver@sha256:5e8963d73b87865e0ad779df4202555fb1b436ee95dc7d7442450308b8ea77f3` |
| R and libraries | R 4.5.0 and the archived Stage 0 library bundle, SHA-256 `ab482e6a48109c476c8d4e6d14231e64f28470259708cdba728b7b8ec0e91027`. The study freeze covers the library versions, package source and installed-file checks. |
| Archive verification | 48 jobs, one per original Task G shard of 30 datasets, with at most 20 active jobs |
| Timing pilot | 20 jobs, one per published Index, with at most 20 active jobs. Each fitting process uses one numerical thread. |
| Projection | Each pilot dataset represents 72 datasets in its own stratum. Measure setup and reporting separately, account for shared candidate work once and deduct the 20 completed pilot datasets from the remaining work. |

## Input verification

From the repository root, run the following to verify the published roster or to materialize the files for one pilot job.

```sh
python3 anchoredRI/tools/subsample/prepare-inputs.py /tmp/task-h-input --metadata-only
python3 anchoredRI/tools/subsample/prepare-inputs.py /tmp/task-h-input --ids 2065
```

The verifier reads immutable revisions from GitHub or accepts local object stores through `--anchorrefine-repo PATH` and `--site-repo PATH`. The two immutable revisions and file digests are in the restoration script. Task G's manifest covers the dataset index, archived fit and check tables, and all 48 original RDS members. The pinned outer archive digest list covers every ZIP before extraction. The restoration script joins the design and Task G records by original Index. The fitting process reads the original case by shard position, verifies the two ordered-vector SHA-256 hashes using R serialization version 2, and checks the observed length, rounding increment and requested probabilities. [Restoration script](prepare-inputs.py) · [R reader and position sampler](inputs.R)

The sampler uses sorted positions from the observed routine vector. Its two seed sets are independent of fitting and reuse the same archived panel of 40 anchors. R's RNG state is restored on return so sampling does not consume draws from the fitter's RNG stream. [Approved protocol](PROTOCOL.md)
