# Compute used by interrupted jobs

Runner time comes from GitHub’s job start and completion timestamps. CPU is the sum of the user and system seconds recorded by completed timed commands. A failed command can exit before its timing receipt is written. Its unrecorded CPU time cannot be recovered from the runner duration.

| Run | Job | Occupied runner seconds | Recorded CPU seconds, lower bound | Completed timing receipts |
|---|---|---:|---:|---|
| 45854626824 | build-single | 27.0 | 2.213 | build |
| 45854822638 | pilot-4384 | 31.0 | 0.300 | unpack, prepare |
| 45854822638 | pilot-4968 | 37.0 | 0.274 | unpack, prepare |
| 45854822638 | pilot-5017 | 32.0 | 0.256 | unpack, prepare |
| 45854822638 | pilot-5510 | 35.0 | 0.276 | unpack, prepare |
| 45854822638 | pilot-3279 | 28.0 | 0.261 | unpack, prepare |
| 45854822638 | pilot-2734 | 27.0 | 0.384 | unpack, prepare |

The 7 interrupted jobs occupied runners for 0.0603 hours. They have 0.0011 recorded CPU-hours plus unmeasured CPU. Total study CPU is therefore a lower bound. Remaining costs use the 20 completed pilot records. Failed jobs enter spent cost once and do not affect per-stratum fitting times or pilot setup averages.

[Job ledger](job-ledger.csv) · [Projection input](jobs-final.csv)
