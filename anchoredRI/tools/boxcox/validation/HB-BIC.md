# Hemoglobin model selection across runtimes

For Hb 100 the replay on the original GitHub Actions image reproduced the archived choice of identity. Identity beat log by 1.60 BIC points on that run. For Hb 148 the replay selected identity, as did the local fit. The original archive selected shifted-log. The original C2 results contain only the selected model and endpoints, so the BIC margin in the original run is unknown for both datasets.

All winning models had three components. The largest endpoint difference from the original C2 interval was 0.0000560 Hb units for the Hb 100 replay (0.00137% of archived width) and 0.167798 for the Hb 148 replay (4.209%). Only Hb 148 exceeded the 0.1% tolerance on replay.

| Hb index | Original C2 winner | Local winner | Local BIC gap | Repeat on C2 image | Repeat BIC gap |
|---|---|---|---:|---|---:|
| 100 | Identity | Log | 5.7511 | Identity | 1.6020 |
| 148 | Shifted-log | Identity | 37.2374 | Identity | 37.2378 |

## BIC for all candidate models

BIC uses the original measurement scale. Lower values favor a candidate. The bold value in each column marks the selected model.

| Transform | Components | Hb 100 local | Hb 100 repeat | Hb 148 local | Hb 148 repeat |
|---|---:|---:|---:|---:|---:|
| Identity | 1 | 6,128.4742 | 6,128.4742 | 23,523.5956 | 23,523.5956 |
| Identity | 2 | 3,542.3868 | 3,542.3868 | 17,444.9140 | 17,444.9140 |
| Identity | 3 | 3,371.7932 | **3,364.4401** | **16,947.4597** | **16,947.4594** |
| Log | 1 | 4,211.0276 | 4,211.0276 | 20,323.3039 | 20,323.3039 |
| Log | 2 | 3,418.1255 | 3,418.1255 | 17,011.7108 | 17,011.7108 |
| Log | 3 | **3,366.0421** | 3,366.0421 | 17,029.6611 | 17,029.6611 |
| Shifted-log | 1 | 4,191.1323 | 4,191.1323 | 20,285.4887 | 20,285.4887 |
| Shifted-log | 2 | 3,387.6673 | 3,387.6673 | 16,984.6972 | 16,984.6972 |
| Shifted-log | 3 | 3,395.4480 | 3,395.4480 | 16,988.7599 | 16,988.7599 |

On Hb 100, the log BIC agreed across the two current fits to four decimals. The identity fit had a lower BIC on Actions than locally. On Hb 148, the archived winner (shifted-log with three components) ranked third in both new fits, about 41.30 BIC points behind identity. The cause of the original Hb 148 choice remains unresolved.

Both new fits used the same routine vectors and independent anchors as C2. The replay verified the runner image version (Ubuntu 26.04, 20261116.202.1), R 4.5.0, BLAS/LAPACK and dependency versions. The five study source files matched the local source hashes. The new candidate tables also include candidate limits and ranks.

[Local candidates](hb-100-148-local-candidates.csv) · [Actions replay candidates](hb-100-148-c2-runtime-candidates.csv) · [Replay and original endpoints](hb-100-148-c2-runtime-vs-original.csv) · [Runtime](hb-100-148-c2-runtime.txt) · [Original C2 run](https://github.com/amromeo/amromeo.github.io/actions/runs/45096614442) · [Replay run](https://github.com/amromeo/amromeo.github.io/actions/runs/45110852352)
