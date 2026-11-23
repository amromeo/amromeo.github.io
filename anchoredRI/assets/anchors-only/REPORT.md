# Anchors-only comparator results

Completed 45 scenarios × 20 paired datasets. Routine data met the specified rule in 28 of 30 separated/moderate cells.

Negative paired differences favor anchoredRI. The rule is mean difference + 2 × MCSE < 0, applied separately to each separated or moderate cell. Heavy-overlap cells are descriptive only.

| Family | Separated: cells meeting rule | Moderate: cells meeting rule |
|---|---:|---:|
| normal | 5/5 | 5/5 |
| lognormal | 5/5 | 5/5 |
| shifted | 4/5 | 4/5 |

Each table cell reports **paired mean difference ± 2 MCSE; anchoredRI wins / 20; decision**. Differences are percentage points of error normalized by the true reference interval width. The ±2 MCSE display describes simulation precision; it is not a clinical confidence interval.

## normal

| Overlap | 10% affected | 30% affected | 50% affected | 70% affected | 90% affected |
|---|---|---|---|---|---|
| separated | -5.54 ± 1.29; 20/20; adds value | -4.63 ± 1.62; 19/20; adds value | -4.70 ± 1.58; 19/20; adds value | -4.71 ± 1.53; 18/20; adds value | -3.74 ± 1.47; 17/20; adds value |
| moderate | -4.89 ± 1.33; 20/20; adds value | -3.43 ± 1.48; 18/20; adds value | -2.81 ± 1.86; 16/20; adds value | -3.95 ± 1.04; 19/20; adds value | -1.56 ± 1.53; 12/20; adds value |
| heavy | -2.79 ± 1.72; 16/20; descriptive only | +2.67 ± 3.12; 8/20; descriptive only | +6.92 ± 4.38; 8/20; descriptive only | +9.59 ± 4.72; 6/20; descriptive only | +19.61 ± 4.10; 2/20; descriptive only |

## lognormal

| Overlap | 10% affected | 30% affected | 50% affected | 70% affected | 90% affected |
|---|---|---|---|---|---|
| separated | -14.76 ± 4.47; 19/20; adds value | -14.21 ± 4.74; 19/20; adds value | -13.75 ± 4.45; 19/20; adds value | -12.13 ± 4.76; 19/20; adds value | -9.87 ± 4.55; 18/20; adds value |
| moderate | -13.24 ± 4.44; 20/20; adds value | -10.37 ± 5.11; 16/20; adds value | -9.23 ± 4.43; 16/20; adds value | -6.85 ± 4.42; 13/20; adds value | -5.39 ± 3.97; 15/20; adds value |
| heavy | -6.33 ± 4.77; 15/20; descriptive only | +6.17 ± 6.59; 5/20; descriptive only | +26.47 ± 6.52; 1/20; descriptive only | +33.93 ± 8.37; 1/20; descriptive only | +40.58 ± 11.79; 3/20; descriptive only |

## shifted

| Overlap | 10% affected | 30% affected | 50% affected | 70% affected | 90% affected |
|---|---|---|---|---|---|
| separated | -13.31 ± 4.06; 20/20; adds value | -13.17 ± 3.75; 19/20; adds value | -11.45 ± 4.30; 17/20; adds value | -9.51 ± 4.14; 17/20; adds value | -4.50 ± 5.17; 16/20; not demonstrated |
| moderate | -12.24 ± 3.78; 18/20; adds value | -8.35 ± 4.54; 16/20; adds value | -10.60 ± 3.84; 18/20; adds value | -8.48 ± 3.60; 17/20; adds value | -1.18 ± 4.25; 13/20; not demonstrated |
| heavy | -3.11 ± 4.95; 11/20; descriptive only | +6.10 ± 8.52; 9/20; descriptive only | +35.81 ± 8.72; 2/20; descriptive only | +59.14 ± 6.81; 0/20; descriptive only | +72.40 ± 8.12; 0/20; descriptive only |

## Anchors-only context

These means use the 20 distinct anchor samples for each family, without counting their reuse across scenarios as additional independent samples.

| Family | Auto mean error | Oracle mean error | Auto scale choices: identity / log / shifted-log | Auto offset bound hits |
|---|---:|---:|---|---:|
| normal | 6.44% | 5.65% | 9 / 11 / 0 | 0/20 |
| lognormal | 15.87% | 11.58% | 1 / 16 / 3 | 0/20 |
| shifted | 15.57% | 12.72% | 0 / 12 / 8 | 0/20 |

Failed anchors-only estimates: 0 of 1,800 evaluated rows. Failed archived anchoredRI estimates: 0 of 900.

The oracle knows the family and true shifted-log offset, but estimates its mean and sample SD from the same 40 anchors. Auto uses maximum-likelihood SD for BIC and sample SD for limits, exactly as specified in the protocol. Both use 1.96.

No changes were made to R/ or tools/stress-design.R. AnchoredRI estimates were reused from the frozen completed comparison. This is a follow-up on the existing grid, not independent validation. The 30 primary decisions are not a multiplicity-adjusted global superiority test.
