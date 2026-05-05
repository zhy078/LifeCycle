# Grid and Stochastic Validation

## Step 1: Is the weird alpha a grid problem?

Target state: `age=60`, `wealth=20`, `rho in [6,8,10]`.

| grid | rho | alpha | lifecycle utility | status | run_seconds |
|---|---:|---:|---:|---|---:|
| 10x10 | 6 | 0.8936 | 2.606791 | ok | 173.4 |
| 10x10 | 8 | 0.8936 | 2.439302 | ok | 173.7 |
| 10x10 | 10 | 0.8936 | 2.306976 | ok | 174.1 |
| 20x20 | 6 | 0.4520 | 2.794249 | ok | 741.1 |
| 20x20 | 8 | 0.3467 | 2.591854 | ok | 725.9 |
| 20x20 | 10 | 0.2632 | 2.470377 | ok | 728.1 |
| 40x40 | 6 | 0.4476 | 2.796849 | ok | 3214.6 |
| 40x40 | 8 | 0.3450 | 2.594380 | ok | 3090.9 |
| 40x40 | 10 | 0.2681 | 2.473277 | ok | 3191.9 |

Conclusion:

- The 10x10 grid produced an extreme alpha around 0.8936 at age=60, wealth=20.
- The 20x20 and 40x40 grids reduce alpha sharply into a much more moderate range.
- This strongly supports the hypothesis that the original extreme alpha was mainly a coarse-grid artifact.

## Step 2: Age rebound validation

Fixed state: `rho=8`, `wealth=10`.

| grid | age | alpha | lifecycle utility |
|---|---:|---:|---:|
| 10x10 | 30 | 0.1215 | 2.043574 |
| 10x10 | 35 | 0.1215 | 1.896824 |
| 10x10 | 40 | 0.0156 | 1.823222 |
| 10x10 | 45 | 0.0156 | 1.816653 |
| 10x10 | 50 | 0.1215 | 1.878778 |
| 10x10 | 55 | 0.1215 | 2.020938 |
| 10x10 | 60 | 0.2118 | 2.269850 |
| 20x20 | 30 | 0.5325 | 1.188117 |
| 20x20 | 35 | 0.4799 | 1.182445 |
| 20x20 | 40 | 0.4799 | 1.221699 |
| 20x20 | 45 | 0.4799 | 1.305847 |
| 20x20 | 50 | 0.4799 | 1.439983 |
| 20x20 | 55 | 0.4799 | 1.634533 |
| 20x20 | 60 | 0.4799 | 1.904383 |
| 40x40 | 30 | 0.5128 | 1.189490 |
| 40x40 | 35 | 0.4871 | 1.183836 |
| 40x40 | 40 | 0.4615 | 1.223124 |
| 40x40 | 45 | 0.4615 | 1.307383 |
| 40x40 | 50 | 0.4615 | 1.441744 |
| 40x40 | 55 | 0.4615 | 1.636669 |
| 40x40 | 60 | 0.4871 | 1.907484 |

Conclusion:

- In the 10x10 grid, alpha falls from age 30 to 45 and then rebounds by age 60.
- Finer grids smooth part of the pattern, but the age profile is still not perfectly monotone.
- This is no longer just a single bad point; it should be treated as a structural policy-shape diagnostic around the retirement transition.

## Step 3: Stochastic validation

Simulation setup:

- 1000 Monte Carlo paths per scenario and method.
- Risky return: `R_t = r + mu + sigma * eps_t`, `eps_t ~ N(0,1)`.
- Compare `optimal_policy` from the 40x40 model against `naive_60_40_4pct`.

- optimal_policy: avg utility_mean=-0.00243375, avg terminal wealth=240.5065, avg low-wealth probability=0.0000
- naive_60_40_4pct: avg utility_mean=-0.02273486, avg terminal wealth=379.6954, avg low-wealth probability=0.0000

| rho | age | wealth | optimal utility mean | naive utility mean | utility gap | optimal terminal wealth | naive terminal wealth |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 6 | 30 | 5 | -0.01118701 | -0.09286492 | 0.08167790 | 733.3647 | 623.7068 |
| 6 | 30 | 10 | -0.00763352 | -0.04351931 | 0.03588579 | 739.5630 | 628.2611 |
| 6 | 30 | 20 | -0.00327837 | -0.01291360 | 0.00963523 | 775.4218 | 628.3590 |
| 6 | 45 | 5 | -0.00814369 | -0.04991722 | 0.04177352 | 34.5619 | 391.0781 |
| 6 | 45 | 10 | -0.00535143 | -0.02542328 | 0.02007186 | 37.7300 | 387.0452 |
| 6 | 45 | 20 | -0.00257552 | -0.00842104 | 0.00584553 | 41.3496 | 388.1729 |
| 6 | 60 | 5 | -0.00655731 | -0.08930200 | 0.08274469 | 3.0657 | 116.3206 |
| 6 | 60 | 10 | -0.00514760 | -0.04861977 | 0.04347218 | 3.1368 | 119.8907 |
| 6 | 60 | 20 | -0.00250676 | -0.01854453 | 0.01603777 | 3.4301 | 126.2508 |
| 8 | 30 | 5 | -0.00366394 | -0.04587710 | 0.04221317 | 677.0992 | 637.5013 |
| 8 | 30 | 10 | -0.00158698 | -0.01547054 | 0.01388356 | 686.8085 | 626.3207 |
| 8 | 30 | 20 | -0.00059516 | -0.00267623 | 0.00208107 | 706.7898 | 634.7838 |
| 8 | 45 | 5 | -0.00177400 | -0.01918187 | 0.01740787 | 33.3556 | 382.1410 |
| 8 | 45 | 10 | -0.00091101 | -0.00729003 | 0.00637902 | 34.9449 | 381.3977 |
| 8 | 45 | 20 | -0.00028092 | -0.00147686 | 0.00119594 | 38.8143 | 392.2027 |
| 8 | 60 | 5 | -0.00091393 | -0.03990306 | 0.03898913 | 3.1332 | 114.6858 |
| 8 | 60 | 10 | -0.00061977 | -0.01748382 | 0.01686405 | 3.2674 | 120.1578 |
| 8 | 60 | 20 | -0.00018303 | -0.00350092 | 0.00331789 | 3.4857 | 122.9072 |
| 10 | 30 | 5 | -0.00106843 | -0.02540361 | 0.02433518 | 598.7001 | 628.0496 |
| 10 | 30 | 10 | -0.00050052 | -0.00622587 | 0.00572535 | 617.2468 | 644.4234 |
| 10 | 30 | 20 | -0.00007949 | -0.00063725 | 0.00055776 | 612.8263 | 640.3411 |
| 10 | 45 | 5 | -0.00054046 | -0.00827585 | 0.00773539 | 29.4624 | 371.4445 |
| 10 | 45 | 10 | -0.00022089 | -0.00236505 | 0.00214416 | 31.4720 | 388.3590 |
| 10 | 45 | 20 | -0.00004537 | -0.00029798 | 0.00025261 | 34.1761 | 400.3739 |
| 10 | 60 | 5 | -0.00023790 | -0.02125137 | 0.02101348 | 3.2954 | 113.1022 |
| 10 | 60 | 10 | -0.00008374 | -0.00570911 | 0.00562537 | 3.4590 | 115.6517 |
| 10 | 60 | 20 | -0.00002459 | -0.00128900 | 0.00126441 | 3.7146 | 128.8467 |

## Runner stability

- 20x20: all 3 scenarios completed with `status=ok`.
- 40x40: all 3 scenarios completed with `status=ok`.
- No timeout rows were used in this analysis.

## Output files

- This report: `C:\Users\haoyu\Desktop\code\github\LifeCycle\outputs\grid-validation-analysis\grid_and_stochastic_validation.md`
- Grid comparison CSV: `C:\Users\haoyu\Desktop\code\github\LifeCycle\outputs\grid-validation-analysis\grid_refinement_comparison.csv`
- Age slice CSV: `C:\Users\haoyu\Desktop\code\github\LifeCycle\outputs\grid-validation-analysis\age_rebound_slice.csv`
- Stochastic validation CSV: `C:\Users\haoyu\Desktop\code\github\LifeCycle\outputs\grid-validation-analysis\stochastic_validation_1000_paths.csv`
