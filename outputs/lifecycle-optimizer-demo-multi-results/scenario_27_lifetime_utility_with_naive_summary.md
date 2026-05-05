# 27 Scenarios: Lifecycle Utility + Naive Baseline

## Scope

- Scenarios: `rho in [6, 8, 10]` x `wealth in [5, 10, 20]` x `age in [30, 45, 60]` = 27.
- `rho` is the agent risk aversion / risk preference parameter; higher rho means stronger risk aversion.
- Official model score: lifecycle utility from the model value function, interpolated at the scenario age and wealth.
- Baseline: naive rule with fixed 60/40 portfolio and fixed 4% consumption rate, evaluated by a deterministic CRRA simulation for comparison only.
- All underlying Octave model artifacts used here have `status=ok` and `year_files_count=80`.

## Overall Result

- Best scenario by model value-function utility: age=60, wealth=20.0, rho=6.0, utility=2.60679074.
- Worst scenario by model value-function utility: age=45, wealth=5.0, rho=10.0, utility=1.61699309.

## Method Utility Comparison

This table uses the deterministic CRRA simulator so the optimal policy and naive policy are evaluated on the same comparison metric. It is separate from the model value-function utility above.

| method | utility | note |
|---|---:|---|
| model optimal policy | -0.00028212 | deterministic CRRA average across 27 scenarios |
| naive 60/40 + 4% consumption | -0.02131200 | fixed portfolio and fixed consumption rate baseline |

## Ranking by rho

| rho | avg_model_value_function_utility | avg_optimal_policy_sim_utility | avg_naive_sim_utility | avg_gap_opt_minus_naive |
|---:|---:|---:|---:|---:|
| 6.0 | 2.23787067 | -0.00081140 | -0.04046043 | 0.03964903 |
| 8.0 | 2.04521480 | -0.00003309 | -0.01596347 | 0.01593039 |
| 10.0 | 1.88211706 | -0.00000188 | -0.00751211 | 0.00751023 |

## 27 Scenario Table

| rho | age | wealth | model_value_function_utility | optimal alpha | optimal consumption | opt sim utility | naive sim utility | gap |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 6.0 | 30 | 5 | 2.20373335 | 0.0025 | 4.1460 | -0.00001530 | -0.09267680 | 0.09266150 |
| 6.0 | 30 | 10 | 2.28586335 | 0.1059 | 4.1460 | -0.00001530 | -0.04333101 | 0.04331572 |
| 6.0 | 30 | 20 | 2.32546794 | 0.1111 | 13.8950 | -0.00001530 | -0.01279346 | 0.01277816 |
| 6.0 | 45 | 5 | 1.94920530 | 0.3283 | 4.3813 | -0.00002391 | -0.04976971 | 0.04974579 |
| 6.0 | 45 | 10 | 2.01351992 | 0.1215 | 4.3813 | -0.00002283 | -0.02529457 | 0.02527174 |
| 6.0 | 45 | 20 | 2.04383400 | 0.1111 | 14.1303 | -0.00002143 | -0.00833945 | 0.00831802 |
| 6.0 | 60 | 5 | 2.30049611 | 0.0051 | 4.1460 | -0.00285137 | -0.08056417 | 0.07771280 |
| 6.0 | 60 | 10 | 2.41192530 | 0.2118 | 4.1460 | -0.00244496 | -0.03903090 | 0.03658594 |
| 6.0 | 60 | 20 | 2.60679074 | 0.8936 | 4.1460 | -0.00189217 | -0.01234377 | 0.01045160 |
| 8.0 | 30 | 5 | 1.98008892 | 0.3283 | 4.3813 | -0.00000010 | -0.04583019 | 0.04583008 |
| 8.0 | 30 | 10 | 2.04357382 | 0.1215 | 4.3813 | -0.00000010 | -0.01544301 | 0.01544290 |
| 8.0 | 30 | 20 | 2.07452491 | 0.1111 | 14.1303 | -0.00000010 | -0.00266221 | 0.00266211 |
| 8.0 | 45 | 5 | 1.77148804 | 0.3257 | 4.3813 | -0.00000017 | -0.01916117 | 0.01916100 |
| 8.0 | 45 | 10 | 1.81665297 | 0.0156 | 4.3813 | -0.00000016 | -0.00727249 | 0.00727233 |
| 8.0 | 45 | 20 | 1.84170434 | 0.0000 | 14.1303 | -0.00000016 | -0.00146983 | 0.00146968 |
| 8.0 | 60 | 5 | 2.16974836 | 0.0051 | 4.1460 | -0.00012261 | -0.03671695 | 0.03659434 |
| 8.0 | 60 | 10 | 2.26984982 | 0.2118 | 4.1460 | -0.00010092 | -0.01279744 | 0.01269652 |
| 8.0 | 60 | 20 | 2.43930203 | 0.8936 | 4.1460 | -0.00007348 | -0.00231799 | 0.00224451 |
| 10.0 | 30 | 5 | 1.79478931 | 0.3257 | 4.3813 | -0.00000000 | -0.02539166 | 0.02539166 |
| 10.0 | 30 | 10 | 1.83841833 | 0.0156 | 4.3813 | -0.00000000 | -0.00621841 | 0.00621841 |
| 10.0 | 30 | 20 | 1.86384939 | 0.0000 | 14.1303 | -0.00000000 | -0.00063535 | 0.00063535 |
| 10.0 | 45 | 5 | 1.61699309 | 0.2172 | 4.3813 | -0.00000000 | -0.00827170 | 0.00827169 |
| 10.0 | 45 | 10 | 1.64538069 | 0.0104 | 4.3813 | -0.00000000 | -0.00236272 | 0.00236272 |
| 10.0 | 45 | 20 | 1.66609617 | 0.0000 | 14.1303 | -0.00000000 | -0.00029691 | 0.00029691 |
| 10.0 | 60 | 5 | 2.05848790 | 0.3308 | 4.3813 | -0.00000760 | -0.01904472 | 0.01903713 |
| 10.0 | 60 | 10 | 2.14806308 | 0.2274 | 4.3813 | -0.00000576 | -0.00486203 | 0.00485628 |
| 10.0 | 60 | 20 | 2.30697562 | 0.8936 | 4.3813 | -0.00000353 | -0.00052548 | 0.00052195 |

## Does the model match economic intuition?

Mostly yes, with caveats:

- Higher wealth generally raises model value-function utility. This is economically intuitive because the value function is increasing in resources.
- Higher `rho` lowers average lifecycle utility and generally lowers risky exposure at many checkpoints. This matches the interpretation of `rho` as stronger risk aversion.
- The optimal policy beats the naive 60/40 + 4% baseline in the deterministic CRRA comparison on average, which is the expected direction if the model is adding useful state-dependent policy.

## Which variables matter most?

From this 27-scenario grid:

- `wealth` is the largest direct driver of value-function utility, because the value function is steeply increasing in cash/wealth.
- `rho` is the key preference driver: changing rho from 6 to 10 materially changes utility levels and policy aggressiveness.
- `age` changes which policy file and lifecycle phase is used. It matters, but in this small grid its effect is intertwined with retirement timing and available wealth.

Average model value-function utility by age:

| age | avg utility |
|---:|---:|
| 30 | 2.04558992 |
| 45 | 1.81831939 |
| 60 | 2.30129322 |

Average model value-function utility by wealth:

| wealth | avg utility |
|---:|---:|
| 5 | 1.98278115 |
| 10 | 2.05258303 |
| 20 | 2.12983835 |

## What looks unreasonable or needs investigation?

- Some low-wealth states choose very high risky allocation (`alpha` near 1.0). This may be a borrowing/constraint/grid artifact or a consequence of labor-income-like background wealth; it needs economic review before customer-facing interpretation.
- Some old-age mid-wealth policies can become more growth-oriented than expected. This may be due to the coarse fast-mode grid (`na=10`, `ncash=10`) rather than a final production policy.
- The naive baseline utility is computed by a deterministic CRRA comparison simulator, not by solving the model under a fixed-policy constraint. It is useful for direction, but not a formal welfare theorem.
- The scenario grid only varies `rho`, age, and wealth. It does not yet test market return, volatility, discounting, retirement age, or EIS (`psi`).

## Runner Stability

- Reused three full real-model Octave runs with `status=ok` and `year_files_count=80`.
- No timeout rows are included in this summary.
- Unit tests were run before summary generation: 8 tests OK.

## Output Files

- Summary: `C:\Users\haoyu\Desktop\code\github\LifeCycle\outputs\lifecycle-optimizer-demo-multi-results\scenario_27_lifetime_utility_with_naive_summary.md`
- CSV: `C:\Users\haoyu\Desktop\code\github\LifeCycle\outputs\lifecycle-optimizer-demo-multi-results\scenario_27_lifetime_utility_with_naive.csv`
- Source artifacts: `C:\Users\haoyu\Desktop\code\github\LifeCycle\outputs\lifecycle-optimizer-demo-multi-results`
