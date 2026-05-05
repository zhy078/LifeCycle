# 100-Client Lifecycle vs Naive 60/40 Comparison

## Setup

- Clients: 100 random clients
- Simulation: 1000 Monte Carlo paths per client, 30 years
- Lifecycle policy: alpha interpolated from completed 20x20 lifecycle skill artifacts
- Naive trading baseline: fixed 60% risky / 40% stabilizing assets
- Consumption assumption: fixed 4% of initial wealth per year
- Utility metric: simulated CRRA utility using each client's mapped rho and a 25% minimum-consumption floor after depletion; this is a product-layer comparison metric, not the model value function itself.
- Detail CSV: `C:\Users\haoyu\Desktop\code\github\LifeCycle\outputs\client-batch-comparison\client_100_comparison.csv`

## Aggregate Result

| Metric | Lifecycle | Naive 60/40 | Lift (Lifecycle - Naive) |
|---|---:|---:|---:|
| Utility win rate | 55.0% | 45.0% | 55/100 clients |
| Failure probability | 31.55% | 35.68% | -4.12 pp |
| Max drawdown | -20.01% | -26.84% | 6.83 pp |
| Annual volatility | 7.07% | 8.99% | -1.92 pp |
| Terminal wealth | 3.5627 | 4.8080 | -1.2452 |

## Win Counts

- Utility higher: 55/100 clients
- Failure probability lower: 76/100 clients
- Drawdown less severe: 82/100 clients
- Terminal wealth higher: 18/100 clients

## By Risk Preference

| Risk preference | Clients | Avg lifecycle alpha | Utility win rate | Failure lift | Volatility lift |
|---|---:|---:|---:|---:|---:|
| conservative | 45 | 0.348 | 57.8% | -8.24 pp | -3.78 pp |
| balanced | 21 | 0.482 | 57.1% | -3.08 pp | -1.76 pp |
| aggressive | 34 | 0.629 | 50.0% | 0.68 pp | 0.44 pp |

## Interpretation

- Positive utility lift means the lifecycle policy improved simulated client utility versus naive 60/40 under the same return shocks.
- Negative volatility lift or positive max-drawdown lift means the lifecycle policy reduced risk versus naive 60/40.
- Terminal wealth is reported separately because lifecycle utility is not terminal-wealth maximization.
- This is a product-layer stochastic validation, not a replacement for solving a fixed-policy naive baseline inside the dynamic model.
