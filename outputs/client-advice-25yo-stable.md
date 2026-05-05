# Lifecycle Client Advice

## Parsed Client Input

- Age: 25
- Wealth: 0.10 model units
- Risk preference: conservative -> rho=10

## Model Output

- Model status: ok
- Matched rho: 10.0
- Raw model risky allocation alpha: 1.0000
- Recommended consumption: 0.2497
- Lifecycle utility: 0.7044

## Risk Control Decision

- Case type: UNSAFE_MODEL_OUTPUT
- Decision status: adjusted
- Suggested allocation range: 0%-10% risky assets
- Point estimate alpha: 0.0500
- Failure probability limit: 5.0%
- Risk horizon: 30 years (age 25-55; long-term accumulation planning window)

## Decision Trace

| Step | Before | After | Reason |
|---|---:|---:|---|
| failure_risk | 100% | 5% | Reduced until estimated failure probability remains within 5.0% under Monte Carlo uncertainty bounds. |
| wealth_constraint | 5% | 5% | Applied low-wealth protection cap. |

The final point estimate remains at 5% because the failure-risk search is already the binding risk control.

## Failure Risk Search

| Candidate alpha | Failure estimate | 95% MC margin | Upper bound | Result |
|---:|---:|---:|---:|---|
| 100% | 43.2% | +/- 3.1 pp | 46.3% | Reduce |
| 95% | 41.1% | +/- 3.0 pp | 44.1% | Reduce |
| 90% | 41.5% | +/- 3.1 pp | 44.6% | Reduce |
| 85% | 41.1% | +/- 3.0 pp | 44.1% | Reduce |
| 80% | 37.4% | +/- 3.0 pp | 40.4% | Reduce |
| 75% | 37.3% | +/- 3.0 pp | 40.3% | Reduce |
| 70% | 34.8% | +/- 3.0 pp | 37.8% | Reduce |
| 65% | 35.9% | +/- 3.0 pp | 38.9% | Reduce |
| 60% | 34.7% | +/- 3.0 pp | 37.7% | Reduce |
| 55% | 35.3% | +/- 3.0 pp | 38.3% | Reduce |
| 50% | 33.2% | +/- 2.9 pp | 36.1% | Reduce |
| 45% | 31.9% | +/- 2.9 pp | 34.8% | Reduce |
| 40% | 28.3% | +/- 2.8 pp | 31.1% | Reduce |
| 35% | 27.0% | +/- 2.8 pp | 29.8% | Reduce |
| 30% | 23.1% | +/- 2.6 pp | 25.7% | Reduce |
| 25% | 19.3% | +/- 2.4 pp | 21.7% | Reduce |
| 20% | 17.6% | +/- 2.4 pp | 20.0% | Reduce |
| 15% | 15.4% | +/- 2.2 pp | 17.6% | Reduce |
| 10% | 9.1% | +/- 1.8 pp | 10.9% | Reduce |
| 5% | 1.6% | +/- 0.8 pp | 2.4% | Pass |

## Risk Metrics

| Metric | Raw model | Final recommendation |
|---|---:|---:|
| Max drawdown | -44.3% | -0.1% |
| Failure probability | 43.2% (95% MC margin +/- 3.1 pp) | 1.6% (95% MC margin +/- 0.8 pp) |
| Annual volatility | 14.9% | 0.7% |

## Baseline Comparison

| Metric | Standard 60/40 | Final recommendation |
|---|---:|---:|
| Risky allocation | 60% | 5% |
| Max drawdown | -26.5% | -0.1% |
| Failure probability | 34.7% (95% MC margin +/- 3.0 pp) | 1.6% (95% MC margin +/- 0.8 pp) |
| Annual volatility | 9.0% | 0.7% |
| Expected annual return | 2.4% | 1.6% |

## Applicability

This recommendation framework is suitable for early-career clients with limited wealth buffers, where building financial resilience should come before taking large investment risk.

## Customer Explanation

Model output adjusted due to risk controls.

The original model suggests a high-risk allocation despite the client's conservative risk preference. This recommendation is adjusted to ensure:

- Failure probability remains below the product risk limit
- Risk exposure is appropriate for the client's lifecycle stage
- Wealth level is protected from large drawdowns

Adjustment reasons:
- Failure probability 43.2% exceeds the 5.0% limit.
- Low wealth requires extra protection from large drawdowns.

Suggested allocation range: 0%-10% risky assets.
Point estimate: 5% risky assets, 95% stabilizing assets. This balances growth and capital preservation.

For a young client with a very small wealth base, building a safety buffer is typically more important than maximizing near-term investment returns.

Compared with a standard 60/40 portfolio, the adjusted recommendation lowers drawdown, failure probability, and volatility, with a modestly lower expected annual return.

## Notes

- Policy source: `C:\Users\haoyu\Desktop\code\github\LifeCycle\outputs\grid-test-20x20-results\scenario_0002\year06.txt`
- Risk metrics are Monte Carlo estimates using the model's return assumptions, a 4% withdrawal assumption, and an age-based projection horizon.
- Failure probabilities include an approximate 95% Monte Carlo margin of error in percentage points.
- Assumed planning horizon: 30 years (age 25-55). Risk horizon is determined based on the client's lifecycle stage.
