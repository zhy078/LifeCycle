# Lifecycle Client Advice

## Parsed Client Input

- Age: 45
- Wealth: 1.00 model units
- Risk preference: aggressive -> rho=6

## Model Output

- Model status: ok
- Matched rho: 6.0
- Raw model risky allocation alpha: 0.8178
- Recommended consumption: 0.6936
- Lifecycle utility: 1.3309

## 生命周期配置路径

| 年龄区间 | 风险资产比例 | 策略目标 |
|---|---:|---|
| 45-50 | 2% | 保本为主，低波动，保留少量增长敞口 |
| 50-60 | 0%-2% | 逐步降低风险，减少权益波动对退休准备金的影响 |
| 60-75 | 0%-3% | 现金流稳定，优先控制本金回撤和提前耗尽风险 |

这不是一次性的单点建议，而是一条随年龄自动下调风险暴露的 glide path。随着年龄接近退休阶段，模型会提高对本金安全和现金流稳定的权重；即使客户风险偏好较高，也会逐步降低风险资产配置。

## 年度消费与取现建议

- 当前建议支取: 约 3%-4%/年 (0.03-0.04 model units)
- 模型 consumption 参考值: 0.6936，在客户经理场景中转化为可执行的年度支取区间，而不是要求客户按模型单位直接消费。
- 建议避免过度提取，防止资产在退休前后提前耗尽。
- 在当前配置和取现假设下，资产可维持至 75 岁的模拟置信下界约为 100.0%，满足 95% 风控口径。

## 未来资产路径

| 情景 | 目标年龄 | 预计资产 |
|---|---:|---:|
| 保守情景 P10 | 75 | 0.43 |
| 中位情景 P50 | 75 | 0.44 |
| 乐观情景 P90 | 75 | 0.46 |

这张表按 3% 年度支取口径展示，把 failure probability 翻译成客户语言：在当前生命周期配置和支取纪律下，资产大概率不会在规划期内耗尽；4% 支取作为上方风控压力测试。

## 分阶段产品映射

| 阶段 | 产品方向 | 执行重点 |
|---|---|---|
| 当前阶段 | 低波动固收产品 + 少量宽基指数基金 | 控制权益仓位，优先降低回撤和波动 |
| 50岁以后 | 提高现金类和短久期固收比例 | 逐步减少权益类暴露，强化退休现金流准备 |
| 60岁以后 | 现金管理、存款、稳健固收为主 | 稳定支取来源，避免大幅净值波动影响生活安排 |

## 自动调整规则

- 年龄增长: 按生命周期路径年度自动下调风险资产比例。
- 市场下跌超过 20%: 暂停加风险，优先检查现金流和本金安全。
- 财富水平显著上升: 若资产翻倍且现金流压力下降，可将风险资产比例小幅提高，但上限建议不超过 10%。
- 财富水平下降或支取压力上升: 降低风险资产比例，优先保留 12-24 个月流动性。

## Risk Control Decision

- Case type: UNSAFE_MODEL_OUTPUT
- Decision status: adjusted
- Suggested allocation range: 0%-7% risky assets
- Point estimate alpha: 0.0178
- Failure probability limit: 5.0%
- Risk horizon: 30 years (age 45-75; pre-retirement planning window)

## Decision Trace

| Step | Before | After | Reason |
|---|---:|---:|---|
| failure_risk | 82% | 2% | Reduced until estimated failure probability remains within 5.0% under Monte Carlo uncertainty bounds. |
| wealth_constraint | 2% | 2% | Applied low-wealth protection cap. |

The final point estimate remains at 2% because the failure-risk search is already the binding risk control.

## Failure Risk Search

| Candidate alpha | Failure estimate | 95% MC margin | Upper bound | Result |
|---:|---:|---:|---:|---|
| 82% | 37.4% | +/- 3.0 pp | 40.4% | Reduce |
| 77% | 36.4% | +/- 3.0 pp | 39.4% | Reduce |
| 72% | 34.5% | +/- 2.9 pp | 37.4% | Reduce |
| 67% | 35.8% | +/- 3.0 pp | 38.8% | Reduce |
| 62% | 34.8% | +/- 3.0 pp | 37.8% | Reduce |
| 57% | 35.2% | +/- 3.0 pp | 38.2% | Reduce |
| 52% | 32.8% | +/- 2.9 pp | 35.7% | Reduce |
| 47% | 32.4% | +/- 2.9 pp | 35.3% | Reduce |
| 42% | 29.0% | +/- 2.8 pp | 31.8% | Reduce |
| 37% | 28.7% | +/- 2.8 pp | 31.5% | Reduce |
| 32% | 24.4% | +/- 2.7 pp | 27.1% | Reduce |
| 27% | 22.3% | +/- 2.6 pp | 24.9% | Reduce |
| 22% | 20.3% | +/- 2.5 pp | 22.8% | Reduce |
| 17% | 18.0% | +/- 2.4 pp | 20.4% | Reduce |
| 12% | 10.7% | +/- 1.9 pp | 12.6% | Reduce |
| 7% | 4.1% | +/- 1.2 pp | 5.3% | Reduce |
| 2% | 0.0% | +/- 0.0 pp | 0.0% | Pass |

## Risk Metrics

| Metric | Raw model | Final recommendation |
|---|---:|---:|
| Max drawdown | -36.7% | 0.0% |
| Failure probability | 37.4% (95% MC margin +/- 3.0 pp) | 0.0% (95% MC margin +/- 0.0 pp) |
| Annual volatility | 12.2% | 0.3% |

## Baseline Comparison

| Metric | Standard 60/40 | Final recommendation |
|---|---:|---:|
| Risky allocation | 60% | 2% |
| Max drawdown | -26.5% | 0.0% |
| Failure probability | 34.7% (95% MC margin +/- 3.0 pp) | 0.0% (95% MC margin +/- 0.0 pp) |
| Annual volatility | 9.0% | 0.3% |
| Expected annual return | 2.4% | 1.5% |

## Applicability

This recommendation framework is suitable for pre-retirement clients who need to balance continued growth with rising capital-preservation needs.

## 行为约束与再平衡

生命周期模型假设客户能坚持长期计划，但真实客户容易被短期波动影响。若市场短期下跌导致账面亏损，建议避免频繁追涨杀跌，以年度周期进行再平衡；除非触发年龄、市场跌幅或财富水平变化等规则，不建议临时大幅调整配置。

## Customer Explanation

Model output adjusted due to risk controls.

The original model suggests a high-risk allocation despite the client's aggressive risk preference. This recommendation is adjusted to ensure:

- Failure probability remains below the product risk limit
- Risk exposure is appropriate for the client's lifecycle stage
- Wealth level is protected from large drawdowns
- Withdrawal needs can be supported through the planning horizon

Adjustment reasons:
- Failure probability 37.4% exceeds the 5.0% limit.
- Low wealth requires extra protection from large drawdowns.

Suggested allocation range: 0%-7% risky assets.
Point estimate: 2% risky assets, 98% stabilizing assets. This balances growth and capital preservation.

As the client moves closer to retirement, the lifecycle constraint automatically increases the weight on principal safety. Even with an aggressive risk preference, risky allocation should decline along the glide path.

For a client at this stage, preserving accumulated savings while keeping measured growth exposure is usually more important than maximizing risk-taking.

Compared with a standard 60/40 portfolio, the adjusted recommendation lowers drawdown, failure probability, and volatility, with a modestly lower expected annual return.

## Notes

- Policy source: `outputs\client-real-smoke\scenario_0000\year26.txt`
- Risk metrics are Monte Carlo estimates using the model's return assumptions, a 4% withdrawal assumption, and an age-based projection horizon.
- Failure probabilities include an approximate 95% Monte Carlo margin of error in percentage points.
- Assumed planning horizon: 30 years (age 45-75). Risk horizon is determined based on the client's lifecycle stage.
