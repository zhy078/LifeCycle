# Weird Cases

筛选条件：

- `wealth <= 5 AND alpha > 0.5`
- `age >= 60 AND alpha > 0.5`

## 异常点表

| case | 问题 | 可能原因 |
|---|---|---|
| rho=6, age=60, wealth=20, alpha=0.8936, utility=2.6068 | 接近退休/退休前后但高风险配置 | fast-mode 网格太粗（na=10, ncash=10）造成 alpha 跳变；退休附近 value function 曲率或插值造成局部非单调；高财富客户有更多风险缓冲，但该解释需要复核；需要用更细 wealth grid 验证 |
| rho=8, age=60, wealth=20, alpha=0.8936, utility=2.4393 | 接近退休/退休前后但高风险配置 | fast-mode 网格太粗（na=10, ncash=10）造成 alpha 跳变；退休附近 value function 曲率或插值造成局部非单调；高财富客户有更多风险缓冲，但该解释需要复核；需要用更细 wealth grid 验证 |
| rho=10, age=60, wealth=20, alpha=0.8936, utility=2.3070 | 接近退休/退休前后但高风险配置 | fast-mode 网格太粗（na=10, ncash=10）造成 alpha 跳变；退休附近 value function 曲率或插值造成局部非单调；高财富客户有更多风险缓冲，但该解释需要复核；需要用更细 wealth grid 验证 |

## 解释

本轮没有发现 `wealth <= 5 AND alpha > 0.5` 的低财富高风险异常。

异常主要来自 `age >= 60 AND alpha > 0.5`，尤其 wealth=20 时 alpha 接近 0.89。这个结果可能有两种解释：

- 合理解释：高财富客户有更强风险缓冲，模型允许承担更多权益风险来提升生命周期效用。
- 不合理解释：fast-mode 的网格太粗，导致 alpha 在 wealth=20 附近跳到极端值。

需要进一步检查：

- 增大 `na` 和 `ncash`，看 alpha 是否变平滑。
- 加密 wealth grid，尤其 wealth=10 到 20 附近。
- 对退休前后 age=55-70 做更细年龄切片。
- 做 stochastic validation，检查 naive/optimal 的收益、回撤和失败概率分布，而不是只看均值。
