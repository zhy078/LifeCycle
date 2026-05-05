# Policy Shape Analysis

## 1. Alpha vs Wealth

固定：`age=45, rho=8`

| wealth | alpha | interpretation |
|---:|---:|---|
| 5 | 0.3257 | baseline |
| 10 | 0.0156 | down |
| 20 | 0.0000 | down |

结论：

- 这个切片中 alpha 随 wealth 上升反而下降：5 -> 10 -> 20 时约 0.326 -> 0.016 -> 0。
- 这不符合简单的“财富越高越能承担风险”直觉。
- 可能原因：value function 下的消费/储蓄权衡、现金状态点离散、粗网格，以及当前 wealth 区间的插值。
- 这是论文/面试里可以强调的亮点：模型结构整体可用，但局部 policy shape 需要诊断。

## 2. Alpha vs Age

固定：`wealth=10, rho=8`

| age | alpha | interpretation |
|---:|---:|---|
| 30 | 0.1215 | baseline |
| 45 | 0.0156 | down |
| 60 | 0.2118 | up |

结论：

- alpha 从 age=30 到 45 明显下降，但到 age=60 又上升。
- 简单生命周期直觉通常要求接近退休时风险下降，因此 age=60 的反弹需要解释。
- 可能原因包括：退休前劳动收入/退休收入切换、网格过粗、value function 非线性，以及模型中消费与投资同时优化导致的局部替代效应。

## Product Explanation Example

```text
Recommended portfolio: 89% risky assets, 11% stabilizing assets.
Lifecycle utility: 2.6068.

Explanation:
由于接近退休，模型通常应降低风险暴露，优先保护退休前后的现金流稳定性；财富较高时，客户有更多风险缓冲，但仍需要结合年龄阶段控制回撤；rho 较低表示风险承受意愿更高，模型更可能接受权益波动。

Caution:
接近退休却高风险配置，可能违背常识，需要检查 fast-mode 粗网格、边界条件或 value function 插值。
```
