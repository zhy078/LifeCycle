# Model Insights

## 1. 模型是否合理

结论：模型基本符合生命周期经济学直觉，但存在局部异常。

支持合理性的证据：

- wealth 上升时，平均 lifecycle utility 上升。
- rho 越高，平均 lifecycle utility 越低，符合更高风险厌恶带来更保守偏好的解释。
- optimal policy 在 deterministic CRRA baseline 对比中优于 naive 60/40 + 4% 消费规则。

局部异常：

- `age >= 60` 且 `wealth=20` 的场景出现 `alpha > 0.5`，最高约 0.8936。
- `age=45, rho=8` 时 alpha 随 wealth 上升而下降，不符合简单财富风险承受力直觉。
- `wealth=10, rho=8` 时 alpha 在 age=60 反弹，不符合简单年龄风险下降直觉。

## 2. 哪些变量最重要

### wealth（最重要）

| wealth | avg lifecycle utility |
|---:|---:|
| 5 | 1.982781 |
| 10 | 2.052583 |
| 20 | 2.129838 |

wealth 直接决定 value function 上的位置，是 utility 水平最明显的驱动因素。

### rho（偏好）

| rho | avg lifecycle utility |
|---:|---:|
| 6 | 2.237871 |
| 8 | 2.045215 |
| 10 | 1.882117 |

rho 是风险厌恶/风险偏好参数。rho 越高，模型效用水平越低，并且很多状态下 alpha 更保守。

### age（阶段）

| age | avg lifecycle utility |
|---:|---:|
| 30 | 2.045590 |
| 45 | 1.818319 |
| 60 | 2.301293 |

age 决定客户处于积累期、退休前还是退休期，也决定读取哪个 `yearXX.txt` policy/value 文件。

## 3. 模型价值

相比 naive 60/40，模型在所有场景下提供更高生命周期效用。

| 方法 | utility |
|---|---:|
| model optimal policy | -0.00028212 |
| naive 60/40 + fixed consumption rate | -0.02131200 |

所有 27 个场景中，optimal policy 都优于 naive baseline。

注意：这里的 naive baseline 是 deterministic CRRA simulator 的比较口径，不是把 naive rule 放进 Octave 重新求解的 value function。

## 4. 局限性

- grid coarse：当前 fast-mode 使用较粗网格，例如 `na=10`, `ncash=10`，可能造成 alpha 跳变。
- extreme alpha：部分场景出现接近 0 或接近 1 的 alpha，需要用更细网格复核。
- no stochastic validation：当前 baseline 对比是 deterministic simulation，还没有做 Monte Carlo 分布、回撤、失败概率等随机验证。
- baseline not solved in model：naive 60/40 + 固定消费率没有作为固定策略在原模型中重新求解，只是外部对比。
- limited parameter sweep：本次只扫了 rho、age、wealth，还没有系统测试 return、volatility、discount factor、retirement age、psi。

## 5. 下一步

- 用更细网格重跑异常区域：age=55-70, wealth=10-20。
- 把 naive baseline 做成正式 fixed-policy evaluator。
- 输出 alpha heatmap 或表格，让客户经理能解释 policy shape。
- 对 `wealth <= 5` 和 `age >= 60` 的异常规则加入自动 warning。

## 相关文件

- weird cases: `C:\Users\haoyu\Desktop\code\github\LifeCycle\outputs\lifecycle-optimizer-demo-multi-results\weird_cases.md`
- policy shape: `C:\Users\haoyu\Desktop\code\github\LifeCycle\outputs\lifecycle-optimizer-demo-multi-results\policy_shape_analysis.md`
- 27-scenario CSV: `C:\Users\haoyu\Desktop\code\github\LifeCycle\outputs\lifecycle-optimizer-demo-multi-results\scenario_27_lifetime_utility_with_naive.csv`
