# LifeCycle 多场景生命周期效用汇总

## 重要更正

- LifeCycle 模型比较的应是生命周期 utility/value function，不是 terminal wealth。
- 本文件用客户当前年龄 45 岁、当前财富 10.25，在每个场景的 year26.txt 中对 value function 插值得到 lifecycle utility。
- 下面的 `old_terminal_metric` 是旧口径留下的对照，不作为正式排序依据。

## 运行结论

- 目标函数: maximize_lifetime_utility
- 最优风险偏好参数 rho: 6.0
- 最优生命周期效用: 2.01673565
- 模型状态: ok
- 已完成场景数: 3

## 场景排名

| 排名 | scenario | status | lifecycle_utility | rho（风险厌恶/偏好） | delta | psi | mu | sigr | old_terminal_metric |
|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0 | ok | 2.01673565 | 6.0 | 0.97 | 0.6 | 0.03 | 0.15 | 0.05542695 |
| 2 | 1 | ok | 1.81891122 | 8.0 | 0.97 | 0.6 | 0.03 | 0.15 | 0.05541207 |
| 3 | 2 | ok | 1.64680007 | 10.0 | 0.97 | 0.6 | 0.03 | 0.15 | 0.05528873 |

## 解释

- `rho` 是 agent 的风险厌恶/风险偏好参数；rho 越高，风险厌恶越强。
- 在这三个候选点里，rho=6 的生命周期效用最高。
- 这个结论来自完整真实模型已生成的 value function，不是根据终端财富排序。

## 文件

- 本汇总: C:\Users\haoyu\Desktop\code\github\LifeCycle\outputs\lifecycle-optimizer-demo-multi-results\multi_results_lifetime_utility_summary.md
- 原始旧口径结果: C:\Users\haoyu\Desktop\code\github\LifeCycle\outputs\lifecycle-optimizer-demo-multi-results\results.jsonl
- 客户经理报告: C:\Users\haoyu\Desktop\code\github\LifeCycle\outputs\lifecycle-optimizer-demo-multi-results\customer_manager_report.md
