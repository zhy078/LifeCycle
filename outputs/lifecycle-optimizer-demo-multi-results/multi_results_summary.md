# LifeCycle 多场景优化汇总

## 运行结论

- 目标函数: maximize_terminal_wealth
- 最优参数: rho=6.0, delta=0.97, psi=0.6, mu=0.03, sigr=0.15
- 最优场景分数: 0.05542695
- 模型状态: ok
- 已完成场景数: 3
- 总模型运行秒数: 521.222

## 场景排名

| 排名 | scenario | status | score | rho | delta | psi | mu | sigr | run_seconds | year_files |
|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0 | ok | 0.05542695 | 6.0 | 0.97 | 0.6 | 0.03 | 0.15 | 173.392 | 80 |
| 2 | 1 | ok | 0.05541207 | 8.0 | 0.97 | 0.6 | 0.03 | 0.15 | 173.723 | 80 |
| 3 | 2 | ok | 0.05528873 | 10.0 | 0.97 | 0.6 | 0.03 | 0.15 | 174.107 | 80 |

## 最优场景客户建议摘录

- 当前年龄: 45岁
- 当前财富: current_wealth=10.25
- 当前建议组合: 稳健型组合
- 当前建议配比: 权益类约11%，稳健类约89%
- 当前建议消费强度参考: 4.3813

## 最优场景生命周期检查点

| 年龄 | 阶段 | 中财富权益比例 alpha | 中财富消费参考 |
|---:|---|---:|---:|
| 20 | 工作期 | 33.33% | 4.1460 |
| 64 | 工作期 | 44.44% | 4.0772 |
| 65 | 退休期 | 44.44% | 3.8639 |
| 99 | 退休期 | 22.22% | 5.5329 |

## 输出文件

- 汇总文件: C:\Users\haoyu\Desktop\code\github\LifeCycle\outputs\lifecycle-optimizer-demo-multi-results\multi_results_summary.md
- 客户经理报告: C:\Users\haoyu\Desktop\code\github\LifeCycle\outputs\lifecycle-optimizer-demo-multi-results\customer_manager_report.md
- 原始结果: C:\Users\haoyu\Desktop\code\github\LifeCycle\outputs\lifecycle-optimizer-demo-multi-results\results.jsonl
- 最优参数: C:\Users\haoyu\Desktop\code\github\LifeCycle\outputs\lifecycle-optimizer-demo-multi-results\best_params.json

## 备注

- 本次只测试 rho = 6、8、10 三个候选点，其它参数固定。
- 三个场景均为真实 Octave 模型完整跑完，status 均为 ok。
