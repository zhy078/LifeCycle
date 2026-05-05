# Visual Explanation

## 1. 粗网格是否导致异常 alpha？

![Grid refinement alpha](C:/Users/haoyu/Desktop/code/github/LifeCycle/outputs/grid-validation-analysis/figures/grid_refinement_alpha.png)

解释：10x10 网格下，`age=60, wealth=20` 的 alpha 接近 0.89；20x20 和 40x40 后下降到约 0.26-0.45。这个图强烈支持：原始极端风险配置主要是 coarse grid artifact。

![Grid refinement utility](C:/Users/haoyu/Desktop/code/github/LifeCycle/outputs/grid-validation-analysis/figures/grid_refinement_utility.png)

解释：utility 在 20x20 到 40x40 之间变化很小，说明更细网格后结果开始稳定。

## 2. age 反弹问题

![Age rebound alpha](C:/Users/haoyu/Desktop/code/github/LifeCycle/outputs/grid-validation-analysis/figures/age_rebound_alpha.png)

解释：10x10 下 alpha 先降后升；20x20/40x40 更平滑，但 age=60 附近仍有轻微反弹。这个点不是完全消失了，应该作为 retirement transition 附近的 policy-shape 诊断。

## 3. 异常点柱状图

![Problem state alpha bars](C:/Users/haoyu/Desktop/code/github/LifeCycle/outputs/grid-validation-analysis/figures/problem_state_alpha_bars.png)

解释：同一个异常状态下，10x10 明显高于 warning line，而 20x20/40x40 回到较合理区间。

## 4. Stochastic validation

![Stochastic utility](C:/Users/haoyu/Desktop/code/github/LifeCycle/outputs/grid-validation-analysis/figures/stochastic_utility_methods.png)

解释：1000 paths Monte Carlo 下，model optimal policy 的平均 utility 高于 naive 60/40 + 4% consumption。

![Stochastic terminal wealth](C:/Users/haoyu/Desktop/code/github/LifeCycle/outputs/grid-validation-analysis/figures/stochastic_terminal_wealth_methods.png)

解释：naive baseline 的 terminal wealth 反而更高，但 utility 更低。这正好说明 LifeCycle 模型不是最大化终端财富，而是在消费、风险和生命周期效用之间权衡。

## Files

- Figures directory: `C:/Users/haoyu/Desktop/code/github/LifeCycle/outputs/grid-validation-analysis/figures`
- Source report: `C:/Users/haoyu/Desktop/code/github/LifeCycle/outputs/grid-validation-analysis/grid_and_stochastic_validation.md`
