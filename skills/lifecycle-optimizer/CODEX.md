# Codex Context: lifecycle-optimizer

This file is persistent context for Codex/subagents working on the `lifecycle-optimizer` skill.
Read it before changing scripts, reports, or examples in this skill.

## Repository Layout

The project is now organized as a reusable system:

- `model/`: core MATLAB/Octave model source (`life_cycle.m` and helper `f_*.m` files).
- `data/`: saved root-level model output snapshots (`yearXX.txt`, `CWY.txt`, `CWYs.txt`, `SB.txt`).
- `outputs/`: run-specific experiment outputs, scenario artifacts, reports, CSVs, and figures.
- `skills/lifecycle-optimizer/`: Codex skill layer, scripts, tests, config assets, reports, and agent prompts.
- `app/`: project-level CLI / agent entrypoints.

The real-model runner copies source files from `model/` into each scenario artifact directory before running Octave.
Do not put new runtime artifacts in the repository root.

## Product Layer

The current product shell has five required modules:

- Client input layer: `app/client_input.py`
  - parses text such as `45岁，100万，偏稳健`
  - maps `100万` to model wealth `10.0` using the demo convention `1 model unit = 100k`
- Risk preference layer:
  - `conservative -> rho=10`
  - `balanced -> rho=8`
  - `aggressive -> rho=6`
- Decision explanation layer: `skills/lifecycle-optimizer/scripts/policy_explainer.py`
- Risk-control layer: `app/risk_metrics.py`
  - max drawdown
  - failure probability
  - annual volatility
- Demo entrypoints:
  - CLI: `py app/run_agent.py "45岁，100万，偏稳健"`
  - Streamlit: `streamlit run app/app.py`
- Batch validation:
  - `py app/batch_compare.py --clients 100 --paths 1000`
  - compares interpolated lifecycle policy against naive fixed 60/40 using shared Monte Carlo shocks
  - writes `outputs/client-batch-comparison/client_100_comparison.csv`
  - writes `outputs/client-batch-comparison/client_100_comparison_summary.md`

The CLI advisor defaults to the completed `20x20` artifacts under `outputs/grid-test-20x20-results`.

## Role Setup

When testing this skill, use two roles:

- Customer manager: uses the skill with client-provided assumptions and explains results to the client.
- Backend engineer: maintains the skill implementation, tests, docs, and runner behavior.

The customer manager must not present incomplete model output as advice. Timeout output is diagnostic only.

## Core Model Meaning

This LifeCycle model is not primarily a terminal-wealth optimizer.

Correct business interpretation:

- The model objective is lifecycle utility / value function.
- Use `maximize_lifetime_utility` as the default real-model objective.
- `rho` is the agent risk aversion / risk preference parameter.
- Higher `rho` means stronger risk aversion and more conservative preference.
- Lower `rho` means greater willingness to accept risk.

Do not describe `rho` as a generic parameter without explaining risk preference.

## Real-Model Scoring

For real Octave model runs, rank scenarios by lifecycle utility:

1. Determine the client age from `client_profile.current_age`; fall back to `fixed.tb`.
2. Map age to the policy file:
   - `file_idx = current_age - tb + 1`
   - example: `tb=20`, `current_age=45` -> `year26.txt`
3. Parse `yearXX.txt` as three equal blocks:
   - first block: portfolio share `alpha`
   - second block: consumption
   - third block: value function / utility
4. Interpolate the third block at `client_profile.current_wealth`.
5. Use that interpolated value as `score` / lifecycle utility.

The older terminal metric from `CWY.txt` can be kept only as a diagnostic comparison. It must not drive official ranking.

## Report Language

Customer-facing or customer-manager-facing reports should say:

- "生命周期效用" or "lifecycle utility"
- "value function"
- "rho（风险厌恶/风险偏好）"
- "模型状态: ok" only when the scenario fully completed

Avoid saying:

- "maximize_terminal_wealth" as the business objective
- "terminal wealth" as the official optimization target
- any recommendation when `status != ok`, except clearly marked as diagnostic/preview

## Runner Requirements

On Windows, prefer `octave-cli` over `octave` to avoid launching `octave-gui`.

The real-model runner must:

- write Octave stdout/stderr to files, not a captured pipe;
- kill the process tree on timeout;
- write `runner_diagnostics.json`;
- flush each `results.jsonl` row immediately;
- clear stale scenario artifacts before rerunning a scenario;
- keep `score=null` for incomplete strict real-model scenarios;
- write valid JSON with no `Infinity` or `-Infinity`.

## Report Artifacts

Each completed run should provide:

- `results.jsonl`
- `best_params.json`
- `report.md`
- `customer_manager_report.md`
- scenario-level `runner_diagnostics.json`
- scenario-level `policy_summary.json` when year files exist

For multi-scenario comparisons, create a single readable summary file, for example:

- `multi_results_lifetime_utility_summary.md`

That file should include:

- objective definition;
- ranking by lifecycle utility;
- `rho` explained as risk aversion / preference;
- full parameter vector;
- status per scenario;
- a note if any old terminal metric is shown only for comparison.

## Known Example

Example client profile used during testing:

```json
{
  "current_age": 45,
  "current_wealth": 10.25
}
```

Example three-scenario comparison with `rho = [6, 8, 10]` ranked by lifecycle utility:

| rank | rho | lifecycle_utility |
|---:|---:|---:|
| 1 | 6.0 | 2.01673565 |
| 2 | 8.0 | 1.81891122 |
| 3 | 10.0 | 1.64680007 |

These values came from complete real-model artifacts by reading `year26.txt` value functions for a 45-year-old client with wealth `10.25`.

## 27-Scenario Analysis Pattern

For client-state sensitivity, expand completed model artifacts across:

- `rho = [6, 8, 10]`
- `wealth = [5, 10, 20]`
- `age = [30, 45, 60]`

This is 27 scored client states, but only 3 unique Octave model solves if only `rho` changes. Age and wealth are scoring/query dimensions against the policy/value files.

Expected summary artifacts:

- `scenario_27_lifetime_utility_with_naive_summary.md`
- `scenario_27_lifetime_utility_with_naive.csv`

The summary should include:

- whether the model matches economic intuition;
- which variables matter most;
- unreasonable/suspicious model behavior;
- a baseline comparison table:
  - `model optimal policy`
  - `naive 60/40 + 4% consumption`

Important: the naive baseline utility currently uses a deterministic CRRA comparison simulator. It is useful as a baseline direction check, but it is not the exact model value-function utility under a fixed-policy constraint.

## Grid Resolution Finding

Latest grid validation found that `20x20` is enough for current diagnostics.

Problem state tested:

- `age = 60`
- `wealth = 20`
- `rho = [6, 8, 10]`

Key result:

| grid | alpha behavior | runtime implication |
|---|---|---|
| `10x10` | produced extreme alpha around `0.8936` | fast but too coarse |
| `20x20` | alpha dropped to a moderate range around `0.26-0.45` | good default diagnostic grid |
| `40x40` | very close to `20x20` but much slower | use only for final robustness checks |

Conclusion:

- Treat the original high-alpha anomaly as mostly a coarse-grid artifact.
- Use `na=20`, `ncash=20` as the default nontrivial validation grid.
- Do not default to `40x40`; it is expensive and did not materially change the conclusion relative to `20x20`.
- Keep `10x10` only for smoke tests and quick UI/report checks.

## Useful Commands

Run unit tests:

```powershell
py -m unittest discover -s skills/lifecycle-optimizer/scripts -p test_*.py -v
```

Project-level dry run:

```powershell
py app/cli.py --dry-run --max-evals 2
```

Fast dry run:

```powershell
py skills/lifecycle-optimizer/scripts/optimize.py --config skills/lifecycle-optimizer/assets/sample-case.json --output-dir outputs/lifecycle-optimizer-dry --dry-run --max-evals 2 --progress-every 1
```

Real-model smoke with explicit timeout:

```powershell
py app/cli.py --real --fast-mode --max-evals 1 --timeout-sec 1200
```

## Current Artifact Locations

Root model outputs were moved out of the repository root:

- model source is in `model/`
- saved policy/value snapshots are in `data/`
- `outputs/`

Do not revert generated outputs unless the user explicitly asks. Treat them as user or generated state.
