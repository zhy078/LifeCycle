---
name: lifecycle-optimizer
description: Optimize life-cycle consumption/portfolio model settings from user-provided assumptions and constraints. Use when users provide preconditions (risk preference, retirement age, market return/volatility, simulation size, objective and bounds) and ask for the best parameter set, best run summary, or scenario comparison for life-cycle model outputs.
---

# lifecycle-optimizer

Collect assumptions, generate case files, run model evaluations, rank outcomes by lifecycle utility, and return a concise optimization report.

Repository layout:

- Core Octave model source lives in `model/`.
- Saved policy/value snapshots live in `data/`.
- Run-specific artifacts live in `outputs/`.
- Project-level entrypoints live in `app/`.

## Workflow

1. Parse user inputs into a case JSON using `assets/input-template.json`.
2. Validate ranges and required fields (see `references/model-params.md`).
3. Start with a fast real-model smoke run.
4. Execute optimization search (`grid` or `random`) through `scripts/optimize.py`.
5. Return:
   - best parameter set,
   - objective score,
   - top-K alternatives,
   - output directory and artifacts.

## Run Commands

### Fast dry-run check

```bash
py skills/lifecycle-optimizer/scripts/optimize.py \
  --config skills/lifecycle-optimizer/assets/sample-case.json \
  --output-dir outputs/lifecycle-optimizer \
  --dry-run
```

### Fast real-model smoke run (recommended first)

```bash
py skills/lifecycle-optimizer/scripts/optimize.py \
  --config skills/lifecycle-optimizer/assets/sample-case.json \
  --output-dir outputs/lifecycle-optimizer-real \
  --use-real-model --fast-mode --max-evals 2 --timeout-sec 180 --progress-every 1
```

## Reporting Format

Always provide:

- objective definition and value. For the real model, the default objective is `maximize_lifetime_utility`: score is the value-function utility at the client's current age and wealth path, not terminal wealth,
- chosen parameter vector,
- validation notes,
- paths to `best_params.json`, `results.jsonl`, and `report.md`.


For real-model verification, check `results.jsonl` fields per scenario:
- `status` must be `ok`
- `year_files_count` should be > 0
- `run_seconds` should be non-trivial
- `octave_command` should be present

## Notes

- Keep objective explicitly declared; never assume.
- Keep reproducibility by setting `seed`.
- If model execution fails, return actionable diagnostics and suggest dry-run verification.

## Common CLI Pitfall

Do not type literal `\n` in shell commands. Use real line breaks with `\` continuation, or run a single-line command.

## Output Path Behavior

Relative `--output-dir` values are anchored to the repository root (`LifeCycle/`), not your current terminal folder.

Example on Windows clone path:

```bash
py skills/lifecycle-optimizer/scripts/optimize.py --config skills/lifecycle-optimizer/assets/sample-case.json --output-dir outputs/lifecycle-optimizer --dry-run
```

This writes to: `C:\Users\haoyu\Desktop\code\github\LifeCycle\outputs\lifecycle-optimizer`.

Project-level equivalent:

```bash
py app/cli.py --dry-run --max-evals 2
```


## Real-Model Integrity

When `--use-real-model` is enabled, the runner now requires Octave by default.
If Octave is missing, it exits with an error instead of silently returning proxy scores.
The real-model runner now executes inside each scenario artifact directory so generated `year*.txt`, `CWY.txt`, `CWYs.txt`, and `SB.txt` stay attached to the scenario that produced them.
The runner copies `life_cycle.m` and helper `f_*.m` files from `model/` into each scenario artifact directory before patching parameters.
Fast mode keeps the model's Gaussian quadrature size intact and only reduces lightweight dimensions such as `nsim`, so it remains compatible with the current `life_cycle.m`.

Use `--allow-proxy-fallback` only when you explicitly want fallback behavior for debugging.


## Integration Test (Skill + Octave)

Use this to validate end-to-end integration:

```bash
py skills/lifecycle-optimizer/scripts/test_real_integration.py --max-evals 2 --timeout-sec 180
```

For strict real-model validation (require all scenarios `status=ok` and year files present):

```bash
py skills/lifecycle-optimizer/scripts/test_real_integration.py --max-evals 2 --timeout-sec 180 --strict-real
```

Timeout diagnostics are written per scenario to `runner_diagnostics.json`.
When a strict real-model scenario times out, `results.jsonl` still receives a row with `status=error_timeout`, `timed_out=true`, and `score=null`.
