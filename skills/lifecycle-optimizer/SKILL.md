---
name: lifecycle-optimizer
description: Optimize life-cycle consumption/portfolio model settings from user-provided assumptions and constraints. Use when users provide preconditions (risk preference, retirement age, market return/volatility, simulation size, objective and bounds) and ask for the best parameter set, best run summary, or scenario comparison for life-cycle model outputs.
---

# lifecycle-optimizer

Collect assumptions, generate case files, run model evaluations, rank outcomes by objective, and return a concise optimization report.

## Workflow

1. Parse user inputs into a case JSON using `assets/input-template.json`.
2. Validate ranges and required fields (see `references/model-params.md`).
3. Execute one baseline run to ensure runtime wiring works.
4. Execute optimization search (`grid` or `random`) through `scripts/optimize.py`.
5. Return:
   - best parameter set,
   - objective score,
   - top-K alternatives,
   - output directory and artifacts.

## Run Commands

Use:

```bash
python3 skills/lifecycle-optimizer/scripts/optimize.py \
  --config skills/lifecycle-optimizer/assets/sample-case.json \
  --output-dir outputs/lifecycle-optimizer
```

If GNU Octave is unavailable, run dry-run mode:

```bash
python3 skills/lifecycle-optimizer/scripts/optimize.py \
  --config skills/lifecycle-optimizer/assets/sample-case.json \
  --output-dir outputs/lifecycle-optimizer \
  --dry-run
```

## Reporting Format

Always provide:

- objective definition and value,
- chosen parameter vector,
- validation notes,
- paths to `best_params.json`, `results.jsonl`, and `report.md`.

## Notes

- Keep objective explicitly declared; never assume.
- Keep reproducibility by setting `seed`.
- If model execution fails, return actionable diagnostics and suggest dry-run verification.


## Common CLI Pitfall

Do not type literal `\n` in shell commands. Use real line breaks with `\` continuation, or run a single-line command.
