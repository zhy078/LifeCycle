# Output Spec

Expected output directory files:

- `results.jsonl`: one line per evaluated scenario.
- `best_params.json`: winning parameter set and score.
- `report.md`: human-readable optimization summary.

Each scenario object should include:

- `params`
- `objective`
- `score`
- `status`
- `artifact_dir`
