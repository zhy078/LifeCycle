# Output Spec

Expected output directory files:

- `results.jsonl`: one line per evaluated scenario.
- `best_params.json`: winning parameter set and score.
- `report.md`: human-readable optimization summary.
- `customer_manager_report.md`: advisor-facing report with client-ready talking points.

Each scenario object should include:

- `params`
- `objective`
- `score`
- `status`
- `artifact_dir`
- `policy_summary_path` when real model artifacts exist
- `lifecycle_checkpoints` for concise age-based policy differences

Optional config input:

- `client_profile.current_wealth` to choose the nearest wealth path when generating advisor-facing recommendations.
