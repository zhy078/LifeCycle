# LifeCycle Optimization Report

{
  "objective": "maximize_terminal_wealth",
  "search_method": "grid",
  "use_real_model": true,
  "fast_mode": true,
  "total_scenarios": 2,
  "failed_scenarios": 2,
  "total_model_runtime_seconds": 555.088,
  "elapsed_seconds": 555.11,
  "best": {
    "id": 0,
    "params": {
      "rho": 6.0,
      "delta": 0.97,
      "psi": 0.6,
      "mu": 0.03,
      "sigr": 0.15
    },
    "objective": "maximize_terminal_wealth",
    "score": -Infinity,
    "status": "error",
    "artifact_dir": "C:\\Users\\haoyu\\Desktop\\code\\github\\LifeCycle\\outputs\\lifecycle-integration-test-2\\scenario_0000",
    "metric": null,
    "run_seconds": 280.313,
    "year_files_count": 0,
    "octave_command": "octave --quiet life_cycle.m"
  },
  "top_k": [
    {
      "id": 0,
      "params": {
        "rho": 6.0,
        "delta": 0.97,
        "psi": 0.6,
        "mu": 0.03,
        "sigr": 0.15
      },
      "objective": "maximize_terminal_wealth",
      "score": -Infinity,
      "status": "error",
      "artifact_dir": "C:\\Users\\haoyu\\Desktop\\code\\github\\LifeCycle\\outputs\\lifecycle-integration-test-2\\scenario_0000",
      "metric": null,
      "run_seconds": 280.313,
      "year_files_count": 0,
      "octave_command": "octave --quiet life_cycle.m"
    },
    {
      "id": 1,
      "params": {
        "rho": 8.0,
        "delta": 0.97,
        "psi": 0.4,
        "mu": 0.04,
        "sigr": 0.25
      },
      "objective": "maximize_terminal_wealth",
      "score": -Infinity,
      "status": "error",
      "artifact_dir": "C:\\Users\\haoyu\\Desktop\\code\\github\\LifeCycle\\outputs\\lifecycle-integration-test-2\\scenario_0001",
      "metric": null,
      "run_seconds": 274.775,
      "year_files_count": 0,
      "octave_command": "octave --quiet life_cycle.m"
    }
  ]
}