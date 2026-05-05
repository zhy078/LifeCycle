# LifeCycle Optimization Report

{
  "objective": "maximize_terminal_wealth",
  "search_method": "grid",
  "use_real_model": true,
  "fast_mode": true,
  "total_scenarios": 2,
  "failed_scenarios": 0,
  "total_model_runtime_seconds": 533.557,
  "elapsed_seconds": 533.575,
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
    "score": 200.0,
    "status": "ok",
    "artifact_dir": "C:\\Users\\haoyu\\Desktop\\code\\github\\LifeCycle\\outputs\\lifecycle-integration-test-6\\scenario_0000",
    "metric": 200.0,
    "run_seconds": 272.081,
    "year_files_count": 80,
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
      "score": 200.0,
      "status": "ok",
      "artifact_dir": "C:\\Users\\haoyu\\Desktop\\code\\github\\LifeCycle\\outputs\\lifecycle-integration-test-6\\scenario_0000",
      "metric": 200.0,
      "run_seconds": 272.081,
      "year_files_count": 80,
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
      "score": 200.0,
      "status": "ok",
      "artifact_dir": "C:\\Users\\haoyu\\Desktop\\code\\github\\LifeCycle\\outputs\\lifecycle-integration-test-6\\scenario_0001",
      "metric": 200.0,
      "run_seconds": 261.476,
      "year_files_count": 80,
      "octave_command": "octave --quiet life_cycle.m"
    }
  ]
}