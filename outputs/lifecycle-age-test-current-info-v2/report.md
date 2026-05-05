# LifeCycle Optimization Report

{
  "objective": "maximize_terminal_wealth",
  "search_method": "grid",
  "use_real_model": true,
  "fast_mode": true,
  "total_scenarios": 1,
  "failed_scenarios": 0,
  "total_model_runtime_seconds": 183.838,
  "elapsed_seconds": 183.858,
  "best": {
    "id": 0,
    "params": {
      "rho": 8.0,
      "delta": 0.95,
      "psi": 0.6,
      "mu": 0.03,
      "sigr": 0.15
    },
    "objective": "maximize_terminal_wealth",
    "score": 0.05544901,
    "status": "ok",
    "artifact_dir": "C:\\Users\\haoyu\\Desktop\\code\\github\\LifeCycle\\outputs\\lifecycle-age-test-current-info-v2\\scenario_0000",
    "metric": 0.05544901,
    "run_seconds": 183.838,
    "year_files_count": 80,
    "octave_command": "octave --quiet life_cycle.m",
    "policy_summary_path": "C:\\Users\\haoyu\\Desktop\\code\\github\\LifeCycle\\outputs\\lifecycle-age-test-current-info-v2\\scenario_0000\\policy_summary.json",
    "lifecycle_checkpoints": {
      "age_span": {
        "start": 20,
        "end": 99
      },
      "checkpoints": [
        {
          "age": 20,
          "phase": "working",
          "mid_wealth_alpha": 0.44444444,
          "mid_wealth_consumption": 4.1460421,
          "wealth_bands": {
            "low_wealth": {
              "alpha": 1.0,
              "consumption": 0.24975,
              "cash": 0.25
            },
            "mid_wealth": {
              "alpha": 0.44444444,
              "consumption": 4.1460421,
              "cash": 10.25103648
            },
            "high_wealth": {
              "alpha": 0.11111111,
              "consumption": 15.50875911,
              "cash": 200.0
            }
          }
        },
        {
          "age": 64,
          "phase": "working",
          "mid_wealth_alpha": 0.44444444,
          "mid_wealth_consumption": 4.07718741,
          "wealth_bands": {
            "low_wealth": {
              "alpha": 1.0,
              "consumption": 0.24975,
              "cash": 0.25
            },
            "mid_wealth": {
              "alpha": 0.44444444,
              "consumption": 4.07718741,
              "cash": 10.25103648
            },
            "high_wealth": {
              "alpha": 0.11111111,
              "consumption": 14.56096812,
              "cash": 200.0
            }
          }
        },
        {
          "age": 65,
          "phase": "retired",
          "mid_wealth_alpha": 0.44444444,
          "mid_wealth_consumption": 3.86387255,
          "wealth_bands": {
            "low_wealth": {
              "alpha": 1.0,
              "consumption": 0.25,
              "cash": 0.25
            },
            "mid_wealth": {
              "alpha": 0.44444444,
              "consumption": 3.86387255,
              "cash": 10.25103648
            },
            "high_wealth": {
              "alpha": 0.11111111,
              "consumption": 9.10576291,
              "cash": 200.0
            }
          }
        },
        {
          "age": 99,
          "phase": "retired",
          "mid_wealth_alpha": 0.22222222,
          "mid_wealth_consumption": 5.53293779,
          "wealth_bands": {
            "low_wealth": {
              "alpha": 1.0,
              "consumption": 0.25,
              "cash": 0.25
            },
            "mid_wealth": {
              "alpha": 0.22222222,
              "consumption": 5.53293779,
              "cash": 10.25103648
            },
            "high_wealth": {
              "alpha": 0.22222222,
              "consumption": 100.71937362,
              "cash": 200.0
            }
          }
        }
      ]
    }
  },
  "top_k": [
    {
      "id": 0,
      "params": {
        "rho": 8.0,
        "delta": 0.95,
        "psi": 0.6,
        "mu": 0.03,
        "sigr": 0.15
      },
      "objective": "maximize_terminal_wealth",
      "score": 0.05544901,
      "status": "ok",
      "artifact_dir": "C:\\Users\\haoyu\\Desktop\\code\\github\\LifeCycle\\outputs\\lifecycle-age-test-current-info-v2\\scenario_0000",
      "metric": 0.05544901,
      "run_seconds": 183.838,
      "year_files_count": 80,
      "octave_command": "octave --quiet life_cycle.m",
      "policy_summary_path": "C:\\Users\\haoyu\\Desktop\\code\\github\\LifeCycle\\outputs\\lifecycle-age-test-current-info-v2\\scenario_0000\\policy_summary.json",
      "lifecycle_checkpoints": {
        "age_span": {
          "start": 20,
          "end": 99
        },
        "checkpoints": [
          {
            "age": 20,
            "phase": "working",
            "mid_wealth_alpha": 0.44444444,
            "mid_wealth_consumption": 4.1460421,
            "wealth_bands": {
              "low_wealth": {
                "alpha": 1.0,
                "consumption": 0.24975,
                "cash": 0.25
              },
              "mid_wealth": {
                "alpha": 0.44444444,
                "consumption": 4.1460421,
                "cash": 10.25103648
              },
              "high_wealth": {
                "alpha": 0.11111111,
                "consumption": 15.50875911,
                "cash": 200.0
              }
            }
          },
          {
            "age": 64,
            "phase": "working",
            "mid_wealth_alpha": 0.44444444,
            "mid_wealth_consumption": 4.07718741,
            "wealth_bands": {
              "low_wealth": {
                "alpha": 1.0,
                "consumption": 0.24975,
                "cash": 0.25
              },
              "mid_wealth": {
                "alpha": 0.44444444,
                "consumption": 4.07718741,
                "cash": 10.25103648
              },
              "high_wealth": {
                "alpha": 0.11111111,
                "consumption": 14.56096812,
                "cash": 200.0
              }
            }
          },
          {
            "age": 65,
            "phase": "retired",
            "mid_wealth_alpha": 0.44444444,
            "mid_wealth_consumption": 3.86387255,
            "wealth_bands": {
              "low_wealth": {
                "alpha": 1.0,
                "consumption": 0.25,
                "cash": 0.25
              },
              "mid_wealth": {
                "alpha": 0.44444444,
                "consumption": 3.86387255,
                "cash": 10.25103648
              },
              "high_wealth": {
                "alpha": 0.11111111,
                "consumption": 9.10576291,
                "cash": 200.0
              }
            }
          },
          {
            "age": 99,
            "phase": "retired",
            "mid_wealth_alpha": 0.22222222,
            "mid_wealth_consumption": 5.53293779,
            "wealth_bands": {
              "low_wealth": {
                "alpha": 1.0,
                "consumption": 0.25,
                "cash": 0.25
              },
              "mid_wealth": {
                "alpha": 0.22222222,
                "consumption": 5.53293779,
                "cash": 10.25103648
              },
              "high_wealth": {
                "alpha": 0.22222222,
                "consumption": 100.71937362,
                "cash": 200.0
              }
            }
          }
        ]
      }
    }
  ]
}