# LifeCycle Optimization Report

{
  "objective": "maximize_lifetime_utility",
  "search_method": "grid",
  "use_real_model": true,
  "fast_mode": true,
  "total_scenarios": 3,
  "failed_scenarios": 0,
  "total_model_runtime_seconds": 9497.51,
  "elapsed_seconds": 9497.595,
  "best": {
    "id": 0,
    "params": {
      "rho": 6.0,
      "delta": 0.97,
      "psi": 0.6,
      "mu": 0.03,
      "sigr": 0.15
    },
    "objective": "maximize_lifetime_utility",
    "score": 2.796848847579729,
    "status": "ok",
    "artifact_dir": "C:\\Users\\haoyu\\Desktop\\code\\github\\LifeCycle\\outputs\\grid-test-40x40-results\\scenario_0000",
    "metric": 2.796848847579729,
    "run_seconds": 3214.619,
    "year_files_count": 80,
    "octave_command": "E:\\Octave-11.1.0\\mingw64\\bin\\octave-cli.EXE --quiet life_cycle.m",
    "code": 0,
    "timed_out": false,
    "kill_meta": null,
    "diagnostic_path": "C:\\Users\\haoyu\\Desktop\\code\\github\\LifeCycle\\outputs\\grid-test-40x40-results\\scenario_0000\\runner_diagnostics.json",
    "policy_summary_path": "C:\\Users\\haoyu\\Desktop\\code\\github\\LifeCycle\\outputs\\grid-test-40x40-results\\scenario_0000\\policy_summary.json",
    "lifecycle_checkpoints": {
      "age_span": {
        "start": 20,
        "end": 99
      },
      "checkpoints": [
        {
          "age": 20,
          "phase": "working",
          "mid_wealth_alpha": 1.0,
          "mid_wealth_consumption": 1.20688351,
          "wealth_bands": {
            "low_wealth": {
              "alpha": 1.0,
              "consumption": 0.24975,
              "cash": 0.25
            },
            "mid_wealth": {
              "alpha": 1.0,
              "consumption": 1.20688351,
              "cash": 7.70378413
            },
            "high_wealth": {
              "alpha": 0.20512821,
              "consumption": 10.57725394,
              "cash": 200.0
            }
          }
        },
        {
          "age": 64,
          "phase": "working",
          "mid_wealth_alpha": 0.76923077,
          "mid_wealth_consumption": 1.1408215,
          "wealth_bands": {
            "low_wealth": {
              "alpha": 1.0,
              "consumption": 0.24975,
              "cash": 0.25
            },
            "mid_wealth": {
              "alpha": 0.76923077,
              "consumption": 1.1408215,
              "cash": 7.70378413
            },
            "high_wealth": {
              "alpha": 0.15384615,
              "consumption": 11.92963978,
              "cash": 200.0
            }
          }
        },
        {
          "age": 65,
          "phase": "retired",
          "mid_wealth_alpha": 0.76923077,
          "mid_wealth_consumption": 1.14270744,
          "wealth_bands": {
            "low_wealth": {
              "alpha": 1.0,
              "consumption": 0.25,
              "cash": 0.25
            },
            "mid_wealth": {
              "alpha": 0.76923077,
              "consumption": 1.14270744,
              "cash": 7.70378413
            },
            "high_wealth": {
              "alpha": 0.15384615,
              "consumption": 10.11816019,
              "cash": 200.0
            }
          }
        },
        {
          "age": 99,
          "phase": "retired",
          "mid_wealth_alpha": 0.28205128,
          "mid_wealth_consumption": 4.29816318,
          "wealth_bands": {
            "low_wealth": {
              "alpha": 1.0,
              "consumption": 0.25,
              "cash": 0.25
            },
            "mid_wealth": {
              "alpha": 0.28205128,
              "consumption": 4.29816318,
              "cash": 7.70378413
            },
            "high_wealth": {
              "alpha": 0.23076923,
              "consumption": 102.82477277,
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
        "rho": 6.0,
        "delta": 0.97,
        "psi": 0.6,
        "mu": 0.03,
        "sigr": 0.15
      },
      "objective": "maximize_lifetime_utility",
      "score": 2.796848847579729,
      "status": "ok",
      "artifact_dir": "C:\\Users\\haoyu\\Desktop\\code\\github\\LifeCycle\\outputs\\grid-test-40x40-results\\scenario_0000",
      "metric": 2.796848847579729,
      "run_seconds": 3214.619,
      "year_files_count": 80,
      "octave_command": "E:\\Octave-11.1.0\\mingw64\\bin\\octave-cli.EXE --quiet life_cycle.m",
      "code": 0,
      "timed_out": false,
      "kill_meta": null,
      "diagnostic_path": "C:\\Users\\haoyu\\Desktop\\code\\github\\LifeCycle\\outputs\\grid-test-40x40-results\\scenario_0000\\runner_diagnostics.json",
      "policy_summary_path": "C:\\Users\\haoyu\\Desktop\\code\\github\\LifeCycle\\outputs\\grid-test-40x40-results\\scenario_0000\\policy_summary.json",
      "lifecycle_checkpoints": {
        "age_span": {
          "start": 20,
          "end": 99
        },
        "checkpoints": [
          {
            "age": 20,
            "phase": "working",
            "mid_wealth_alpha": 1.0,
            "mid_wealth_consumption": 1.20688351,
            "wealth_bands": {
              "low_wealth": {
                "alpha": 1.0,
                "consumption": 0.24975,
                "cash": 0.25
              },
              "mid_wealth": {
                "alpha": 1.0,
                "consumption": 1.20688351,
                "cash": 7.70378413
              },
              "high_wealth": {
                "alpha": 0.20512821,
                "consumption": 10.57725394,
                "cash": 200.0
              }
            }
          },
          {
            "age": 64,
            "phase": "working",
            "mid_wealth_alpha": 0.76923077,
            "mid_wealth_consumption": 1.1408215,
            "wealth_bands": {
              "low_wealth": {
                "alpha": 1.0,
                "consumption": 0.24975,
                "cash": 0.25
              },
              "mid_wealth": {
                "alpha": 0.76923077,
                "consumption": 1.1408215,
                "cash": 7.70378413
              },
              "high_wealth": {
                "alpha": 0.15384615,
                "consumption": 11.92963978,
                "cash": 200.0
              }
            }
          },
          {
            "age": 65,
            "phase": "retired",
            "mid_wealth_alpha": 0.76923077,
            "mid_wealth_consumption": 1.14270744,
            "wealth_bands": {
              "low_wealth": {
                "alpha": 1.0,
                "consumption": 0.25,
                "cash": 0.25
              },
              "mid_wealth": {
                "alpha": 0.76923077,
                "consumption": 1.14270744,
                "cash": 7.70378413
              },
              "high_wealth": {
                "alpha": 0.15384615,
                "consumption": 10.11816019,
                "cash": 200.0
              }
            }
          },
          {
            "age": 99,
            "phase": "retired",
            "mid_wealth_alpha": 0.28205128,
            "mid_wealth_consumption": 4.29816318,
            "wealth_bands": {
              "low_wealth": {
                "alpha": 1.0,
                "consumption": 0.25,
                "cash": 0.25
              },
              "mid_wealth": {
                "alpha": 0.28205128,
                "consumption": 4.29816318,
                "cash": 7.70378413
              },
              "high_wealth": {
                "alpha": 0.23076923,
                "consumption": 102.82477277,
                "cash": 200.0
              }
            }
          }
        ]
      }
    },
    {
      "id": 1,
      "params": {
        "rho": 8.0,
        "delta": 0.97,
        "psi": 0.6,
        "mu": 0.03,
        "sigr": 0.15
      },
      "objective": "maximize_lifetime_utility",
      "score": 2.594380368908515,
      "status": "ok",
      "artifact_dir": "C:\\Users\\haoyu\\Desktop\\code\\github\\LifeCycle\\outputs\\grid-test-40x40-results\\scenario_0001",
      "metric": 2.594380368908515,
      "run_seconds": 3090.943,
      "year_files_count": 80,
      "octave_command": "E:\\Octave-11.1.0\\mingw64\\bin\\octave-cli.EXE --quiet life_cycle.m",
      "code": 0,
      "timed_out": false,
      "kill_meta": null,
      "diagnostic_path": "C:\\Users\\haoyu\\Desktop\\code\\github\\LifeCycle\\outputs\\grid-test-40x40-results\\scenario_0001\\runner_diagnostics.json",
      "policy_summary_path": "C:\\Users\\haoyu\\Desktop\\code\\github\\LifeCycle\\outputs\\grid-test-40x40-results\\scenario_0001\\policy_summary.json",
      "lifecycle_checkpoints": {
        "age_span": {
          "start": 20,
          "end": 99
        },
        "checkpoints": [
          {
            "age": 20,
            "phase": "working",
            "mid_wealth_alpha": 0.74358974,
            "mid_wealth_consumption": 1.0842095,
            "wealth_bands": {
              "low_wealth": {
                "alpha": 1.0,
                "consumption": 0.24975,
                "cash": 0.25
              },
              "mid_wealth": {
                "alpha": 0.74358974,
                "consumption": 1.0842095,
                "cash": 7.70378413
              },
              "high_wealth": {
                "alpha": 0.20512821,
                "consumption": 10.99404152,
                "cash": 200.0
              }
            }
          },
          {
            "age": 64,
            "phase": "working",
            "mid_wealth_alpha": 0.58974359,
            "mid_wealth_consumption": 1.10989126,
            "wealth_bands": {
              "low_wealth": {
                "alpha": 1.0,
                "consumption": 0.24975,
                "cash": 0.25
              },
              "mid_wealth": {
                "alpha": 0.58974359,
                "consumption": 1.10989126,
                "cash": 7.70378413
              },
              "high_wealth": {
                "alpha": 0.15384615,
                "consumption": 11.71234274,
                "cash": 200.0
              }
            }
          },
          {
            "age": 65,
            "phase": "retired",
            "mid_wealth_alpha": 0.61538462,
            "mid_wealth_consumption": 1.13783747,
            "wealth_bands": {
              "low_wealth": {
                "alpha": 1.0,
                "consumption": 0.25,
                "cash": 0.25
              },
              "mid_wealth": {
                "alpha": 0.61538462,
                "consumption": 1.13783747,
                "cash": 7.70378413
              },
              "high_wealth": {
                "alpha": 0.12820513,
                "consumption": 10.10263368,
                "cash": 200.0
              }
            }
          },
          {
            "age": 99,
            "phase": "retired",
            "mid_wealth_alpha": 0.20512821,
            "mid_wealth_consumption": 4.27171898,
            "wealth_bands": {
              "low_wealth": {
                "alpha": 1.0,
                "consumption": 0.25,
                "cash": 0.25
              },
              "mid_wealth": {
                "alpha": 0.20512821,
                "consumption": 4.27171898,
                "cash": 7.70378413
              },
              "high_wealth": {
                "alpha": 0.17948718,
                "consumption": 102.22530482,
                "cash": 200.0
              }
            }
          }
        ]
      }
    },
    {
      "id": 2,
      "params": {
        "rho": 10.0,
        "delta": 0.97,
        "psi": 0.6,
        "mu": 0.03,
        "sigr": 0.15
      },
      "objective": "maximize_lifetime_utility",
      "score": 2.4732772076076097,
      "status": "ok",
      "artifact_dir": "C:\\Users\\haoyu\\Desktop\\code\\github\\LifeCycle\\outputs\\grid-test-40x40-results\\scenario_0002",
      "metric": 2.4732772076076097,
      "run_seconds": 3191.948,
      "year_files_count": 80,
      "octave_command": "E:\\Octave-11.1.0\\mingw64\\bin\\octave-cli.EXE --quiet life_cycle.m",
      "code": 0,
      "timed_out": false,
      "kill_meta": null,
      "diagnostic_path": "C:\\Users\\haoyu\\Desktop\\code\\github\\LifeCycle\\outputs\\grid-test-40x40-results\\scenario_0002\\runner_diagnostics.json",
      "policy_summary_path": "C:\\Users\\haoyu\\Desktop\\code\\github\\LifeCycle\\outputs\\grid-test-40x40-results\\scenario_0002\\policy_summary.json",
      "lifecycle_checkpoints": {
        "age_span": {
          "start": 20,
          "end": 99
        },
        "checkpoints": [
          {
            "age": 20,
            "phase": "working",
            "mid_wealth_alpha": 0.53846154,
            "mid_wealth_consumption": 0.95872319,
            "wealth_bands": {
              "low_wealth": {
                "alpha": 1.0,
                "consumption": 0.24975,
                "cash": 0.25
              },
              "mid_wealth": {
                "alpha": 0.53846154,
                "consumption": 0.95872319,
                "cash": 7.70378413
              },
              "high_wealth": {
                "alpha": 0.15384615,
                "consumption": 12.23266879,
                "cash": 200.0
              }
            }
          },
          {
            "age": 64,
            "phase": "working",
            "mid_wealth_alpha": 0.48717949,
            "mid_wealth_consumption": 1.09175058,
            "wealth_bands": {
              "low_wealth": {
                "alpha": 1.0,
                "consumption": 0.24975,
                "cash": 0.25
              },
              "mid_wealth": {
                "alpha": 0.48717949,
                "consumption": 1.09175058,
                "cash": 7.70378413
              },
              "high_wealth": {
                "alpha": 0.12820513,
                "consumption": 11.31579026,
                "cash": 200.0
              }
            }
          },
          {
            "age": 65,
            "phase": "retired",
            "mid_wealth_alpha": 0.48717949,
            "mid_wealth_consumption": 1.1294862,
            "wealth_bands": {
              "low_wealth": {
                "alpha": 1.0,
                "consumption": 0.25,
                "cash": 0.25
              },
              "mid_wealth": {
                "alpha": 0.48717949,
                "consumption": 1.1294862,
                "cash": 7.70378413
              },
              "high_wealth": {
                "alpha": 0.12820513,
                "consumption": 9.8563346,
                "cash": 200.0
              }
            }
          },
          {
            "age": 99,
            "phase": "retired",
            "mid_wealth_alpha": 0.15384615,
            "mid_wealth_consumption": 4.25900556,
            "wealth_bands": {
              "low_wealth": {
                "alpha": 1.0,
                "consumption": 0.25,
                "cash": 0.25
              },
              "mid_wealth": {
                "alpha": 0.15384615,
                "consumption": 4.25900556,
                "cash": 7.70378413
              },
              "high_wealth": {
                "alpha": 0.12820513,
                "consumption": 102.0062429,
                "cash": 200.0
              }
            }
          }
        ]
      }
    }
  ]
}