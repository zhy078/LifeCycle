import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

import optimize


def _write_policy_file(path: Path, alpha_mid: float, consumption_mid: float):
    lines = [
        f"{alpha_mid - 0.1} 1.0",
        f"{alpha_mid} 2.0",
        f"{consumption_mid - 0.2} 1.0",
        f"{consumption_mid} 2.0",
        "0.1 1.0",
        "0.2 2.0",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


class LifecycleAgeDifferentiationTests(unittest.TestCase):
    def test_run_real_model_exposes_age_differentiated_lifecycle_checkpoints(self):
        params = {"rho": 6.0, "delta": 0.97, "psi": 0.6, "mu": 0.03, "sigr": 0.15}
        fixed = {"tb": 20, "tr": 65, "td": 100, "r": 1.015, "nsim": 1000}

        with tempfile.TemporaryDirectory() as tmp:
            artifact_dir = Path(tmp) / "scenario_0001"

            def fake_run_octave(cmd, cwd, timeout):
                cwd_path = Path(cwd)
                Path(cwd, "CWY.txt").write_text("0.1 42.5 0.3 1.0\n", encoding="utf-8")
                for age in range(20, 100):
                    if age < 40:
                        alpha_mid = 0.80
                        consumption_mid = 1.10
                    elif age < 65:
                        alpha_mid = 0.55
                        consumption_mid = 1.00
                    elif age < 80:
                        alpha_mid = 0.30
                        consumption_mid = 0.90
                    else:
                        alpha_mid = 0.15
                        consumption_mid = 0.75
                    file_idx = age - fixed["tb"] + 1
                    _write_policy_file(cwd_path / f"year{file_idx:02d}.txt", alpha_mid=alpha_mid, consumption_mid=consumption_mid)
                Path(cwd, "octave_stdout.txt").write_text("ok", encoding="utf-8")
                Path(cwd, "octave_stderr.txt").write_text("", encoding="utf-8")
                return {
                    "code": 0,
                    "timed_out": False,
                    "run_seconds": 0.123,
                    "stdout_path": str(cwd_path / "octave_stdout.txt"),
                    "stderr_path": str(cwd_path / "octave_stderr.txt"),
                }

            with patch("optimize.octave_executable", return_value="octave-cli"), patch("optimize.run_octave_command", side_effect=fake_run_octave):
                result = optimize.run_real_model(params, fixed, artifact_dir, fast_mode=True, timeout_sec=60)

        self.assertEqual(result["status"], "ok")
        self.assertEqual(result["year_files_count"], 80)
        self.assertIsNotNone(result["policy_summary_path"])
        self.assertIsNotNone(result["lifecycle_checkpoints"])

        checkpoints = result["lifecycle_checkpoints"]["checkpoints"]
        self.assertEqual([checkpoint["age"] for checkpoint in checkpoints], [20, 64, 65, 99])
        self.assertEqual(checkpoints[0]["age"], 20)
        self.assertEqual(checkpoints[0]["phase"], "working")
        self.assertEqual(checkpoints[1]["phase"], "working")
        self.assertEqual(checkpoints[2]["phase"], "retired")
        self.assertEqual(checkpoints[-1]["phase"], "retired")
        self.assertGreater(checkpoints[0]["mid_wealth_alpha"], checkpoints[-1]["mid_wealth_alpha"])
        self.assertGreater(checkpoints[0]["mid_wealth_consumption"], checkpoints[-1]["mid_wealth_consumption"])

    def test_partial_year_files_keep_file_index_age_mapping(self):
        fixed = {"tb": 20, "tr": 65, "td": 100}
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            _write_policy_file(run_dir / "year47.txt", alpha_mid=0.30, consumption_mid=0.90)
            _write_policy_file(run_dir / "year80.txt", alpha_mid=0.15, consumption_mid=0.75)

            summary = optimize.build_policy_summary(run_dir, fixed)

        self.assertEqual(summary["age_span"], {"start": 66, "end": 99})
        self.assertEqual([row["age"] for row in summary["per_age"]], [66, 99])
        self.assertTrue(all(row["phase"] == "retired" for row in summary["per_age"]))


if __name__ == "__main__":
    unittest.main()
