import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

import optimize


class RunRealModelTests(unittest.TestCase):
    def test_run_real_model_executes_inside_artifact_dir(self):
        params = {"rho": 6.0, "delta": 0.97, "psi": 0.6, "mu": 0.03, "sigr": 0.15}
        fixed = {"tb": 20, "tr": 66, "td": 100, "r": 1.015, "nsim": 1000}

        with tempfile.TemporaryDirectory() as tmp:
            artifact_dir = Path(tmp) / "scenario_0000"

            def fake_run_octave(cmd, cwd, timeout):
                Path(cwd, "CWY.txt").write_text("0.1 42.5 0.3 1.0\n", encoding="utf-8")
                Path(cwd, "year01.txt").write_text(
                    "\n".join([
                        "0.1 1.0",
                        "0.2 2.0",
                        "1.1 1.0",
                        "1.2 2.0",
                        "2.1 1.0",
                        "2.2 2.0",
                    ]) + "\n",
                    encoding="utf-8",
                )
                Path(cwd, "octave_stdout.txt").write_text("ok", encoding="utf-8")
                Path(cwd, "octave_stderr.txt").write_text("", encoding="utf-8")
                return {
                    "code": 0,
                    "timed_out": False,
                    "run_seconds": 0.123,
                    "stdout_path": str(Path(cwd, "octave_stdout.txt")),
                    "stderr_path": str(Path(cwd, "octave_stderr.txt")),
                }

            with patch("optimize.octave_executable", return_value="octave-cli"), patch("optimize.run_octave_command", side_effect=fake_run_octave) as run_mock:
                result = optimize.run_real_model(params, fixed, artifact_dir, fast_mode=True, timeout_sec=321)

            called_cmd, called_dir, called_timeout = run_mock.call_args.args
            self.assertTrue(str(called_cmd[0]).lower().endswith(("octave-cli", "octave-cli.exe", "octave", "octave.exe")))
            self.assertEqual(called_dir, artifact_dir)
            self.assertEqual(called_timeout, 321)
            self.assertEqual(result["status"], "ok")
            self.assertEqual(result["year_files_count"], 1)
            self.assertEqual(result["metric"], 2.2)
            self.assertTrue((artifact_dir / "life_cycle.m").exists())
            self.assertTrue((artifact_dir / "year01.txt").exists())
            self.assertTrue((artifact_dir / "octave_stdout.txt").exists())
            self.assertTrue((artifact_dir / "runner_diagnostics.json").exists())
            self.assertTrue((artifact_dir / "policy_summary.json").exists())
            rendered = (artifact_dir / "life_cycle.m").read_text(encoding="utf-8")
            summary = (artifact_dir / "policy_summary.json").read_text(encoding="utf-8")
            self.assertIn("na = 10;", rendered)
            self.assertIn("ncash = 10;", rendered)
            self.assertIn("nsim = 1000;", rendered)
            self.assertIn("n      = 5;", rendered)
            self.assertIn('"age": 20', summary)
            self.assertIn('"average_mid_wealth_alpha": 0.2', summary)

    def test_run_real_model_returns_timeout_with_diagnostics(self):
        params = {"rho": 6.0, "delta": 0.97, "psi": 0.6, "mu": 0.03, "sigr": 0.15}
        fixed = {"tb": 20, "tr": 66, "td": 100, "r": 1.015, "nsim": 1000}

        with tempfile.TemporaryDirectory() as tmp:
            artifact_dir = Path(tmp) / "scenario_0000"

            def fake_timeout(cmd, cwd, timeout):
                Path(cwd, "year80.txt").write_text(
                    "\n".join([
                        "0.1 1.0",
                        "0.2 2.0",
                        "1.1 1.0",
                        "1.2 2.0",
                        "2.1 1.0",
                        "2.2 2.0",
                    ]) + "\n",
                    encoding="utf-8",
                )
                Path(cwd, "octave_stdout.txt").write_text("", encoding="utf-8")
                Path(cwd, "octave_stderr.txt").write_text("timeout", encoding="utf-8")
                return {
                    "code": None,
                    "timed_out": True,
                    "run_seconds": 1.0,
                    "stdout_path": str(Path(cwd, "octave_stdout.txt")),
                    "stderr_path": str(Path(cwd, "octave_stderr.txt")),
                }

            with patch("optimize.octave_executable", return_value="octave-cli"), patch("optimize.run_octave_command", side_effect=fake_timeout):
                result = optimize.run_real_model(params, fixed, artifact_dir, fast_mode=True, timeout_sec=1)

            self.assertEqual(result["status"], "error_timeout")
            self.assertTrue(result["timed_out"])
            self.assertEqual(result["year_files_count"], 1)
            self.assertTrue((artifact_dir / "runner_diagnostics.json").exists())
            self.assertTrue((artifact_dir / "policy_summary.json").exists())

    def test_fast_mode_caps_large_config_nsim(self):
        params = {"rho": 6.0, "delta": 0.97, "psi": 0.6, "mu": 0.03, "sigr": 0.15}
        fixed = {"tb": 20, "tr": 66, "td": 100, "r": 1.015, "nsim": 10000}

        with tempfile.TemporaryDirectory() as tmp:
            dst = Path(tmp) / "life_cycle.m"
            optimize.patch_lifecycle_script(optimize.model_dir() / "life_cycle.m", dst, params, fixed, fast_mode=True)

            rendered = dst.read_text(encoding="utf-8")
            self.assertIn("na = 10;", rendered)
            self.assertIn("ncash = 10;", rendered)
            self.assertIn("nsim = 1000;", rendered)

    def test_run_model_clears_stale_scenario_artifacts(self):
        params = {"rho": 6.0, "delta": 0.97, "psi": 0.6, "mu": 0.03, "sigr": 0.15}
        fixed = {"tb": 20, "tr": 66, "td": 100, "r": 1.015, "nsim": 1000}

        with tempfile.TemporaryDirectory() as tmp:
            artifact_dir = Path(tmp) / "scenario_0000"
            artifact_dir.mkdir()
            (artifact_dir / "year80.txt").write_text("stale\n", encoding="utf-8")

            result = optimize.run_model(params, fixed, artifact_dir, dry_run=True, use_real_model=False, fast_mode=False, timeout_sec=1)

            self.assertEqual(result["status"], "dry_run")
            self.assertFalse((artifact_dir / "year80.txt").exists())


if __name__ == "__main__":
    unittest.main()
