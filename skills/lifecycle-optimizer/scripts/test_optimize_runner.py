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

            def fake_run(cmd, cwd, capture_output, text, encoding, errors, timeout):
                Path(cwd, "CWY.txt").write_text("0.1 42.5 0.3 1.0\n", encoding="utf-8")
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

                class Result:
                    returncode = 0
                    stdout = "ok"
                    stderr = ""

                return Result()

            with patch("optimize.subprocess.run", side_effect=fake_run) as run_mock:
                result = optimize.run_real_model(params, fixed, artifact_dir, fast_mode=True, timeout_sec=321)

            called_kwargs = run_mock.call_args.kwargs
            self.assertEqual(called_kwargs["cwd"], str(artifact_dir))
            self.assertEqual(called_kwargs["timeout"], 321)
            self.assertEqual(result["status"], "ok")
            self.assertEqual(result["year_files_count"], 1)
            self.assertEqual(result["metric"], 42.5)
            self.assertTrue((artifact_dir / "life_cycle.m").exists())
            self.assertTrue((artifact_dir / "year80.txt").exists())
            self.assertTrue((artifact_dir / "octave_stdout.txt").exists())
            self.assertTrue((artifact_dir / "policy_summary.json").exists())
            rendered = (artifact_dir / "life_cycle.m").read_text(encoding="utf-8")
            summary = (artifact_dir / "policy_summary.json").read_text(encoding="utf-8")
            self.assertIn("na = 10;", rendered)
            self.assertIn("ncash = 10;", rendered)
            self.assertIn("nsim = 1000;", rendered)
            self.assertIn("n      = 5;", rendered)
            self.assertIn('"age": 20', summary)
            self.assertIn('"average_mid_wealth_alpha": 0.2', summary)


if __name__ == "__main__":
    unittest.main()
