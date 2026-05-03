import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import optimize


class CustomerManagerReportTests(unittest.TestCase):
    def test_render_customer_manager_report_uses_current_info_and_decade_suggestions(self):
        best = {
            "objective": "maximize_terminal_wealth",
            "score": 0.05543808,
            "params": {"rho": 8.0, "delta": 0.95, "psi": 0.6, "mu": 0.03, "sigr": 0.15},
            "policy_summary_path": None,
            "lifecycle_checkpoints": {
                "age_span": {"start": 20, "end": 99},
                "checkpoints": [
                    {
                        "age": 20,
                        "phase": "working",
                        "mid_wealth_alpha": 0.44,
                        "mid_wealth_consumption": 4.14,
                        "wealth_bands": {
                            "low_wealth": {"alpha": 1.0, "consumption": 0.25, "cash": 0.25},
                            "mid_wealth": {"alpha": 0.44, "consumption": 4.14, "cash": 10.25},
                            "high_wealth": {"alpha": 0.11, "consumption": 15.5, "cash": 200.0},
                        },
                    },
                    {
                        "age": 30,
                        "phase": "working",
                        "mid_wealth_alpha": 0.33,
                        "mid_wealth_consumption": 4.00,
                        "wealth_bands": {
                            "low_wealth": {"alpha": 0.89, "consumption": 0.30, "cash": 0.25},
                            "mid_wealth": {"alpha": 0.33, "consumption": 4.00, "cash": 10.25},
                            "high_wealth": {"alpha": 0.22, "consumption": 14.2, "cash": 200.0},
                        },
                    },
                    {
                        "age": 40,
                        "phase": "working",
                        "mid_wealth_alpha": 0.33,
                        "mid_wealth_consumption": 3.80,
                        "wealth_bands": {
                            "low_wealth": {"alpha": 0.78, "consumption": 0.35, "cash": 0.25},
                            "mid_wealth": {"alpha": 0.33, "consumption": 3.80, "cash": 10.25},
                            "high_wealth": {"alpha": 0.22, "consumption": 13.1, "cash": 200.0},
                        },
                    },
                ],
            },
        }
        report_meta = {"objective": "maximize_terminal_wealth", "total_scenarios": 1, "elapsed_seconds": 183.7}
        fixed = {"tb": 20, "tr": 65, "td": 100}
        client_profile = {"current_age": 31, "current_wealth": 10.25}

        rendered = optimize.render_customer_manager_report(best, report_meta, fixed, client_profile)

        self.assertIn("## 客户当前信息", rendered)
        self.assertIn("## 后续建议", rendered)
        self.assertIn("31岁", rendered)
        self.assertIn("current_wealth=10.25", rendered)
        self.assertIn("20岁到30岁", rendered)
        self.assertIn("30岁到40岁", rendered)
        self.assertIn("平衡型组合", rendered)
        self.assertNotIn("退休", rendered)


if __name__ == "__main__":
    unittest.main()
