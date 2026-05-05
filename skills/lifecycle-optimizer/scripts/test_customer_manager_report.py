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
            "metric": 0.05543808,
            "status": "ok",
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
        self.assertNotIn("瀹", rendered)
        self.assertNotIn("閫", rendered)

    def test_timeout_report_is_readable_and_marks_preview(self):
        best = {
            "objective": "maximize_terminal_wealth",
            "score": None,
            "metric": 7313.93913985,
            "status": "error_timeout",
            "params": {"rho": 6.0, "delta": 0.97, "psi": 0.6, "mu": 0.03, "sigr": 0.15},
            "policy_summary_path": None,
            "lifecycle_checkpoints": {
                "age_span": {"start": 66, "end": 99},
                "checkpoints": [
                    {
                        "age": 66,
                        "phase": "retired",
                        "mid_wealth_alpha": 0.44,
                        "mid_wealth_consumption": 3.86,
                        "wealth_bands": {
                            "low_wealth": {"alpha": 1.0, "consumption": 0.25, "cash": 0.25},
                            "mid_wealth": {"alpha": 0.44, "consumption": 3.86, "cash": 10.25},
                            "high_wealth": {"alpha": 0.11, "consumption": 9.1, "cash": 200.0},
                        },
                    }
                ],
            },
        }
        report_meta = {"objective": "maximize_terminal_wealth", "total_scenarios": 1, "elapsed_seconds": 5.1}
        fixed = {"tb": 20, "tr": 65, "td": 100}

        rendered = optimize.render_customer_manager_report(best, report_meta, fixed)

        self.assertIn("模型状态: error_timeout", rendered)
        self.assertIn("最优场景分数: NA", rendered)
        self.assertIn("部分模型 metric: 7313.939140", rendered)
        self.assertIn("诊断/预览", rendered)
        self.assertIn("66岁到100岁", rendered)
        self.assertNotIn("瀹", rendered)
        self.assertNotIn("閫", rendered)


if __name__ == "__main__":
    unittest.main()
