import unittest

from advisor import build_advice
from client_input import map_risk_preference, parse_client_input
from policy_guardrails import apply_policy_guardrails
from risk_metrics import estimate_risk_metrics


class ClientInputTests(unittest.TestCase):
    def test_parse_chinese_client_input(self):
        parsed = parse_client_input("45岁，100万，偏稳健")
        self.assertEqual(parsed["age"], 45)
        self.assertEqual(parsed["wealth"], 10.0)
        self.assertEqual(parsed["risk_preference"], "conservative")
        self.assertEqual(parsed["rho"], 10)

    def test_parse_english_stable_returns(self):
        label, rho = map_risk_preference("I want stable returns")
        self.assertEqual(label, "conservative")
        self.assertEqual(rho, 10)

    def test_parse_english_thousands(self):
        parsed = parse_client_input("65 years old, 100 thousands, risky returns")
        self.assertEqual(parsed["age"], 65)
        self.assertEqual(parsed["wealth"], 1.0)
        self.assertEqual(parsed["risk_preference"], "aggressive")


class RiskMetricsTests(unittest.TestCase):
    def test_risk_metrics_are_valid(self):
        metrics = estimate_risk_metrics(alpha=0.4, wealth=10, paths=100, years=5)
        self.assertLessEqual(metrics["max_drawdown_pct"], 0)
        self.assertGreaterEqual(metrics["failure_probability_pct"], 0)
        self.assertGreaterEqual(metrics["failure_probability_margin_pct"], 0)
        self.assertGreater(metrics["annual_volatility_pct"], 0)


class PolicyGuardrailTests(unittest.TestCase):
    def test_low_wealth_retirement_case_is_overridden(self):
        client = {"age": 65, "wealth": 1.0}

        def risk_for_alpha(alpha):
            return estimate_risk_metrics(alpha=alpha, wealth=1.0, years=20)

        model_risk = risk_for_alpha(1.0)
        decision = apply_policy_guardrails(client, 1.0, model_risk, risk_for_alpha)
        self.assertEqual(decision["case_type"], "UNSAFE_MODEL_OUTPUT")
        self.assertEqual(decision["status"], "adjusted")
        self.assertAlmostEqual(decision["final_alpha"], 0.4)
        self.assertLess(decision["final_alpha"], decision["model_alpha"])
        self.assertLessEqual(decision["final_risk"]["failure_probability_pct"], decision["failure_limit_pct"])
        self.assertGreater(len(decision["failure_search"]), 1)
        self.assertAlmostEqual(decision["failure_search"][0]["alpha"], 1.0)
        self.assertAlmostEqual(decision["failure_search"][-1]["alpha"], 0.55)
        self.assertTrue(decision["failure_search"][-1]["passed"])
        self.assertLessEqual(
            decision["failure_search"][-1]["failure_probability_upper_bound_pct"],
            decision["failure_limit_pct"],
        )
        self.assertEqual([step["name"] for step in decision["decision_trace"]], [
            "failure_risk",
            "age_constraint",
            "wealth_constraint",
        ])


class AdvisorTests(unittest.TestCase):
    def test_build_advice_from_existing_20x20_results(self):
        result = build_advice("45岁，100万，偏稳健")
        self.assertTrue(result["ok"])
        self.assertIn("Risk Control Decision", result["markdown"])
        self.assertIn("rho=10", result["markdown"])

    def test_unsafe_client_output_is_marked_and_overridden(self):
        result = build_advice("65 years old, 100 thousands, risky returns")
        self.assertTrue(result["ok"])
        self.assertEqual(result["decision"]["case_type"], "UNSAFE_MODEL_OUTPUT")
        self.assertEqual(result["decision"]["status"], "adjusted")
        self.assertLess(result["policy"]["alpha"], result["policy"]["model_alpha"])
        self.assertIn("Decision Trace", result["markdown"])
        self.assertIn("Failure Risk Search", result["markdown"])
        self.assertIn("Baseline Comparison", result["markdown"])
        self.assertIn("Assumed planning horizon: 20 years (age 65-85)", result["markdown"])
        self.assertIn("Although 55% satisfies the failure constraint", result["markdown"])
        self.assertIn("Suggested allocation range: 30%-45%", result["markdown"])
        self.assertIn("Applicability", result["markdown"])
        self.assertIn("95% MC margin", result["markdown"])
        self.assertIn("avoiding large losses", result["markdown"])
        self.assertIn("Model output adjusted", result["markdown"])


if __name__ == "__main__":
    unittest.main()
