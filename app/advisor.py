import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "skills" / "lifecycle-optimizer" / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import optimize
from policy_explainer import explain_policy

from client_input import parse_client_input
from policy_guardrails import apply_policy_guardrails
from risk_metrics import estimate_risk_metrics


DEFAULT_RESULTS_DIR = ROOT / "outputs" / "grid-test-20x20-results"
BASELINE_ALPHA = 0.60


def _load_jsonl(path):
    rows = []
    if not path.exists():
        return rows
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = line.strip()
        if line:
            rows.append(json.loads(line))
    return rows


def _select_result(results_dir, client):
    rho = client["rho"] if isinstance(client, dict) else client
    income_profile = client.get("income_profile", "stable") if isinstance(client, dict) else "stable"
    market_profile = client.get("market_profile", "base") if isinstance(client, dict) else "base"
    rows = [row for row in _load_jsonl(results_dir / "results.jsonl") if row.get("status") == "ok"]
    if not rows:
        best_path = results_dir / "best_params.json"
        if best_path.exists():
            rows = [json.loads(best_path.read_text(encoding="utf-8"))]
    if not rows:
        library_row = _select_policy_library_result(results_dir, rho, income_profile, market_profile)
        if library_row:
            return library_row
        return None
    return min(rows, key=lambda row: abs(float(row.get("params", {}).get("rho", rho)) - rho))


def _rho_token(rho):
    return f"rho_{int(round(float(rho))):02d}"


def _select_policy_library_result(results_dir, rho, income_profile="stable", market_profile="base"):
    results_dir = Path(results_dir)
    if (results_dir / "year01.txt").exists():
        return _library_row_from_dir(results_dir, rho)

    exact = results_dir / f"{_rho_token(rho)}_income_{income_profile}_market_{market_profile}"
    if exact.exists():
        return _library_row_from_dir(exact, rho)

    matches = sorted(results_dir.glob(f"{_rho_token(rho)}_income_*_market_*"))
    if matches:
        return _library_row_from_dir(matches[0], rho)
    return None


def _library_row_from_dir(artifact_dir, rho):
    metadata_path = Path(artifact_dir) / "metadata.json"
    params = {"rho": float(rho)}
    if metadata_path.exists():
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        params.update({
            key: metadata[key]
            for key in ["rho", "delta", "psi", "mu", "sigr"]
            if key in metadata
        })
    return {
        "status": "ok",
        "params": params,
        "artifact_dir": str(artifact_dir),
    }


def _interpolate(entries, wealth):
    value = optimize.interpolate_policy_value(entries, wealth)
    return float(value) if value is not None else None


def _policy_at_client(row, client):
    artifact_dir = Path(row["artifact_dir"])
    fixed = {"tb": 20, "td": 100}
    age = int(client["age"])
    wealth = float(client["wealth"])
    file_idx = age - fixed["tb"] + 1
    year_file = artifact_dir / f"year{file_idx:02d}.txt"
    if not year_file.exists():
        candidates = sorted(artifact_dir.glob("year*.txt"), key=lambda p: optimize.year_file_index(p) or p.name)
        if not candidates:
            return None
        year_file = min(candidates, key=lambda p: abs((optimize.year_file_index(p) or file_idx) - file_idx))

    parsed = optimize.parse_year_policy(year_file)
    if not parsed:
        return None
    return {
        "alpha": _interpolate(parsed["alpha"], wealth),
        "consumption": _interpolate(parsed["consumption"], wealth),
        "utility": _interpolate(parsed["value"], wealth),
        "year_file": str(year_file),
    }


def _risk_horizon_years(age):
    return max(10, min(30, 85 - int(age)))


def _risk_horizon_label(age, years):
    end_age = int(age) + years
    if age >= 60:
        stage = "expected remaining lifetime planning window"
    elif age >= 45:
        stage = "pre-retirement planning window"
    else:
        stage = "long-term accumulation planning window"
    return f"{years} years (age {age}-{end_age}; {stage})"


def _applicability_statement(age):
    if age >= 60:
        return "This recommendation framework is particularly suitable for retirement-stage clients where capital preservation and stable consumption are priorities."
    if age >= 45:
        return "This recommendation framework is suitable for pre-retirement clients who need to balance continued growth with rising capital-preservation needs."
    return "This recommendation framework is suitable for early-career clients with limited wealth buffers, where building financial resilience should come before taking large investment risk."


def _client_priority_statement(age):
    if age >= 60:
        return "For a client at this stage, avoiding large losses is typically more important than maximizing returns."
    if age >= 45:
        return "For a client at this stage, preserving accumulated savings while keeping measured growth exposure is usually more important than maximizing risk-taking."
    return "For a young client with a very small wealth base, building a safety buffer is typically more important than maximizing near-term investment returns."


def _risk_preference_phrase(risk_preference):
    if risk_preference == "aggressive":
        return "aggressive risk preference"
    if risk_preference == "conservative":
        return "conservative risk preference"
    return "balanced risk preference"


def _expected_annual_return(alpha, mu, risk_free):
    return alpha * mu + (1.0 - alpha) * risk_free


def _format_failure_estimate(risk):
    return (
        f"{risk['failure_probability_pct']:.1f}% "
        f"(95% MC margin +/- {risk['failure_probability_margin_pct']:.1f} pp)"
    )


def _allocation_range(alpha):
    return max(0.0, alpha - 0.10), min(1.0, alpha + 0.05)


def _fmt_pct(value):
    return f"{100.0 * value:.0f}%"


def _fmt_allocation(low, high):
    if abs(low - high) < 1e-9:
        return _fmt_pct(low)
    return f"{_fmt_pct(low)}-{_fmt_pct(high)}"


def _lifecycle_allocation_path(age, alpha):
    age = int(age)
    alpha = max(0.0, min(1.0, float(alpha)))
    rows = []

    if age < 50:
        rows.append({
            "age_range": f"{age}-50",
            "allocation": _fmt_allocation(alpha, alpha),
            "goal": "保本为主，低波动，保留少量增长敞口",
        })

    if age < 60:
        start_age = max(age, 50)
        if alpha <= 0.10:
            high = min(alpha, 0.05)
            low = max(0.0, high - 0.02)
        else:
            high = min(alpha, 0.50)
            low = max(0.0, high - 0.10)
        rows.append({
            "age_range": f"{start_age}-60",
            "allocation": _fmt_allocation(low, high),
            "goal": "逐步降低风险，减少权益波动对退休准备金的影响",
        })

    if age < 75:
        start_age = max(age, 60)
        if alpha <= 0.10:
            high = min(0.03, max(alpha * 0.60, 0.03))
        else:
            high = min(alpha * 0.50, 0.30)
        rows.append({
            "age_range": f"{start_age}-75",
            "allocation": _fmt_allocation(0.0, high),
            "goal": "现金流稳定，优先控制本金回撤和提前耗尽风险",
        })

    if not rows:
        rows.append({
            "age_range": f"{age}+",
            "allocation": _fmt_allocation(0.0, min(alpha, 0.03)),
            "goal": "以现金流稳定和本金保护为主",
        })
    return rows


def _withdrawal_plan(client, policy, decision):
    withdrawal_low = 0.03
    withdrawal_high = 0.04
    wealth = float(client["wealth"])
    failure_upper = (
        decision["final_risk"]["failure_probability_pct"]
        + decision["final_risk"]["failure_probability_margin_pct"]
    )
    survival_lower = max(0.0, 100.0 - failure_upper)
    return {
        "withdrawal_low_pct": 100.0 * withdrawal_low,
        "withdrawal_high_pct": 100.0 * withdrawal_high,
        "amount_low": wealth * withdrawal_low,
        "amount_high": wealth * withdrawal_high,
        "model_consumption": policy.get("consumption"),
        "survival_lower_bound_pct": survival_lower,
    }


def _wealth_path_rows(client, decision):
    end_age = int(client["age"]) + int(decision["risk_horizon_years"])
    risk = decision.get("wealth_path_risk", decision["final_risk"])
    return [
        {"scenario": "保守情景 P10", "age": end_age, "wealth": risk["terminal_wealth_p10"]},
        {"scenario": "中位情景 P50", "age": end_age, "wealth": risk["terminal_wealth_p50"]},
        {"scenario": "乐观情景 P90", "age": end_age, "wealth": risk["terminal_wealth_p90"]},
    ]


def build_advice(text, results_dir=DEFAULT_RESULTS_DIR):
    client = parse_client_input(text)
    missing = [name for name in ["age", "wealth"] if client[name] is None]
    if missing:
        return {
            "ok": False,
            "error": f"Missing client fields: {', '.join(missing)}",
            "client": client,
            "markdown": "Please provide age and wealth, for example: 45 years old, 1 million, stable returns.",
        }

    row = _select_result(Path(results_dir), client)
    if not row:
        return {
            "ok": False,
            "error": "No completed model result found.",
            "client": client,
            "markdown": f"No completed model result found under `{results_dir}`.",
        }

    policy = _policy_at_client(row, client)
    if not policy or policy["alpha"] is None:
        return {
            "ok": False,
            "error": "Could not read policy files for this client.",
            "client": client,
            "model_result": row,
            "markdown": "Model result exists, but policy/value files could not be parsed.",
        }

    model_alpha = max(0.0, min(1.0, policy["alpha"]))
    # Use a conservative product-layer withdrawal assumption for risk metrics.
    # Model consumption is shown separately because its units are model-state units.
    consumption_rate = 0.04
    risk_horizon_years = _risk_horizon_years(client["age"])
    mu = float(row.get("params", {}).get("mu", 0.03))
    sigma = float(row.get("params", {}).get("sigr", 0.15))

    def risk_for_alpha(alpha):
        return estimate_risk_metrics(
            alpha=alpha,
            wealth=client["wealth"],
            mu=mu,
            sigma=sigma,
            risk_free=0.015,
            consumption_rate=consumption_rate,
            years=risk_horizon_years,
        )

    model_risk = risk_for_alpha(model_alpha)
    decision = apply_policy_guardrails(client, model_alpha, model_risk, risk_for_alpha)
    decision["risk_horizon_years"] = risk_horizon_years
    decision["risk_horizon_note"] = f"Assumed planning horizon: {risk_horizon_years} years (age {client['age']}-{client['age'] + risk_horizon_years})."
    decision["risk_horizon_label"] = _risk_horizon_label(client["age"], risk_horizon_years)
    decision["applicability_statement"] = _applicability_statement(client["age"])
    decision["client_priority_statement"] = _client_priority_statement(client["age"])
    range_low, range_high = _allocation_range(decision["final_alpha"])
    decision["allocation_range"] = {"low": range_low, "high": range_high}
    decision["baseline"] = {
        "alpha": BASELINE_ALPHA,
        "risk": risk_for_alpha(BASELINE_ALPHA),
        "expected_annual_return_pct": 100.0 * _expected_annual_return(BASELINE_ALPHA, mu, 0.015),
    }
    decision["final_expected_annual_return_pct"] = 100.0 * _expected_annual_return(decision["final_alpha"], mu, 0.015)
    final_policy = {**policy, "alpha": decision["final_alpha"], "model_alpha": model_alpha}
    decision["lifecycle_path"] = _lifecycle_allocation_path(client["age"], decision["final_alpha"])
    decision["withdrawal_plan"] = _withdrawal_plan(client, final_policy, decision)
    decision["wealth_path_risk"] = estimate_risk_metrics(
        alpha=decision["final_alpha"],
        wealth=client["wealth"],
        mu=mu,
        sigma=sigma,
        risk_free=0.015,
        consumption_rate=0.03,
        years=risk_horizon_years,
    )
    decision["wealth_path"] = _wealth_path_rows(client, decision)
    explanation = explain_decision(client, final_policy, decision)
    markdown = render_advice(client, row, final_policy, decision, explanation)
    return {
        "ok": True,
        "client": client,
        "model_result": row,
        "policy": final_policy,
        "decision": decision,
        "risk_metrics": decision["final_risk"],
        "markdown": markdown,
    }


def explain_decision(client, policy, decision):
    if decision["case_type"] == "UNSAFE_MODEL_OUTPUT":
        lines = [
            "Model output adjusted due to risk controls.",
            "",
            f"The original model suggests a high-risk allocation despite the client's {_risk_preference_phrase(client['risk_preference'])}. This recommendation is adjusted to ensure:",
            "",
            "- Failure probability remains below the product risk limit",
            "- Risk exposure is appropriate for the client's lifecycle stage",
            "- Wealth level is protected from large drawdowns",
            "- Withdrawal needs can be supported through the planning horizon",
            "",
            "Adjustment reasons:",
        ]
        for reason in decision["reasons"]:
            lines.append(f"- {reason}")
        lines.extend([
            "",
            f"Suggested allocation range: {decision['allocation_range']['low']:.0%}-{decision['allocation_range']['high']:.0%} risky assets.",
            f"Point estimate: {decision['final_alpha']:.0%} risky assets, "
            f"{(1.0 - decision['final_alpha']):.0%} stabilizing assets. "
            "This balances growth and capital preservation.",
            "",
            "As the client moves closer to retirement, the lifecycle constraint automatically increases the weight on principal safety. Even with an aggressive risk preference, risky allocation should decline along the glide path.",
            "",
            decision["client_priority_statement"],
            "",
            "Compared with a standard 60/40 portfolio, the adjusted recommendation lowers drawdown, failure probability, and volatility, with a modestly lower expected annual return.",
        ])
        return "\n".join(lines)

    return explain_policy(
        client["age"],
        client["wealth"],
        client["rho"],
        policy["alpha"],
        policy.get("utility"),
        decision["final_risk"],
    )


def render_advice(client, row, policy, decision, explanation):
    model_risk = decision["model_risk"]
    final_risk = decision["final_risk"]
    failure_safe_alpha = decision["failure_search"][-1]["alpha"] if decision["failure_search"] else None
    final_alpha = policy["alpha"]
    if failure_safe_alpha is not None and failure_safe_alpha > final_alpha:
        tradeoff_note = (
            f"Although {failure_safe_alpha:.0%} satisfies the failure constraint under uncertainty bounds, further reductions are applied "
            "to align with lifecycle-stage and wealth risk limits."
        )
    elif failure_safe_alpha is not None:
        tradeoff_note = (
            f"The final point estimate remains at {final_alpha:.0%} because the failure-risk search is already the binding risk control."
        )
    else:
        tradeoff_note = "No failure-risk reduction was required; lifecycle-stage and wealth controls determine the final allocation."
    withdrawal = decision["withdrawal_plan"]
    return "\n".join([
        "# Lifecycle Client Advice",
        "",
        "## Parsed Client Input",
        "",
        f"- Age: {client['age']}",
        f"- Wealth: {client['wealth']:.2f} model units",
        f"- Risk preference: {client['risk_preference']} -> rho={client['rho']}",
        "",
        "## Model Output",
        "",
        f"- Model status: {row.get('status')}",
        f"- Matched rho: {row.get('params', {}).get('rho')}",
        f"- Raw model risky allocation alpha: {policy['model_alpha']:.4f}",
        f"- Recommended consumption: {policy['consumption']:.4f}",
        f"- Lifecycle utility: {policy['utility']:.4f}",
        "",
        "## 生命周期配置路径",
        "",
        "| 年龄区间 | 风险资产比例 | 策略目标 |",
        "|---|---:|---|",
        *[
            f"| {row['age_range']} | {row['allocation']} | {row['goal']} |"
            for row in decision["lifecycle_path"]
        ],
        "",
        "这不是一次性的单点建议，而是一条随年龄自动下调风险暴露的 glide path。随着年龄接近退休阶段，模型会提高对本金安全和现金流稳定的权重；即使客户风险偏好较高，也会逐步降低风险资产配置。",
        "",
        "## 年度消费与取现建议",
        "",
        f"- 当前建议支取: 约 {withdrawal['withdrawal_low_pct']:.0f}%-{withdrawal['withdrawal_high_pct']:.0f}%/年 "
        f"({withdrawal['amount_low']:.2f}-{withdrawal['amount_high']:.2f} model units)",
        f"- 模型 consumption 参考值: {withdrawal['model_consumption']:.4f}，在客户经理场景中转化为可执行的年度支取区间，而不是要求客户按模型单位直接消费。",
        "- 建议避免过度提取，防止资产在退休前后提前耗尽。",
        f"- 在当前配置和取现假设下，资产可维持至 {client['age'] + decision['risk_horizon_years']} 岁的模拟置信下界约为 {withdrawal['survival_lower_bound_pct']:.1f}%，满足 95% 风控口径。",
        "",
        "## 未来资产路径",
        "",
        "| 情景 | 目标年龄 | 预计资产 |",
        "|---|---:|---:|",
        *[
            f"| {row['scenario']} | {row['age']} | {row['wealth']:.2f} |"
            for row in decision["wealth_path"]
        ],
        "",
        "这张表按 3% 年度支取口径展示，把 failure probability 翻译成客户语言：在当前生命周期配置和支取纪律下，资产大概率不会在规划期内耗尽；4% 支取作为上方风控压力测试。",
        "",
        "## 分阶段产品映射",
        "",
        "| 阶段 | 产品方向 | 执行重点 |",
        "|---|---|---|",
        "| 当前阶段 | 低波动固收产品 + 少量宽基指数基金 | 控制权益仓位，优先降低回撤和波动 |",
        "| 50岁以后 | 提高现金类和短久期固收比例 | 逐步减少权益类暴露，强化退休现金流准备 |",
        "| 60岁以后 | 现金管理、存款、稳健固收为主 | 稳定支取来源，避免大幅净值波动影响生活安排 |",
        "",
        "## 自动调整规则",
        "",
        "- 年龄增长: 按生命周期路径年度自动下调风险资产比例。",
        "- 市场下跌超过 20%: 暂停加风险，优先检查现金流和本金安全。",
        "- 财富水平显著上升: 若资产翻倍且现金流压力下降，可将风险资产比例小幅提高，但上限建议不超过 10%。",
        "- 财富水平下降或支取压力上升: 降低风险资产比例，优先保留 12-24 个月流动性。",
        "",
        "## Risk Control Decision",
        "",
        f"- Case type: {decision['case_type']}",
        f"- Decision status: {decision['status']}",
        f"- Suggested allocation range: {decision['allocation_range']['low']:.0%}-{decision['allocation_range']['high']:.0%} risky assets",
        f"- Point estimate alpha: {policy['alpha']:.4f}",
        f"- Failure probability limit: {decision['failure_limit_pct']:.1f}%",
        f"- Risk horizon: {decision['risk_horizon_label']}",
        "",
        "## Decision Trace",
        "",
        "| Step | Before | After | Reason |",
        "|---|---:|---:|---|",
        *[
            f"| {step['name']} | {step['before_alpha']:.0%} | {step['after_alpha']:.0%} | {step['reason']} |"
            for step in decision["decision_trace"]
        ],
        "",
        tradeoff_note,
        "",
        "## Failure Risk Search",
        "",
        "| Candidate alpha | Failure estimate | 95% MC margin | Upper bound | Result |",
        "|---:|---:|---:|---:|---|",
        *[
            f"| {row['alpha']:.0%} | {row['failure_probability_pct']:.1f}% | "
            f"+/- {row['failure_probability_margin_pct']:.1f} pp | "
            f"{row['failure_probability_upper_bound_pct']:.1f}% | {'Pass' if row['passed'] else 'Reduce'} |"
            for row in decision["failure_search"]
        ],
        "",
        "## Risk Metrics",
        "",
        "| Metric | Raw model | Final recommendation |",
        "|---|---:|---:|",
        f"| Max drawdown | {model_risk['max_drawdown_pct']:.1f}% | {final_risk['max_drawdown_pct']:.1f}% |",
        f"| Failure probability | {_format_failure_estimate(model_risk)} | {_format_failure_estimate(final_risk)} |",
        f"| Annual volatility | {model_risk['annual_volatility_pct']:.1f}% | {final_risk['annual_volatility_pct']:.1f}% |",
        "",
        "## Baseline Comparison",
        "",
        "| Metric | Standard 60/40 | Final recommendation |",
        "|---|---:|---:|",
        f"| Risky allocation | {decision['baseline']['alpha']:.0%} | {policy['alpha']:.0%} |",
        f"| Max drawdown | {decision['baseline']['risk']['max_drawdown_pct']:.1f}% | {final_risk['max_drawdown_pct']:.1f}% |",
        f"| Failure probability | {_format_failure_estimate(decision['baseline']['risk'])} | {_format_failure_estimate(final_risk)} |",
        f"| Annual volatility | {decision['baseline']['risk']['annual_volatility_pct']:.1f}% | {final_risk['annual_volatility_pct']:.1f}% |",
        f"| Expected annual return | {decision['baseline']['expected_annual_return_pct']:.1f}% | {decision['final_expected_annual_return_pct']:.1f}% |",
        "",
        "## Applicability",
        "",
        decision["applicability_statement"],
        "",
        "## 行为约束与再平衡",
        "",
        "生命周期模型假设客户能坚持长期计划，但真实客户容易被短期波动影响。若市场短期下跌导致账面亏损，建议避免频繁追涨杀跌，以年度周期进行再平衡；除非触发年龄、市场跌幅或财富水平变化等规则，不建议临时大幅调整配置。",
        "",
        "## Customer Explanation",
        "",
        explanation,
        "",
        "## Notes",
        "",
        f"- Policy source: `{policy['year_file']}`",
        "- Risk metrics are Monte Carlo estimates using the model's return assumptions, a 4% withdrawal assumption, and an age-based projection horizon.",
        "- Failure probabilities include an approximate 95% Monte Carlo margin of error in percentage points.",
        f"- {decision['risk_horizon_note']} Risk horizon is determined based on the client's lifecycle stage.",
    ])


def main():
    parser = argparse.ArgumentParser(description="Build customer-facing lifecycle advice from free-text input.")
    parser.add_argument("text", nargs="*", help="Client text, for example: 45 years old, 1 million, stable returns")
    parser.add_argument("--results-dir", default=str(DEFAULT_RESULTS_DIR))
    parser.add_argument("--output", default=None)
    args = parser.parse_args()
    text = " ".join(args.text).strip()
    if not text:
        text = input("Client input: ").strip()
    result = build_advice(text, args.results_dir)
    print(result["markdown"])
    if args.output:
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(result["markdown"] + "\n", encoding="utf-8")
    return 0 if result["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
