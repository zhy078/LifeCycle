def explain_policy(age, wealth, rho, alpha, lifecycle_utility=None, risk_metrics=None):
    risky_pct = round(alpha * 100)
    stable_pct = max(0, 100 - risky_pct)
    lines = [
        f"Recommended portfolio: {risky_pct}% risky assets, {stable_pct}% stabilizing assets."
    ]
    if lifecycle_utility is not None:
        lines.append(f"Lifecycle utility: {lifecycle_utility:.4f}.")
    if risk_metrics:
        lines.extend([
            "",
            "Risk metrics:",
            f"- Max drawdown: {risk_metrics['max_drawdown_pct']:.1f}%",
            f"- Failure probability: {risk_metrics['failure_probability_pct']:.1f}%",
            f"- Annual volatility: {risk_metrics['annual_volatility_pct']:.1f}%",
        ])

    reasons = []
    if age >= 60:
        reasons.append("Because you are close to retirement, the recommendation should protect retirement cash-flow stability and avoid large late-stage drawdowns.")
    elif age >= 40:
        reasons.append("Because you are in the wealth accumulation stage, the model can keep some risk exposure to improve lifecycle utility while controlling retirement risk.")
    else:
        reasons.append("Because retirement is still far away, the model can tolerate more market risk when it improves long-term utility.")

    if wealth <= 5:
        reasons.append("Your current wealth is low in the model grid, so results may be more sensitive to borrowing constraints, labor-income assumptions, and grid resolution.")
    elif wealth >= 20:
        reasons.append("Your wealth buffer is relatively high, so the model may allow more flexibility, but age and risk preference still matter.")
    else:
        reasons.append("Your wealth level is in the middle of the tested range, so the result is less likely to be an edge-grid artifact.")

    if rho >= 10:
        reasons.append("Your mapped rho is conservative, which means the model should penalize risky outcomes more heavily.")
    elif rho <= 6:
        reasons.append("Your mapped rho is aggressive, which means the model accepts more volatility when expected utility improves.")
    else:
        reasons.append("Your mapped rho is balanced, so the policy should trade off growth and stability.")

    cautions = []
    if wealth <= 5 and alpha > 0.5:
        cautions.append("Low wealth with high risky allocation is a suspicious case; check borrowing constraints, implicit labor income, and grid resolution.")
    if age >= 60 and alpha > 0.5:
        cautions.append("Near-retirement high risky allocation may violate economic intuition; use at least a 20x20 grid and inspect retirement-boundary behavior.")

    lines.extend(["", "Explanation:", " ".join(reasons)])
    if cautions:
        lines.extend(["", "Caution:", " ".join(cautions)])
    return "\n".join(lines)
