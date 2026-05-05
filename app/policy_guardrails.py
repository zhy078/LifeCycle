HARD_ALPHA_CAP = 0.60
FAILURE_REJECT_THRESHOLD_PCT = 5.0
ALPHA_STEP = 0.05


def alpha_cap_for_client(age, wealth):
    cap = HARD_ALPHA_CAP
    reasons = []

    if age >= 60:
        cap = min(cap, 0.50)
        reasons.append("Age >= 60 caps risky allocation at 50%.")
    elif age >= 50:
        cap = min(cap, 0.60)
        reasons.append("Age >= 50 caps risky allocation at 60%.")

    if wealth < 2:
        cap = min(cap, 0.40)
        reasons.append("Low wealth caps risky allocation at 40%.")

    return cap, reasons


def _safe_alpha(alpha, risk_estimator):
    candidate = max(0.0, min(1.0, float(alpha)))
    search_path = []
    while candidate > 0.0:
        risk = risk_estimator(candidate)
        upper_bound = risk["failure_probability_pct"] + risk["failure_probability_margin_pct"]
        passed = upper_bound <= FAILURE_REJECT_THRESHOLD_PCT
        search_path.append({
            "alpha": candidate,
            "failure_probability_pct": risk["failure_probability_pct"],
            "failure_probability_margin_pct": risk["failure_probability_margin_pct"],
            "failure_probability_upper_bound_pct": upper_bound,
            "passed": passed,
        })
        if passed:
            return candidate, risk, search_path
        candidate = round(max(0.0, candidate - ALPHA_STEP), 10)
    risk = risk_estimator(0.0)
    search_path.append({
        "alpha": 0.0,
        "failure_probability_pct": risk["failure_probability_pct"],
        "failure_probability_margin_pct": risk["failure_probability_margin_pct"],
        "failure_probability_upper_bound_pct": risk["failure_probability_pct"] + risk["failure_probability_margin_pct"],
        "passed": True,
    })
    return 0.0, risk, search_path


def _trace_step(name, before, after, reason):
    return {
        "name": name,
        "before_alpha": before,
        "after_alpha": after,
        "reason": reason,
    }


def apply_policy_guardrails(client, model_alpha, model_risk, risk_estimator):
    """Convert model alpha into a customer-safe decision with an auditable trace."""
    model_alpha = max(0.0, min(1.0, float(model_alpha)))
    alpha = model_alpha
    reasons = []
    trace = []
    failure_search = []

    if model_risk["failure_probability_pct"] > FAILURE_REJECT_THRESHOLD_PCT:
        before = alpha
        alpha, adjusted_risk, failure_search = _safe_alpha(alpha, risk_estimator)
        reasons.append(
            f"Failure probability {model_risk['failure_probability_pct']:.1f}% exceeds "
            f"the {FAILURE_REJECT_THRESHOLD_PCT:.1f}% limit."
        )
        trace.append(_trace_step(
            "failure_risk",
            before,
            alpha,
            f"Reduced until estimated failure probability remains within {FAILURE_REJECT_THRESHOLD_PCT:.1f}% under Monte Carlo uncertainty bounds.",
        ))
    else:
        adjusted_risk = model_risk

    before = alpha
    if client["age"] >= 60:
        alpha = min(alpha, 0.50)
        reasons.append("Retirement-stage client caps risky allocation at 50%.")
    elif client["age"] >= 50:
        alpha = min(alpha, 0.60)
        reasons.append("Pre-retirement client caps risky allocation at 60%.")
    if alpha != before or client["age"] >= 50:
        trace.append(_trace_step(
            "age_constraint",
            before,
            alpha,
            "Applied age-based risky-allocation cap.",
        ))

    before = alpha
    if client["wealth"] < 2:
        alpha = min(alpha, 0.40)
        reasons.append("Low wealth requires extra protection from large drawdowns.")
    if alpha != before or client["wealth"] < 2:
        trace.append(_trace_step(
            "wealth_constraint",
            before,
            alpha,
            "Applied low-wealth protection cap.",
        ))

    adjusted_risk = risk_estimator(alpha)

    status = "approved"
    case_type = "SAFE_MODEL_OUTPUT"
    if abs(alpha - model_alpha) > 1e-9:
        status = "adjusted"
        case_type = "UNSAFE_MODEL_OUTPUT"
    if model_risk["failure_probability_pct"] > FAILURE_REJECT_THRESHOLD_PCT:
        status = "adjusted"
        case_type = "UNSAFE_MODEL_OUTPUT"

    return {
        "status": status,
        "case_type": case_type,
        "model_alpha": model_alpha,
        "final_alpha": alpha,
        "reasons": reasons,
        "decision_trace": trace,
        "failure_search": failure_search,
        "model_risk": model_risk,
        "final_risk": adjusted_risk,
        "failure_limit_pct": FAILURE_REJECT_THRESHOLD_PCT,
    }
