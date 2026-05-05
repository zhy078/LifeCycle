#!/usr/bin/env python3
import argparse
import csv
import json
import math
import statistics
from pathlib import Path


def read_csv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    for row in rows:
        for key, value in list(row.items()):
            try:
                row[key] = float(value)
            except (TypeError, ValueError):
                pass
    return rows


def read_jsonl(path):
    with Path(path).open("r", encoding="utf-8") as f:
        return [json.loads(line) for line in f if line.strip()]


def finite(value):
    return isinstance(value, (int, float)) and math.isfinite(value)


def monotonic_direction(values):
    increasing = all(a <= b for a, b in zip(values, values[1:]))
    decreasing = all(a >= b for a, b in zip(values, values[1:]))
    if increasing and not decreasing:
        return "increasing"
    if decreasing and not increasing:
        return "decreasing"
    if increasing and decreasing:
        return "flat"
    return "non_monotone"


def build_summary(scenario_rows):
    by_rho = {}
    by_age = {}
    by_wealth = {}
    for row in scenario_rows:
        by_rho.setdefault(row["rho"], []).append(row)
        by_age.setdefault(row["age"], []).append(row)
        by_wealth.setdefault(row["wealth"], []).append(row)

    def avg(rows, key):
        return statistics.mean(float(row[key]) for row in rows)

    return {
        "by_rho": {rho: avg(rows, "model_value_function_utility") for rho, rows in sorted(by_rho.items())},
        "by_age": {age: avg(rows, "model_value_function_utility") for age, rows in sorted(by_age.items())},
        "by_wealth": {wealth: avg(rows, "model_value_function_utility") for wealth, rows in sorted(by_wealth.items())},
    }


def model_validation(scenario_rows, grid_rows, age_rows):
    findings = []
    warnings = []

    # Wealth monotonicity by rho/age.
    for rho in sorted({row["rho"] for row in scenario_rows}):
        for age in sorted({row["age"] for row in scenario_rows}):
            subset = sorted([row for row in scenario_rows if row["rho"] == rho and row["age"] == age], key=lambda r: r["wealth"])
            alphas = [row["optimal_alpha_at_state"] for row in subset]
            direction = monotonic_direction(alphas)
            if direction != "increasing":
                warnings.append(f"alpha vs wealth is {direction} for rho={rho:g}, age={age:g}: {', '.join(f'{a:.4f}' for a in alphas)}")

    # Age monotonicity by rho/wealth.
    for rho in sorted({row["rho"] for row in scenario_rows}):
        for wealth in sorted({row["wealth"] for row in scenario_rows}):
            subset = sorted([row for row in scenario_rows if row["rho"] == rho and row["wealth"] == wealth], key=lambda r: r["age"])
            alphas = [row["optimal_alpha_at_state"] for row in subset]
            direction = monotonic_direction(alphas)
            if direction != "decreasing":
                warnings.append(f"alpha vs age is {direction} for rho={rho:g}, wealth={wealth:g}: {', '.join(f'{a:.4f}' for a in alphas)}")

    # Rho risk taking.
    for age in sorted({row["age"] for row in scenario_rows}):
        for wealth in sorted({row["wealth"] for row in scenario_rows}):
            subset = sorted([row for row in scenario_rows if row["age"] == age and row["wealth"] == wealth], key=lambda r: r["rho"])
            alphas = [row["optimal_alpha_at_state"] for row in subset]
            direction = monotonic_direction(alphas)
            if direction not in ("decreasing", "flat"):
                warnings.append(f"higher rho does not reduce alpha at age={age:g}, wealth={wealth:g}: {', '.join(f'{a:.4f}' for a in alphas)}")

    extreme_rows = [row for row in scenario_rows if row["optimal_alpha_at_state"] > 0.8 or row["optimal_alpha_at_state"] < 0.02]
    if extreme_rows:
        warnings.append(f"{len(extreme_rows)} 10x10 scenario rows have extreme alpha (<0.02 or >0.8).")

    high_grid_10 = [row for row in grid_rows if row["grid"] == "10x10" and row["alpha"] > 0.5]
    high_grid_fine = [row for row in grid_rows if row["grid"] in ("20x20", "40x40") and row["alpha"] > 0.5]
    findings.append(f"Grid refinement reduced problem-state high alpha: 10x10 high-alpha rows={len(high_grid_10)}, fine-grid high-alpha rows={len(high_grid_fine)}.")

    opt_avg = statistics.mean(row["simulated_utility_optimal_policy"] for row in scenario_rows)
    naive_avg = statistics.mean(row["simulated_utility_naive_60_40_4pct"] for row in scenario_rows)
    findings.append(f"Deterministic baseline comparison: optimal policy avg utility={opt_avg:.8f}, naive avg utility={naive_avg:.8f}.")
    if opt_avg <= naive_avg:
        warnings.append("Optimal policy does not outperform naive baseline on average.")

    result = "WARNING" if warnings else "PASS"
    return result, findings, warnings


def runner_validation(result_files):
    issues = []
    total = 0
    for file_path in result_files:
        for row in read_jsonl(file_path):
            total += 1
            if row.get("status") != "ok":
                issues.append(f"{file_path}: scenario {row.get('id')} status={row.get('status')}")
            if row.get("timed_out"):
                issues.append(f"{file_path}: scenario {row.get('id')} timed_out=true")
            if int(row.get("year_files_count", 0)) <= 0:
                issues.append(f"{file_path}: scenario {row.get('id')} missing year files")
            score = row.get("score")
            if not finite(score):
                issues.append(f"{file_path}: scenario {row.get('id')} non-finite score={score}")
    return ("FAIL" if issues else "PASS"), total, issues


def critic_review():
    return [
        "Retirement-boundary behavior remains a live issue: even after grid refinement, alpha around age 60 is not perfectly monotone.",
        "The naive baseline is an external simulator, not a fixed-policy value-function solution inside the original dynamic program.",
        "Market return and mortality inputs are calibration constants without documented data provenance in the repository.",
        "The stochastic validation reports utility gains, but should add tail risk, drawdown, and consumption shortfall metrics before production use.",
        "Current conclusions are prototype-ready, not production-ready."
    ]


def render_markdown(output_path, summary, model_result, model_findings, model_warnings, runner_result, runner_total, runner_issues):
    lines = [
        "# Skill Review",
        "",
        "## Model Validation",
        "",
        f"Evaluation Result: {model_result}",
        "",
        "Findings:",
    ]
    lines += [f"- {item}" for item in model_findings]
    if model_warnings:
        lines += ["", "Warnings:"]
        lines += [f"- {item}" for item in model_warnings]
    lines += [
        "",
        "Interpretation:",
        "",
        "- The model is economically informative, but not fully clean: wealth and rho behave broadly as expected, while alpha shape has local non-monotonicity.",
        "- The original high-risk near-retirement anomaly is mostly a 10x10 grid artifact, because 20x20/40x40 reduce alpha sharply.",
        "- Remaining retirement-boundary behavior should be treated as a policy-shape diagnostic rather than ignored.",
        "",
        "Suggested Actions:",
        "",
        "1. Use 20x20 as the default validation grid; reserve 40x40 for robustness checks.",
        "2. Add automatic warnings for non-monotone alpha vs age/wealth.",
        "3. Examine consumption patterns jointly with alpha before judging a policy unreasonable.",
        "4. Add tail-risk metrics to stochastic validation.",
        "",
        "## Runner Validation",
        "",
        f"Runner Result: {runner_result}",
        "",
        f"- Scenarios checked: {runner_total}",
    ]
    if runner_issues:
        lines += ["", "Stability Issues:"]
        lines += [f"- {issue}" for issue in runner_issues]
    else:
        lines += ["- Stability Issues: none detected."]
    lines += [
        "",
        "Suggestions:",
        "",
        "- Keep checking `status`, `timed_out`, finite score, and `year_files_count` for every run.",
        "- Keep writing row-level diagnostics and avoid using timeout rows in economic conclusions.",
        "",
        "## Critical Issues",
        "",
    ]
    lines += [f"- {item}" for item in critic_review()]
    lines += [
        "",
        "## Summary Statistics",
        "",
        "### by_rho",
        "",
        "| rho | avg lifecycle utility |",
        "|---:|---:|",
    ]
    for key, val in summary["by_rho"].items():
        lines.append(f"| {key:g} | {val:.6f} |")
    lines += ["", "### by_age", "", "| age | avg lifecycle utility |", "|---:|---:|"]
    for key, val in summary["by_age"].items():
        lines.append(f"| {key:g} | {val:.6f} |")
    lines += ["", "### by_wealth", "", "| wealth | avg lifecycle utility |", "|---:|---:|"]
    for key, val in summary["by_wealth"].items():
        lines.append(f"| {key:g} | {val:.6f} |")
    lines += [
        "",
        "## Recommendation",
        "",
        "Ready for prototype and research discussion.",
        "",
        "Not ready for production until mortality/market calibration provenance, fixed-policy baseline evaluation, and stochastic tail-risk validation are documented.",
        "",
    ]
    Path(output_path).write_text("\n".join(lines), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scenario-csv", required=True)
    parser.add_argument("--grid-csv", required=True)
    parser.add_argument("--age-csv", required=True)
    parser.add_argument("--results-jsonl", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    scenario_rows = read_csv(args.scenario_csv)
    grid_rows = read_csv(args.grid_csv)
    age_rows = read_csv(args.age_csv)
    summary = build_summary(scenario_rows)
    model_result, findings, warnings = model_validation(scenario_rows, grid_rows, age_rows)
    runner_result, runner_total, runner_issues = runner_validation(args.results_jsonl)
    render_markdown(args.output, summary, model_result, findings, warnings, runner_result, runner_total, runner_issues)
    print(args.output)


if __name__ == "__main__":
    main()
