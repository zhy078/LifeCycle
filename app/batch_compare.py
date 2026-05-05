import argparse
import csv
import math
import random
import sys
from pathlib import Path


APP_DIR = Path(__file__).resolve().parent
if str(APP_DIR) not in sys.path:
    sys.path.insert(0, str(APP_DIR))

from advisor import DEFAULT_RESULTS_DIR, _policy_at_client, _select_result
from client_input import RISK_MAP


ROOT = Path(__file__).resolve().parents[1]
RISK_LABELS = ["conservative", "balanced", "aggressive"]
NAIVE_ALPHA = 0.60


def crra_utility(x, rho):
    x = max(float(x), 1e-9)
    if abs(rho - 1.0) < 1e-9:
        return math.log(x)
    return (x ** (1.0 - rho)) / (1.0 - rho)


def generate_clients(n=100, seed=42):
    rng = random.Random(seed)
    clients = []
    for idx in range(1, n + 1):
        risk = rng.choice(RISK_LABELS)
        wealth = round(rng.uniform(5.0, 20.0), 2)
        age = rng.randint(30, 65)
        clients.append({
            "client_id": idx,
            "age": age,
            "wealth": wealth,
            "risk_preference": risk,
            "rho": RISK_MAP[risk],
            "raw_text": f"{age} years old, wealth={wealth:.2f}, {risk}",
        })
    return clients


def simulate_strategy(alpha, wealth, rho, mu, sigma, risk_free=0.015, consumption_rate=0.04, years=30, paths=1000, seed=1):
    rng = random.Random(seed)
    terminal = []
    drawdowns = []
    returns = []
    utilities = []
    failures = 0
    base_consumption = max(float(wealth) * consumption_rate, 1e-9)
    consumption_floor = base_consumption * 0.25
    beta = 0.97

    for _path in range(paths):
        w = max(float(wealth), 1e-9)
        nav = 1.0
        peak = 1.0
        max_drawdown = 0.0
        path_utility = 0.0
        failed = False

        for year in range(years):
            risky_return = rng.gauss(mu, sigma)
            portfolio_return = alpha * risky_return + (1.0 - alpha) * risk_free
            returns.append(portfolio_return)

            nav *= (1.0 + portfolio_return)
            peak = max(peak, nav)
            max_drawdown = min(max_drawdown, max(-1.0, (nav - peak) / peak))

            consumption = min(base_consumption, max(w, 0.0))
            path_utility += (beta ** year) * crra_utility(consumption, rho)
            w = w * (1.0 + portfolio_return) - base_consumption
            if w <= 0:
                failed = True
                w = 0.0
                for future in range(year + 1, years):
                    path_utility += (beta ** future) * crra_utility(consumption_floor, rho)
                break

        path_utility += (beta ** years) * crra_utility(max(w, consumption_floor), rho)
        failures += int(failed)
        terminal.append(w)
        drawdowns.append(max_drawdown)
        utilities.append(path_utility)

    avg_return = sum(returns) / len(returns)
    variance = sum((r - avg_return) ** 2 for r in returns) / max(len(returns) - 1, 1)
    return {
        "alpha": alpha,
        "avg_terminal_wealth": sum(terminal) / len(terminal),
        "failure_probability_pct": 100.0 * failures / paths,
        "max_drawdown_pct": 100.0 * sum(drawdowns) / len(drawdowns),
        "annual_volatility_pct": 100.0 * math.sqrt(variance),
        "simulated_crra_utility": sum(utilities) / len(utilities),
    }


def compare_clients(clients, results_dir=DEFAULT_RESULTS_DIR, paths=1000, years=30, seed=2026):
    rows = []
    for client in clients:
        model_row = _select_result(Path(results_dir), client["rho"])
        if not model_row:
            raise RuntimeError(f"No completed model result for rho={client['rho']}")
        policy = _policy_at_client(model_row, client)
        if not policy or policy.get("alpha") is None:
            raise RuntimeError(f"No policy for client {client['client_id']}")

        alpha = max(0.0, min(1.0, float(policy["alpha"])))
        mu = float(model_row.get("params", {}).get("mu", 0.03))
        sigma = float(model_row.get("params", {}).get("sigr", 0.15))
        shared_seed = seed + client["client_id"] * 1009
        lifecycle = simulate_strategy(alpha, client["wealth"], client["rho"], mu, sigma, years=years, paths=paths, seed=shared_seed)
        naive = simulate_strategy(NAIVE_ALPHA, client["wealth"], client["rho"], mu, sigma, years=years, paths=paths, seed=shared_seed)

        rows.append({
            **client,
            "model_alpha": alpha,
            "naive_alpha": NAIVE_ALPHA,
            "model_lifecycle_utility": policy.get("utility"),
            "model_simulated_crra_utility": lifecycle["simulated_crra_utility"],
            "naive_simulated_crra_utility": naive["simulated_crra_utility"],
            "utility_lift": lifecycle["simulated_crra_utility"] - naive["simulated_crra_utility"],
            "model_failure_probability_pct": lifecycle["failure_probability_pct"],
            "naive_failure_probability_pct": naive["failure_probability_pct"],
            "failure_probability_lift_pct": lifecycle["failure_probability_pct"] - naive["failure_probability_pct"],
            "model_max_drawdown_pct": lifecycle["max_drawdown_pct"],
            "naive_max_drawdown_pct": naive["max_drawdown_pct"],
            "max_drawdown_lift_pct": lifecycle["max_drawdown_pct"] - naive["max_drawdown_pct"],
            "model_volatility_pct": lifecycle["annual_volatility_pct"],
            "naive_volatility_pct": naive["annual_volatility_pct"],
            "volatility_lift_pct": lifecycle["annual_volatility_pct"] - naive["annual_volatility_pct"],
            "model_avg_terminal_wealth": lifecycle["avg_terminal_wealth"],
            "naive_avg_terminal_wealth": naive["avg_terminal_wealth"],
            "terminal_wealth_lift": lifecycle["avg_terminal_wealth"] - naive["avg_terminal_wealth"],
        })
    return rows


def avg(rows, key):
    return sum(float(row[key]) for row in rows) / len(rows)


def pct(count, total):
    return 100.0 * count / total if total else 0.0


def render_summary(rows, output_csv, paths, years):
    utility_wins = sum(1 for row in rows if row["utility_lift"] > 0)
    failure_wins = sum(1 for row in rows if row["failure_probability_lift_pct"] < 0)
    drawdown_wins = sum(1 for row in rows if row["max_drawdown_lift_pct"] > 0)
    terminal_wins = sum(1 for row in rows if row["terminal_wealth_lift"] > 0)

    lines = [
        "# 100-Client Lifecycle vs Naive 60/40 Comparison",
        "",
        "## Setup",
        "",
        f"- Clients: {len(rows)} random clients",
        f"- Simulation: {paths} Monte Carlo paths per client, {years} years",
        "- Lifecycle policy: alpha interpolated from completed 20x20 lifecycle skill artifacts",
        "- Naive trading baseline: fixed 60% risky / 40% stabilizing assets",
        "- Consumption assumption: fixed 4% of initial wealth per year",
        "- Utility metric: simulated CRRA utility using each client's mapped rho and a 25% minimum-consumption floor after depletion; this is a product-layer comparison metric, not the model value function itself.",
        f"- Detail CSV: `{output_csv}`",
        "",
        "## Aggregate Result",
        "",
        "| Metric | Lifecycle | Naive 60/40 | Lift (Lifecycle - Naive) |",
        "|---|---:|---:|---:|",
        f"| Utility win rate | {pct(utility_wins, len(rows)):.1f}% | {(100.0 - pct(utility_wins, len(rows))):.1f}% | {utility_wins}/{len(rows)} clients |",
        f"| Failure probability | {avg(rows, 'model_failure_probability_pct'):.2f}% | {avg(rows, 'naive_failure_probability_pct'):.2f}% | {avg(rows, 'failure_probability_lift_pct'):.2f} pp |",
        f"| Max drawdown | {avg(rows, 'model_max_drawdown_pct'):.2f}% | {avg(rows, 'naive_max_drawdown_pct'):.2f}% | {avg(rows, 'max_drawdown_lift_pct'):.2f} pp |",
        f"| Annual volatility | {avg(rows, 'model_volatility_pct'):.2f}% | {avg(rows, 'naive_volatility_pct'):.2f}% | {avg(rows, 'volatility_lift_pct'):.2f} pp |",
        f"| Terminal wealth | {avg(rows, 'model_avg_terminal_wealth'):.4f} | {avg(rows, 'naive_avg_terminal_wealth'):.4f} | {avg(rows, 'terminal_wealth_lift'):.4f} |",
        "",
        "## Win Counts",
        "",
        f"- Utility higher: {utility_wins}/{len(rows)} clients",
        f"- Failure probability lower: {failure_wins}/{len(rows)} clients",
        f"- Drawdown less severe: {drawdown_wins}/{len(rows)} clients",
        f"- Terminal wealth higher: {terminal_wins}/{len(rows)} clients",
        "",
        "## By Risk Preference",
        "",
        "| Risk preference | Clients | Avg lifecycle alpha | Utility win rate | Failure lift | Volatility lift |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for risk in RISK_LABELS:
        subset = [row for row in rows if row["risk_preference"] == risk]
        if subset:
            subset_utility_wins = sum(1 for row in subset if row["utility_lift"] > 0)
            lines.append(
                f"| {risk} | {len(subset)} | {avg(subset, 'model_alpha'):.3f} | "
                f"{pct(subset_utility_wins, len(subset)):.1f}% | {avg(subset, 'failure_probability_lift_pct'):.2f} pp | "
                f"{avg(subset, 'volatility_lift_pct'):.2f} pp |"
            )

    lines.extend([
        "",
        "## Interpretation",
        "",
        "- Positive utility lift means the lifecycle policy improved simulated client utility versus naive 60/40 under the same return shocks.",
        "- Negative volatility lift or positive max-drawdown lift means the lifecycle policy reduced risk versus naive 60/40.",
        "- Terminal wealth is reported separately because lifecycle utility is not terminal-wealth maximization.",
        "- This is a product-layer stochastic validation, not a replacement for solving a fixed-policy naive baseline inside the dynamic model.",
    ])
    return "\n".join(lines) + "\n"


def write_csv(rows, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description="Compare lifecycle policy against naive 60/40 for random clients.")
    parser.add_argument("--clients", type=int, default=100)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--paths", type=int, default=1000)
    parser.add_argument("--years", type=int, default=30)
    parser.add_argument("--results-dir", default=str(DEFAULT_RESULTS_DIR))
    parser.add_argument("--output-dir", default="outputs/client-batch-comparison")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    if not output_dir.is_absolute():
        output_dir = ROOT / output_dir
    clients = generate_clients(args.clients, args.seed)
    rows = compare_clients(clients, args.results_dir, args.paths, args.years, args.seed)
    csv_path = output_dir / "client_100_comparison.csv"
    summary_path = output_dir / "client_100_comparison_summary.md"
    write_csv(rows, csv_path)
    summary = render_summary(rows, csv_path, args.paths, args.years)
    summary_path.write_text(summary, encoding="utf-8")
    print(summary)
    print(f"Saved CSV: {csv_path}")
    print(f"Saved summary: {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
