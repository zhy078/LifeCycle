import math
import random


def estimate_risk_metrics(alpha, wealth, mu=0.03, sigma=0.15, risk_free=0.015, consumption_rate=0.04, years=30, paths=1000, seed=7):
    rng = random.Random(seed)
    terminal_wealth = []
    drawdowns = []
    annual_returns = []
    failures = 0

    for _ in range(paths):
        w = max(float(wealth), 1e-9)
        nav = 1.0
        peak = nav
        max_drawdown = 0.0
        failed = False

        for _year in range(years):
            risky_return = rng.gauss(mu, sigma)
            portfolio_return = alpha * risky_return + (1 - alpha) * risk_free
            annual_returns.append(portfolio_return)
            nav = nav * (1 + portfolio_return)
            peak = max(peak, nav)
            if peak > 0:
                max_drawdown = min(max_drawdown, max(-1.0, (nav - peak) / peak))

            w = w * (1 + portfolio_return) - max(float(wealth) * consumption_rate, 0.0)
            if w <= 0:
                failed = True
                w = 0.0
                break

        if failed:
            failures += 1
        terminal_wealth.append(w)
        drawdowns.append(max_drawdown)

    avg_return = sum(annual_returns) / len(annual_returns) if annual_returns else 0.0
    variance = sum((r - avg_return) ** 2 for r in annual_returns) / max(len(annual_returns) - 1, 1)
    volatility = math.sqrt(variance)
    failure_probability = failures / paths
    failure_margin_pct = 100.0 * 1.96 * math.sqrt(
        failure_probability * (1.0 - failure_probability) / max(paths, 1)
    )

    terminal_wealth_sorted = sorted(terminal_wealth)

    def percentile(values, q):
        if not values:
            return 0.0
        position = (len(values) - 1) * q
        lower = int(math.floor(position))
        upper = int(math.ceil(position))
        if lower == upper:
            return values[lower]
        weight = position - lower
        return values[lower] * (1.0 - weight) + values[upper] * weight

    return {
        "max_drawdown_pct": 100.0 * (sum(drawdowns) / len(drawdowns)),
        "failure_probability_pct": 100.0 * failure_probability,
        "failure_probability_margin_pct": failure_margin_pct,
        "annual_volatility_pct": 100.0 * volatility,
        "avg_terminal_wealth": sum(terminal_wealth) / len(terminal_wealth),
        "terminal_wealth_p10": percentile(terminal_wealth_sorted, 0.10),
        "terminal_wealth_p50": percentile(terminal_wealth_sorted, 0.50),
        "terminal_wealth_p90": percentile(terminal_wealth_sorted, 0.90),
        "paths": paths,
        "years": years,
    }
