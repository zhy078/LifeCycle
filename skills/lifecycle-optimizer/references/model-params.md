# Model Parameters Reference

## Core model parameters in `life_cycle.m`

- `tb, tr, td`: lifecycle age boundaries.
- `rho, delta, psi`: preferences.
- `r, mu, sigr`: risk-free gross return, risky excess return mean, risky return volatility.
- `nsim`: Monte Carlo simulation count.

## Suggested validation

- `rho > 1`
- `0 < delta < 1`
- `psi > 0`
- `sigr > 0`
- `tb < tr < td`

## Optimization guidance

- Start with coarse grid, then narrow around top region.
- Keep `seed` fixed for reproducibility.
- Compare top-K scenarios, not just the best point.
