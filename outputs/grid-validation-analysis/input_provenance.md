# Input Provenance

This note records where the LifeCycle experiment inputs come from.

## 1. Client/scenario inputs

These were supplied by the experiment design, not by `life_cycle.m`:

| input | value used | source |
|---|---:|---|
| `age` | `[30, 45, 60]` | user experiment request |
| `wealth` | `[5, 10, 20]` | user experiment request |
| `rho` | `[6, 8, 10]` | user experiment request; interpreted as agent risk aversion / risk preference |

These values are used to query the model policy/value files after the Octave model solves.

## 2. Optimizer config inputs

These come from JSON case files under `outputs/` or `skills/lifecycle-optimizer/assets/`:

| input | example value | source |
|---|---:|---|
| `objective` | `maximize_lifetime_utility` | case JSON / skill convention |
| `delta` | `0.97` | case JSON |
| `psi` | `0.6` | case JSON |
| `mu` | `0.03` | case JSON, patches `life_cycle.m` |
| `sigr` | `0.15` | case JSON, patches `life_cycle.m` |
| `r` | `1.015` | case JSON, patches `life_cycle.m` |
| `tb` | `20` | case JSON, patches `life_cycle.m` |
| `tr` | `65` or `66` depending case | case JSON, patches `life_cycle.m` |
| `td` | `100` | case JSON, patches `life_cycle.m` |
| `na`, `ncash` | `10`, `20`, `40` | case JSON for grid tests; patches `life_cycle.m` in fast mode |

The runner copies `life_cycle.m` into each scenario artifact directory and patches these parameter lines before running Octave.

## 3. Market return inputs

The model uses:

| input | meaning | source |
|---|---|---|
| `r` | risk-free gross return | case JSON override; original hardcoded default in `life_cycle.m` is `1.015` |
| `mu` | risky excess return mean | case JSON override; original hardcoded default in `life_cycle.m` is `0.04` |
| `sigr` | risky return volatility | case JSON override; original hardcoded default in `life_cycle.m` is `0.2` |

In `life_cycle.m`, risky returns are represented by:

```matlab
gret(i1,1) = r + mu + grid(i1,1) * sigr;
```

In the stochastic validation script, I used the same interpretation:

```text
R_t = r + mu + sigr * eps_t
eps_t ~ N(0,1)
```

Important: these are assumptions in the case/config, not downloaded market data.

## 4. Cash / wealth grid inputs

The model's state grid is generated inside `life_cycle.m`:

| input | value | source |
|---|---:|---|
| `mincash` | `0.25` | hardcoded in `life_cycle.m` |
| `maxcash` | `200.0` | hardcoded in `life_cycle.m` |
| `ncash` | `10`, `20`, or `40` in experiments | case JSON override / fast-mode patch |

`life_cycle.m` builds a log-spaced cash grid:

```matlab
l_maxcash = log(maxcash);
l_mincash = log(mincash);
stepcash = (l_maxcash-l_mincash)/(ncash-1);
gcash(i1,1)=exp(lgcash(i1,1));
```

So the `wealth = [5, 10, 20]` experiment points are query points interpolated on this model cash grid.

## 5. Survival probability inputs

The survival probabilities are hardcoded directly in `life_cycle.m`:

```matlab
survprob(1,1)  = 0.99845;
...
survprob(80,1) = 0.6809;
```

They are used in the Bellman recursion, for example:

```matlab
auxVV = auxVV + weig(i8,1) * survprob(t,1) * (int_V^(1.0-rho));
```

Current limitation:

- I found the values in the model file, but the repository does not currently cite the mortality table or paper/source these survival probabilities came from.
- Treat them as model calibration constants unless we add a citation or replace them with a documented mortality table.

## 6. Labor income and retirement income inputs

Labor income profile parameters are hardcoded in `life_cycle.m`:

| input | value | source |
|---|---:|---|
| `aa` | `-2.170042 + 2.700381` | hardcoded |
| `b1` | `0.16818` | hardcoded |
| `b2` | `-0.0323371 / 10` | hardcoded |
| `b3` | `0.0019704 / 100` | hardcoded |
| `ret_fac` | `0.68212` | hardcoded |
| `smay` | `0.1` | hardcoded |
| `smav` | `0.1` | hardcoded |
| `corr_y`, `corr_v` | `0.0` | hardcoded |

The income profile is generated in `life_cycle.m` as:

```matlab
f_y(i1-tb+1,1) = exp(aa + b1*i1 + b2*i1^2 + b3*i1^3);
```

Retirement income is represented with:

```matlab
simY(t,:) = ret_fac;
```

Again, these are calibration constants in the source file, not external inputs read by the runner.

## 7. Quadrature / shock approximation

The model hardcodes a 5-point approximation to the normal distribution:

```matlab
grid(1,1) = -2.85697001387280;
...
weig(1,1) = 0.01125741132772;
...
```

These are used for expectation calculations in the dynamic program.

## 8. What should be documented next

The repository should ideally add citations or data sources for:

- survival probabilities;
- labor income profile coefficients;
- retirement income replacement factor `ret_fac`;
- return assumptions `r`, `mu`, `sigr`;
- cash grid bounds `mincash`, `maxcash`.

Right now, the exact technical source is `life_cycle.m`; the economic/data provenance is not documented in the repo.
