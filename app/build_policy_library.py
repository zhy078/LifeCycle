#!/usr/bin/env python3
"""Build an offline lifecycle policy library from repeated life_cycle.m runs."""
import argparse
import json
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "skills" / "lifecycle-optimizer" / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from optimize import run_model


MARKET_PARAMS = {
    "base": {"mu": 0.04, "sigr": 0.20},
    "conservative": {"mu": 0.03, "sigr": 0.15},
}


def scenario_name(rho, income_profile, market_profile):
    return f"rho_{int(round(float(rho))):02d}_income_{income_profile}_market_{market_profile}"


def main() -> int:
    parser = argparse.ArgumentParser(description="Build offline policy folders for customer lookup.")
    parser.add_argument("--output-dir", default="outputs/policy_library")
    parser.add_argument("--rhos", nargs="+", type=float, default=[6.0, 8.0, 10.0])
    parser.add_argument("--income-profiles", nargs="+", default=["stable", "volatile", "retired"])
    parser.add_argument("--market-profiles", nargs="+", default=["base", "conservative"])
    parser.add_argument("--delta", type=float, default=0.97)
    parser.add_argument("--psi", type=float, default=0.5)
    parser.add_argument("--fast-mode", action="store_true")
    parser.add_argument("--timeout-sec", type=int, default=1200)
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    if not output_dir.is_absolute():
        output_dir = ROOT / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    fixed_base = {
        "tb": 20,
        "tr": 66,
        "td": 100,
        "nsim": 10000,
        "r": 1.015,
    }

    rows = []
    start = time.time()
    for rho in args.rhos:
        for income_profile in args.income_profiles:
            for market_profile in args.market_profiles:
                market = MARKET_PARAMS[market_profile]
                params = {
                    "rho": rho,
                    "delta": args.delta,
                    "psi": args.psi,
                    "mu": market["mu"],
                    "sigr": market["sigr"],
                }
                fixed = {
                    **fixed_base,
                    "income_profile": income_profile,
                    "market_profile": market_profile,
                }
                name = scenario_name(rho, income_profile, market_profile)
                scenario_dir = output_dir / name
                print(f"[policy-library] running {name}")
                run = run_model(
                    params=params,
                    fixed=fixed,
                    artifact_dir=scenario_dir,
                    dry_run=False,
                    use_real_model=True,
                    fast_mode=args.fast_mode,
                    timeout_sec=args.timeout_sec,
                    client_profile=None,
                )
                row = {
                    "scenario_name": name,
                    "params": params,
                    "fixed": fixed,
                    "artifact_dir": str(scenario_dir),
                    "status": run.get("status"),
                    "run_seconds": run.get("run_seconds"),
                    "year_files_count": run.get("year_files_count"),
                    "policy_summary_path": run.get("policy_summary_path"),
                    "diagnostic_path": run.get("diagnostic_path"),
                }
                rows.append(row)
                (output_dir / "library_manifest.json").write_text(
                    json.dumps({
                        "elapsed_seconds": round(time.time() - start, 3),
                        "scenarios": rows,
                    }, indent=2, ensure_ascii=False),
                    encoding="utf-8",
                )

    print(f"[policy-library] done scenarios={len(rows)} output_dir={output_dir}")
    return 0 if all(row["status"] == "ok" for row in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
