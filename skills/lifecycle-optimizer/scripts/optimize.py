#!/usr/bin/env python3
import argparse
import itertools
import json
import math
import os
import random
import shutil
import subprocess
from pathlib import Path


def load_config(path: Path):
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def validate(cfg):
    c = cfg["constraints"]
    for k in ["rho", "delta", "psi", "mu", "sigr"]:
        if k not in cfg["search_space"] or not cfg["search_space"][k]:
            raise ValueError(f"missing search_space.{k}")
    def in_range(name, v, lo, hi):
        if not (lo <= v <= hi):
            raise ValueError(f"{name} out of range: {v}")
    for rho in cfg["search_space"]["rho"]:
        in_range("rho", rho, c["rho_min"], 1000)
    for d in cfg["search_space"]["delta"]:
        in_range("delta", d, c["delta_min"], c["delta_max"])
    for psi in cfg["search_space"]["psi"]:
        in_range("psi", psi, c["psi_min"], 1000)
    for mu in cfg["search_space"]["mu"]:
        in_range("mu", mu, c["mu_min"], c["mu_max"])
    for sigr in cfg["search_space"]["sigr"]:
        in_range("sigr", sigr, c["sigr_min"], c["sigr_max"])


def score_proxy(params, objective):
    # Placeholder objective proxy when full MATLAB pipeline is not wired.
    rho = params["rho"]
    delta = params["delta"]
    psi = params["psi"]
    mu = params["mu"]
    sigr = params["sigr"]
    util = 10 * math.log1p(max(mu, -0.99) + 1.0) + 50 * delta - 0.7 * rho - 5 * sigr + 2 * psi
    if objective == "maximize_terminal_wealth":
        return util
    if objective == "maximize_risk_adjusted":
        return util - 3 * sigr
    return util


def octave_available():
    return shutil.which("octave") is not None


def run_model_stub(params, artifact_dir: Path, dry_run: bool):
    artifact_dir.mkdir(parents=True, exist_ok=True)
    if dry_run or not octave_available():
        return {"status": "dry_run", "artifact_dir": str(artifact_dir)}
    cmd = ["octave", "--quiet", "--eval", "disp('Hook your life_cycle wrapper here');"]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    (artifact_dir / "octave_stdout.txt").write_text(proc.stdout, encoding="utf-8")
    (artifact_dir / "octave_stderr.txt").write_text(proc.stderr, encoding="utf-8")
    return {"status": "ok" if proc.returncode == 0 else "error", "artifact_dir": str(artifact_dir), "code": proc.returncode}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    cfg = load_config(Path(args.config))
    validate(cfg)
    random.seed(cfg.get("seed", 0))

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    keys = ["rho", "delta", "psi", "mu", "sigr"]
    grid = list(itertools.product(*[cfg["search_space"][k] for k in keys]))

    results_path = out / "results.jsonl"
    objective = cfg.get("objective", "maximize_terminal_wealth")
    results = []

    with results_path.open("w", encoding="utf-8") as f:
        for i, tup in enumerate(grid):
            params = dict(zip(keys, tup))
            score = score_proxy(params, objective)
            scenario_dir = out / f"scenario_{i:04d}"
            run_meta = run_model_stub(params, scenario_dir, args.dry_run)
            row = {
                "id": i,
                "params": params,
                "objective": objective,
                "score": score,
                "status": run_meta["status"],
                "artifact_dir": run_meta["artifact_dir"],
            }
            f.write(json.dumps(row, ensure_ascii=False) + "\n")
            results.append(row)

    results.sort(key=lambda r: r["score"], reverse=True)
    best = results[0]
    top_k = results[: cfg.get("top_k", 5)]

    (out / "best_params.json").write_text(json.dumps(best, indent=2, ensure_ascii=False), encoding="utf-8")
    report = [
        "# LifeCycle Optimization Report",
        f"Objective: {objective}",
        f"Total scenarios: {len(results)}",
        "",
        "## Best",
        json.dumps(best, indent=2, ensure_ascii=False),
        "",
        "## Top-K",
        json.dumps(top_k, indent=2, ensure_ascii=False),
    ]
    (out / "report.md").write_text("\n".join(report), encoding="utf-8")
    print(f"Done. Best score={best['score']:.6f}. Outputs: {out}")


if __name__ == "__main__":
    main()
