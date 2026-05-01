#!/usr/bin/env python3
import argparse
import itertools
import json
import math
import random
import shutil
import subprocess
import sys
import time
from pathlib import Path


def load_config(path: Path):
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def validate(cfg):
    if "constraints" not in cfg or "search_space" not in cfg:
        raise ValueError("config must include constraints and search_space")
    c = cfg["constraints"]
    required = ["rho", "delta", "psi", "mu", "sigr"]
    for k in required:
        if k not in cfg["search_space"] or not cfg["search_space"][k]:
            raise ValueError(f"missing search_space.{k}")

    def in_range(name, v, lo, hi):
        if not (lo <= v <= hi):
            raise ValueError(f"{name} out of range: {v}. expected [{lo}, {hi}]")

    for rho in cfg["search_space"]["rho"]:
        in_range("rho", rho, c.get("rho_min", 1.0), c.get("rho_max", 1000.0))
    for d in cfg["search_space"]["delta"]:
        in_range("delta", d, c.get("delta_min", 0.0), c.get("delta_max", 1.0))
    for psi in cfg["search_space"]["psi"]:
        in_range("psi", psi, c.get("psi_min", 0.0), c.get("psi_max", 1000.0))
    for mu in cfg["search_space"]["mu"]:
        in_range("mu", mu, c.get("mu_min", -1.0), c.get("mu_max", 1.0))
    for sigr in cfg["search_space"]["sigr"]:
        in_range("sigr", sigr, c.get("sigr_min", 0.0), c.get("sigr_max", 10.0))


def score_proxy(params, objective):
    rho = params["rho"]
    delta = params["delta"]
    psi = params["psi"]
    mu = params["mu"]
    sigr = params["sigr"]
    baseline = 10 * math.log1p(max(mu, -0.99) + 1.0) + 50 * delta - 0.7 * rho - 5 * sigr + 2 * psi
    if objective == "maximize_terminal_wealth":
        return baseline
    if objective == "maximize_risk_adjusted":
        return baseline - 3 * sigr
    if objective == "minimize_volatility_penalty":
        return -sigr + 0.2 * mu + 0.1 * psi
    return baseline


def octave_available():
    return shutil.which("octave") is not None


def run_model_stub(artifact_dir: Path, dry_run: bool):
    artifact_dir.mkdir(parents=True, exist_ok=True)
    if dry_run or not octave_available():
        return {"status": "dry_run", "artifact_dir": str(artifact_dir)}
    cmd = ["octave", "--quiet", "--eval", "disp('Hook your life_cycle wrapper here');"]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    (artifact_dir / "octave_stdout.txt").write_text(proc.stdout, encoding="utf-8")
    (artifact_dir / "octave_stderr.txt").write_text(proc.stderr, encoding="utf-8")
    return {
        "status": "ok" if proc.returncode == 0 else "error",
        "artifact_dir": str(artifact_dir),
        "code": proc.returncode,
    }


def build_candidates(cfg, method, max_evals, seed):
    keys = ["rho", "delta", "psi", "mu", "sigr"]
    spaces = [cfg["search_space"][k] for k in keys]

    if method == "grid":
        candidates = [dict(zip(keys, tup)) for tup in itertools.product(*spaces)]
        if max_evals and max_evals < len(candidates):
            random.Random(seed).shuffle(candidates)
            return candidates[:max_evals]
        return candidates

    # random method
    rng = random.Random(seed)
    all_points = [dict(zip(keys, tup)) for tup in itertools.product(*spaces)]
    if not max_evals:
        max_evals = min(200, len(all_points))
    return [rng.choice(all_points) for _ in range(max_evals)]


def parse_args():
    ap = argparse.ArgumentParser(description="LifeCycle optimizer runner")
    ap.add_argument("--config", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--search-method", choices=["grid", "random"], default=None)
    ap.add_argument("--max-evals", type=int, default=None)
    ap.add_argument("--progress-every", type=int, default=20)
    return ap.parse_args()


def main():
    args = parse_args()
    cfg = load_config(Path(args.config))
    validate(cfg)

    seed = cfg.get("seed", 0)
    random.seed(seed)
    method = args.search_method or cfg.get("search_method", "grid")
    max_evals = args.max_evals or cfg.get("max_evals")

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    objective = cfg.get("objective", "maximize_terminal_wealth")
    candidates = build_candidates(cfg, method, max_evals, seed)
    total = len(candidates)
    if total == 0:
        raise ValueError("no candidates generated")

    print(f"[optimizer] method={method} objective={objective} candidates={total} dry_run={args.dry_run}")
    sys.stdout.flush()

    start = time.time()
    results = []
    results_path = out / "results.jsonl"
    with results_path.open("w", encoding="utf-8") as f:
        for i, params in enumerate(candidates):
            scenario_dir = out / f"scenario_{i:04d}"
            run_meta = run_model_stub(scenario_dir, args.dry_run)
            score = score_proxy(params, objective)
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

            if (i + 1) % max(1, args.progress_every) == 0 or i == total - 1:
                elapsed = time.time() - start
                print(f"[optimizer] progress {i+1}/{total} elapsed={elapsed:.1f}s")
                sys.stdout.flush()

    results.sort(key=lambda r: r["score"], reverse=True)
    best = results[0]
    top_k = results[: cfg.get("top_k", 5)]
    (out / "best_params.json").write_text(json.dumps(best, indent=2, ensure_ascii=False), encoding="utf-8")

    report = {
        "objective": objective,
        "search_method": method,
        "total_scenarios": total,
        "elapsed_seconds": round(time.time() - start, 3),
        "best": best,
        "top_k": top_k,
    }
    (out / "report.md").write_text("# LifeCycle Optimization Report\n\n" + json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"[optimizer] done best_score={best['score']:.6f} output_dir={out}")


if __name__ == "__main__":
    main()
