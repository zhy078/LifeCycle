#!/usr/bin/env python3
import argparse
import itertools
import json
import math
import random
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path


PARAM_LINE_RE = re.compile(r"^(?P<name>\w+)\s*=\s*[^;]+;\s*$")


def resolve_output_dir(raw_path: str) -> Path:
    p = Path(raw_path)
    if p.is_absolute():
        return p
    repo_root = Path(__file__).resolve().parents[3]
    return (repo_root / p).resolve()


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


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


def patch_lifecycle_script(src: Path, dst: Path, params: dict, fixed: dict, fast_mode: bool):
    replacements = {
        "rho": params["rho"],
        "delta": params["delta"],
        "psi": params["psi"],
        "mu": params["mu"],
        "sigr": params["sigr"],
        "r": fixed.get("r", 1.015),
        "tb": fixed.get("tb", 20),
        "tr": fixed.get("tr", 66),
        "td": fixed.get("td", 100),
        "nsim": fixed.get("nsim", 10000),
    }
    if fast_mode:
        replacements.update({
            "na": fixed.get("na", 21),
            "ncash": fixed.get("ncash", 21),
            "n": fixed.get("n", 3),
            "nsim": fixed.get("nsim", 1000),
        })

    lines = src.read_text(encoding="utf-8").splitlines()
    out_lines = []
    for ln in lines:
        m = PARAM_LINE_RE.match(ln.strip())
        if m and m.group("name") in replacements:
            name = m.group("name")
            val = replacements[name]
            if isinstance(val, str):
                out_lines.append(f"{name} = {val};")
            else:
                out_lines.append(f"{name} = {val};")
        else:
            out_lines.append(ln)
    dst.write_text("\n".join(out_lines) + "\n", encoding="utf-8")


def parse_terminal_metric(run_dir: Path):
    year80 = run_dir / "year80.txt"
    if not year80.exists():
        return None
    nums = []
    for token in year80.read_text(encoding="utf-8", errors="ignore").replace(",", " ").split():
        try:
            nums.append(float(token))
        except ValueError:
            pass
    if not nums:
        return None
    return nums[-1]


def run_real_model(params, fixed, artifact_dir: Path, fast_mode: bool, timeout_sec: int):
    artifact_dir.mkdir(parents=True, exist_ok=True)
    root = repo_root()
    needed = ["f_spline.m", "f_sc_splint.m", "f_ntoil.m", "f_randn.m"]
    for fn in needed:
        shutil.copy2(root / fn, artifact_dir / fn)

    patch_lifecycle_script(root / "life_cycle.m", artifact_dir / "life_cycle.m", params, fixed, fast_mode)

    cmd = ["octave", "--quiet", "life_cycle.m"]
    proc = subprocess.run(cmd, cwd=str(artifact_dir), capture_output=True, text=True, timeout=timeout_sec)
    (artifact_dir / "octave_stdout.txt").write_text(proc.stdout, encoding="utf-8")
    (artifact_dir / "octave_stderr.txt").write_text(proc.stderr, encoding="utf-8")
    metric = parse_terminal_metric(artifact_dir)
    return {
        "status": "ok" if proc.returncode == 0 else "error",
        "artifact_dir": str(artifact_dir),
        "code": proc.returncode,
        "metric": metric,
    }


def run_model(params, fixed, artifact_dir: Path, dry_run: bool, use_real_model: bool, fast_mode: bool, timeout_sec: int):
    artifact_dir.mkdir(parents=True, exist_ok=True)
    if dry_run:
        return {"status": "dry_run", "artifact_dir": str(artifact_dir), "metric": None}
    if use_real_model:
        if not octave_available():
            return {"status": "error_no_octave", "artifact_dir": str(artifact_dir), "metric": None}
        try:
            return run_real_model(params, fixed, artifact_dir, fast_mode, timeout_sec)
        except subprocess.TimeoutExpired:
            return {"status": "error_timeout", "artifact_dir": str(artifact_dir), "metric": None}
    return {"status": "simulated", "artifact_dir": str(artifact_dir), "metric": None}


def build_candidates(cfg, method, max_evals, seed):
    keys = ["rho", "delta", "psi", "mu", "sigr"]
    spaces = [cfg["search_space"][k] for k in keys]
    all_points = [dict(zip(keys, tup)) for tup in itertools.product(*spaces)]

    if method == "grid":
        if max_evals and max_evals < len(all_points):
            random.Random(seed).shuffle(all_points)
            return all_points[:max_evals]
        return all_points

    rng = random.Random(seed)
    if not max_evals:
        max_evals = min(200, len(all_points))
    return [rng.choice(all_points) for _ in range(max_evals)]


def parse_args():
    ap = argparse.ArgumentParser(description="LifeCycle optimizer runner")
    ap.add_argument("--config", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--use-real-model", action="store_true", help="run actual life_cycle.m via octave")
    ap.add_argument("--allow-proxy-fallback", action="store_true", help="when real model fails, keep proxy score instead of exiting")
    ap.add_argument("--fast-mode", action="store_true", help="reduce grid/simulation size for quicker runtime")
    ap.add_argument("--timeout-sec", type=int, default=120)
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
    fixed = cfg.get("fixed", {})

    out = resolve_output_dir(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    objective = cfg.get("objective", "maximize_terminal_wealth")
    candidates = build_candidates(cfg, method, max_evals, seed)
    total = len(candidates)
    if total == 0:
        raise ValueError("no candidates generated")

    if args.use_real_model and not octave_available() and not args.allow_proxy_fallback:
        raise RuntimeError("--use-real-model requires Octave. Install Octave or add --allow-proxy-fallback.")

    print(f"[optimizer] method={method} objective={objective} candidates={total} dry_run={args.dry_run} use_real_model={args.use_real_model} fast_mode={args.fast_mode}")
    sys.stdout.flush()

    start = time.time()
    results = []
    results_path = out / "results.jsonl"
    with results_path.open("w", encoding="utf-8") as f:
        for i, params in enumerate(candidates):
            scenario_dir = out / f"scenario_{i:04d}"
            run_meta = run_model(params, fixed, scenario_dir, args.dry_run, args.use_real_model, args.fast_mode, args.timeout_sec)
            if run_meta.get("metric") is not None:
                score = run_meta["metric"]
            elif args.use_real_model and not args.allow_proxy_fallback:
                score = float("-inf")
            else:
                score = score_proxy(params, objective)
            row = {
                "id": i,
                "params": params,
                "objective": objective,
                "score": score,
                "status": run_meta["status"],
                "artifact_dir": run_meta["artifact_dir"],
                "metric": run_meta.get("metric"),
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

    failed = sum(1 for r in results if str(r.get("status","")).startswith("error"))

    report = {
        "objective": objective,
        "search_method": method,
        "use_real_model": args.use_real_model,
        "fast_mode": args.fast_mode,
        "total_scenarios": total,
        "failed_scenarios": failed,
        "elapsed_seconds": round(time.time() - start, 3),
        "best": best,
        "top_k": top_k,
    }
    (out / "report.md").write_text("# LifeCycle Optimization Report\n\n" + json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"[optimizer] done best_score={best['score']:.6f} output_dir={out}")
    if args.use_real_model and failed > 0:
        print(f"[optimizer] warning: {failed} real-model scenarios failed; check scenario logs")
    print(f"[optimizer] absolute_output_dir={out.resolve()}")


if __name__ == "__main__":
    main()
