#!/usr/bin/env python3
"""Project-level CLI entrypoint for lifecycle optimization runs."""
import argparse
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OPTIMIZER = ROOT / "skills" / "lifecycle-optimizer" / "scripts" / "optimize.py"


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the lifecycle optimizer skill from the project root.")
    parser.add_argument("--config", default="skills/lifecycle-optimizer/assets/sample-case.json")
    parser.add_argument("--output-dir", default="outputs/lifecycle-optimizer")
    parser.add_argument("--real", action="store_true", help="Run the Octave lifecycle model.")
    parser.add_argument("--fast-mode", action="store_true", help="Use reduced grid/simulation settings for diagnostics.")
    parser.add_argument("--max-evals", type=int, default=None)
    parser.add_argument("--timeout-sec", type=int, default=1200)
    parser.add_argument("--progress-every", type=int, default=1)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    cmd = [
        sys.executable,
        str(OPTIMIZER),
        "--config",
        args.config,
        "--output-dir",
        args.output_dir,
        "--timeout-sec",
        str(args.timeout_sec),
        "--progress-every",
        str(args.progress_every),
    ]
    if args.real:
        cmd.append("--use-real-model")
    if args.fast_mode:
        cmd.append("--fast-mode")
    if args.max_evals is not None:
        cmd += ["--max-evals", str(args.max_evals)]
    if args.dry_run:
        cmd.append("--dry-run")

    return subprocess.call(cmd, cwd=str(ROOT))


if __name__ == "__main__":
    raise SystemExit(main())
