#!/usr/bin/env python3
"""Customer-manager demo entrypoint."""
import argparse
import subprocess
import sys
from pathlib import Path

from advisor import DEFAULT_RESULTS_DIR, build_advice


ROOT = Path(__file__).resolve().parents[1]
OPTIMIZER = ROOT / "skills" / "lifecycle-optimizer" / "scripts" / "optimize.py"
DEFAULT_CONFIG = ROOT / "skills" / "lifecycle-optimizer" / "assets" / "sample-case.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the lifecycle advisor from natural-language client input.")
    parser.add_argument("text", nargs="*", help="Example: 45 years old, 1 million, stable returns")
    parser.add_argument("--results-dir", default=str(DEFAULT_RESULTS_DIR), help="Completed lifecycle result directory.")
    parser.add_argument("--output", default="outputs/client-advice-demo.md", help="Markdown output path.")
    parser.add_argument("--real", action="store_true", help="Run life_cycle.m through optimizer.py before building advice.")
    parser.add_argument("--config", default=str(DEFAULT_CONFIG), help="Optimizer config used with --real.")
    parser.add_argument("--real-output-dir", default="outputs/client-real-model", help="Optimizer artifact directory used with --real.")
    parser.add_argument("--fast-mode", action="store_true", help="Use smaller grids/simulation counts for real-model diagnostics.")
    parser.add_argument("--timeout-sec", type=int, default=1200, help="Timeout for the real Octave model run.")
    args = parser.parse_args()

    text = " ".join(args.text).strip()
    if not text:
        text = input("Client input, e.g. 45 years old, 1 million, stable returns: ").strip()

    results_dir = args.results_dir
    if args.real:
        real_output_dir = Path(args.real_output_dir)
        if not real_output_dir.is_absolute():
            real_output_dir = ROOT / real_output_dir
        cmd = [
            sys.executable,
            str(OPTIMIZER),
            "--config",
            str(args.config),
            "--output-dir",
            str(real_output_dir),
            "--client-text",
            text,
            "--use-real-model",
            "--timeout-sec",
            str(args.timeout_sec),
            "--progress-every",
            "1",
        ]
        if args.fast_mode:
            cmd.append("--fast-mode")
        code = subprocess.call(cmd, cwd=str(ROOT))
        if code != 0:
            return code
        results_dir = str(real_output_dir)

    result = build_advice(text, results_dir=results_dir)
    print(result["markdown"])

    output = Path(args.output)
    if not output.is_absolute():
        output = ROOT / output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(result["markdown"] + "\n", encoding="utf-8")
    print(f"\nSaved report: {output}")
    return 0 if result["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
