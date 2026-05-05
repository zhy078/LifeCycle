#!/usr/bin/env python3
import argparse
import json
import shutil
import subprocess
from pathlib import Path


def run(cmd, cwd):
    p = subprocess.run(cmd, cwd=str(cwd), capture_output=True, text=True)
    return p.returncode, p.stdout, p.stderr


def main():
    ap = argparse.ArgumentParser(description="Integration test for lifecycle-optimizer real model path")
    ap.add_argument("--config", default="skills/lifecycle-optimizer/assets/sample-case.json")
    ap.add_argument("--output-dir", default="outputs/lifecycle-integration-test")
    ap.add_argument("--max-evals", type=int, default=2)
    ap.add_argument("--timeout-sec", type=int, default=180)
    ap.add_argument("--strict-real", action="store_true", help="fail if octave missing or scenario not ok")
    args = ap.parse_args()

    repo = Path(__file__).resolve().parents[3]
    out = (repo / args.output_dir).resolve()
    if out.exists():
        shutil.rmtree(out)

    octave_exists = shutil.which("octave") is not None

    cmd = [
        "py" if shutil.which("py") else "python3",
        "skills/lifecycle-optimizer/scripts/optimize.py",
        "--config", args.config,
        "--output-dir", args.output_dir,
        "--use-real-model",
        "--fast-mode",
        "--max-evals", str(args.max_evals),
        "--timeout-sec", str(args.timeout_sec),
        "--progress-every", "1",
    ]

    if not args.strict_real:
        cmd.append("--allow-proxy-fallback")

    code, stdout, stderr = run(cmd, repo)
    print("=== optimize stdout ===")
    print(stdout)
    if stderr.strip():
        print("=== optimize stderr ===")
        print(stderr)

    if code != 0:
        raise SystemExit(f"FAIL: optimize command exited with code={code}")

    results_file = out / "results.jsonl"
    best_file = out / "best_params.json"
    report_file = out / "report.md"
    for f in [results_file, best_file, report_file]:
        if not f.exists():
            raise SystemExit(f"FAIL: missing output file: {f}")

    rows = [json.loads(line) for line in results_file.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not rows:
        raise SystemExit("FAIL: results.jsonl has no rows")

    # Core checks
    for i, row in enumerate(rows):
        for k in ["status", "run_seconds", "year_files_count", "octave_command"]:
            if k not in row:
                raise SystemExit(f"FAIL: row {i} missing key: {k}")

    if args.strict_real:
        if not octave_exists:
            raise SystemExit("FAIL: strict-real enabled but octave not found in PATH")
        bad = [r for r in rows if r.get("status") != "ok" or r.get("year_files_count", 0) <= 0]
        if bad:
            raise SystemExit(f"FAIL: strict-real expects all ok rows with year files, got {len(bad)} bad rows")

    print("PASS: integration outputs look valid")
    print(f"Output dir: {out}")


if __name__ == "__main__":
    main()
