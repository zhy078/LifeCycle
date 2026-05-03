# test_skill.py
import sys
import os
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from optimize import optimize

# Use the sample config from assets
script_dir = Path(__file__).resolve().parent
project_root = script_dir.parent.parent.parent
config_path = project_root / "skills" / "lifecycle-optimizer" / "assets" / "input-template.json"
output_dir = project_root / "outputs" / "test_run"

if not config_path.exists():
    print(f"Error: Config file not found at {config_path}")
    sys.exit(1)

result = optimize(
    config_path=str(config_path),
    output_dir=str(output_dir),
    fast_mode=True
)

print("Test passed!")
print(f"Best score: {result['best']['score']}")
print(f"Total scenarios: {result['total']}")
print(f"Runtime: {result['runtime']}s")
print(f"Output directory: {result['output_dir']}")