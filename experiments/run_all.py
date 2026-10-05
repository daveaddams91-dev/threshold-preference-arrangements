"""Run the full experiment suite in order and summarise the outcome.

Each stage is a separate process so that a failure in one does not mask the
others, and each prints its own verdict.  The exit status is non-zero if any
stage fails.

Run:  python experiments/run_all.py
"""
from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

STAGES = [
    ("unit tests", [sys.executable, "-m", "pytest", "tests", "-q"]),
    ("exact intersection counts vs independent numerical search",
     [sys.executable, "experiments/validate_intersections.py"]),
    ("exact region counts for n = 2..6",
     [sys.executable, "experiments/region_counts.py", "--trials", "12"]),
    ("closed form and ranking profiles",
     [sys.executable, "experiments/verify_theorem.py"]),
    ("independent grid chamber counts",
     [sys.executable, "experiments/grid_check.py"]),
    ("figures", [sys.executable, "experiments/make_figures.py"]),
]


def main() -> int:
    failures = []
    for name, cmd in STAGES:
        print("=" * 78)
        print(f"STAGE: {name}")
        print(f"  $ {' '.join(cmd)}")
        print("=" * 78, flush=True)
        t0 = time.time()
        proc = subprocess.run(cmd, cwd=ROOT)
        dt = time.time() - t0
        status = "ok" if proc.returncode == 0 else f"FAILED ({proc.returncode})"
        print(f"\n--- {name}: {status} in {dt:.1f}s ---\n", flush=True)
        if proc.returncode != 0:
            failures.append(name)
    print("=" * 78)
    if failures:
        print("FAILED STAGES:", ", ".join(failures))
        return 1
    print(f"all {len(STAGES)} stages passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())