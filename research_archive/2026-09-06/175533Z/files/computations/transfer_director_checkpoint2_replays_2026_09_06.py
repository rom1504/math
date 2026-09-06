"""Preserve reproducible independent replay output; does not verify proofs."""

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
CHECKS = [
    "transfer_reconstruction_gibbs_cluster_exact_checks_2026_09_06.py",
    "transfer_adversary_cluster_cap_audit_2026_09_06.py",
    "transfer_reconstruction_homogeneous_edge_count_exact_checks_2026_09_06.py",
    "transfer_reconstruction_typical_hadamard_exact_checks_2026_09_06.py",
    "transfer_reconstruction_two_replica_rate_falsifier_2026_09_06.py",
    "transfer_seed_schur_twins_verify_2026_09_06.py",
    "transfer_seed_thin_tensor_verify_2026_09_06.py",
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "computations/results/transfer_director_checkpoint2_replays_2026_09_06.json")
    args = parser.parse_args()
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1",
               TMPDIR=str(ROOT / "tmp"))
    report = {"started_utc": datetime.now(timezone.utc).isoformat(),
              "scope": "Finite regression replays. Some checks are numerical diagnostics, explicitly labeled in their outputs; asymptotic conclusions require written proofs.",
              "checks": []}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for script in CHECKS:
        command = [sys.executable, str(ROOT / "computations" / script)]
        result = subprocess.run(command, cwd=ROOT, env=env, text=True,
                                capture_output=True, timeout=600)
        report["checks"].append({"script": script, "returncode": result.returncode,
                                 "stdout": result.stdout, "stderr": result.stderr})
        args.output.write_text(json.dumps(report, indent=2) + "\n")
        print(script, "PASS" if result.returncode == 0 else "FAIL", flush=True)
        if result.returncode:
            raise SystemExit(result.returncode)
    report["completed_utc"] = datetime.now(timezone.utc).isoformat()
    args.output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
