"""Checkpoint-three independent finite regressions, with preserved output."""

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
CHECKS = [
    "transfer_seed_spectral_gaussian_verify_2026_09_06.py",
    "transfer_seed_hadamard_graph_verify_2026_09_06.py",
    "transfer_adversary_sparse_hadamard_two_over_pi_2026_09_06.py",
    "transfer_seed_intermediate_tensor_moments_2026_09_06.py",
    "transfer_reconstruction_correlated_cluster_exact_checks_2026_09_06.py",
    "transfer_reconstruction_adaptive_cluster_exact_checks_2026_09_06.py",
    "transfer_reconstruction_stationary_logdet_exact_checks_2026_09_06.py",
    "transfer_adversary_complement_exposure_2026_09_06.py",
    "transfer_adversary_seed_only_conditional_tail_2026_09_06.py",
    "transfer_director_selector_convexity_2026_09_06.py",
    "transfer_adversary_fixed_retention_covariance_2026_09_06.py",
    "transfer_seed_exceptional_selector_verify_2026_09_06.py",
    "transfer_seed_retention_transition_verify_2026_09_06.py",
    "transfer_reconstruction_gaussian_threshold_exact_checks_2026_09_06.py",
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "computations/results/transfer_director_checkpoint3_replays_2026_09_06.json")
    args = parser.parse_args()
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1",
               TMPDIR=str(ROOT / "tmp"))
    report = {"started_utc": datetime.now(timezone.utc).isoformat(),
              "scope": "Finite regressions, not formal verification. Numerical parts retain their original labels. Asymptotic statements depend on the written proofs.",
              "checks": []}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for script in CHECKS:
        result = subprocess.run([sys.executable, str(ROOT / "computations" / script)],
                                cwd=ROOT, env=env, text=True,
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
