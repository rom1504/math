"""Replay checkpoint-six exact and diagnostic checks, preserving their scope."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
PROGRAMS = [
    ("transfer_adversary_same_order_orientation_repair_2026_09_06.py", "exact integer finite signing inequalities"),
    ("transfer_adversary_isotropic_response_extension_2026_09_06.py", "exact rational response certificates and extension inequalities"),
    ("transfer_adversary_gaussian_return_rank_2026_09_06.py", "exact integer tensor-rank identities"),
    ("transfer_seed_high_degree_feedback_checks_2026_09_06.py", "mixed: exact finite graph enumeration, floating rank/normalization/ascent diagnostics"),
    ("transfer_director_actual_marked_spin_ascent_2026_09_06.py", "exact integer/rational full-seed-cube checks on fifteen stored signings"),
    ("transfer_fresh_perfect_dyadic_inverse_classify.py", "exact finite Walsh inverse classification, no asymptotic consequence"),
]


def main():
    records = []
    output = ROOT / "computations/results/transfer_director_checkpoint6_replays_2026_09_06.json"
    for name, scope in PROGRAMS:
        source = ROOT / "computations" / name
        completed = subprocess.run([sys.executable, str(source)], cwd=ROOT,
                                   capture_output=True, text=True, check=True)
        records.append({"program": name, "scope": scope, "arguments": [],
                        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                        "stdout": completed.stdout, "stderr": completed.stderr,
                        "returncode": completed.returncode})
        output.write_text(json.dumps({
            "scope": "Replayed finite evidence; mathematical asymptotic proofs are separate",
            "programs_passed": len(records), "records": records}, indent=2)+"\n")
        print("PASS:", name, flush=True)
    print(json.dumps({"programs_passed": len(records), "output": str(output.relative_to(ROOT))}))


if __name__ == "__main__":
    main()
