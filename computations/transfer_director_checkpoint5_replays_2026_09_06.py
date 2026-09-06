"""Director replay of checkpoint-five exact finite checks.

No solver timeout or external catalogue-completeness claim is certified here.
The individual scripts label their imported optimality assumptions explicitly.
"""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
PROGRAMS = [
    ["transfer_seed_clique_orientation_audit_2026_09_06.py"],
    ["transfer_seed_orientation_repair_audit_2026_09_06.py"],
    ["transfer_adversary_minimizer_isotropy_2026_09_06.py"],
    ["transfer_seed_subhalf_partition_diagnostic_2026_09_06.py"],
    ["transfer_seed_balanced_minimizer_screen_2026_09_06.py", "--external-witnesses"],
    ["transfer_seed_balanced_neighbor_probe_2026_09_06.py"],
    ["transfer_adversary_orientation_local_vs_global_2026_09_06.py", "--max-order", "8"],
]


def main():
    records = []
    output = ROOT / "computations/results/transfer_director_checkpoint5_replays_2026_09_06.json"
    for name, *args in PROGRAMS:
        source = ROOT / "computations" / name
        completed = subprocess.run([sys.executable, str(source), *args], cwd=ROOT,
                                   capture_output=True, text=True, check=True)
        records.append({"program": name, "arguments": args,
                        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                        "stdout": completed.stdout, "stderr": completed.stderr,
                        "returncode": completed.returncode})
        output.write_text(json.dumps({
            "scope": "Exact finite replay; no upgrade of imported global minima or asymptotic claims",
            "programs_passed": len(records), "records": records}, indent=2) + "\n")
        print("PASS:", name, flush=True)
    print(json.dumps({"programs_passed": len(records), "output": str(output.relative_to(ROOT))}))


if __name__ == "__main__":
    main()
