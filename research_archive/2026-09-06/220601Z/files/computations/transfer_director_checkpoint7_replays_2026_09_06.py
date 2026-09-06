"""Replay final-checkpoint finite checks, preserving exact/numerical scope."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
PROGRAMS = [
    ('transfer_director_boolean_cubic_chain_exact_2026_09_06.py',
     'Exact integer/Fraction cubic-chain and fourth-moment regressions'),
    ('transfer_director_rich_holes_exact_2026_09_06.py',
     'Exact integer/Fraction full-cube rich-response counterexample'),
    ('transfer_seed_trigonometric_cut_checks_2026_09_06.py',
     'Floating trigonometric identities and matrix-norm regressions; not exact certification'),
    ('transfer_fresh_untruncated_first_marked_verify.py',
     'Exact integer/Fraction first-marked covariance, influence and oriented ascent'),
    ('transfer_adversary_nuclear_alias_counterexample_2026_09_06.py',
     'Finite covariance regressions plus floating bounded-response alias diagnostic'),
    ('transfer_seed_local_sobolev_product_checks_2026_09_06.py',
     'Finite numerical Sobolev/product identities and coefficient checks; no asymptotic certificate'),
]


def main():
    records = []
    output = ROOT / 'computations/results/transfer_director_checkpoint7_replays_2026_09_06.json'
    for name, scope in PROGRAMS:
        source = ROOT / 'computations' / name
        completed = subprocess.run([sys.executable, str(source)], cwd=ROOT,
                                   capture_output=True, text=True, check=True)
        records.append({'program': name, 'scope': scope, 'arguments': [],
                        'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
                        'stdout': completed.stdout, 'stderr': completed.stderr,
                        'returncode': completed.returncode})
        output.write_text(json.dumps({
            'scope': 'Finite checks only; asymptotic proofs are independent',
            'programs_passed': len(records), 'records': records}, indent=2) + '\n')
        print('PASS:', name, flush=True)


if __name__ == '__main__':
    main()
