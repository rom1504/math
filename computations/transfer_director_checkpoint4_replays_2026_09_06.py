"""Reproduce checkpoint-four finite checks; retain evidentiary labels."""
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
PROGRAMS = [
    'transfer_reconstruction_threshold_entropy_exact_checks_2026_09_06.py',
    'transfer_seed_threshold_entropy_resampling_verify_2026_09_06.py',
    'transfer_director_weighted_information_checks_2026_09_06.py',
    'transfer_director_threshold_cap_integration_2026_09_06.py',
    'transfer_fresh_paley_isotropic_extension_verify.py',
]


def main():
    outputs = []
    for name in PROGRAMS:
        result = subprocess.run([sys.executable, str(ROOT / 'computations' / name)],
                                cwd=ROOT, check=True, capture_output=True, text=True)
        outputs.append({'program': name, 'result': json.loads(result.stdout),
                        'stderr': result.stderr})
    output = ROOT / 'computations/results/transfer_director_checkpoint4_replays_2026_09_06.json'
    output.write_text(json.dumps({'status': 'all programs passed; numerical portions remain numerical',
                                  'outputs': outputs}, indent=2) + '\n')
    print(json.dumps({'programs_passed': len(outputs), 'saved': str(output.relative_to(ROOT))}))


if __name__ == '__main__':
    main()
