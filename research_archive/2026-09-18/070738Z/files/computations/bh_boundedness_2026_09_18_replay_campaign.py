#!/usr/bin/env python3
"""Reproducible closing replay, with explicit certificate/diagnostic labels.

Run using the project venv. This preserves stdout/stderr and hashes the
tracked proof-support scripts; it does not turn numerical checks into proofs.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time


REPLAYS = [
    ("bh_boundedness_barrier_2026_09_18.py", "exact_and_analytic_certificate_checks"),
    ("bh_boundedness_counterexamples_2026_09_18_tangent.py", "exact_Eisenstein_and_rational_certificate"),
    ("bh_boundedness_barrier_2026_09_18_tangent.py", "independent_exact_Eisenstein_certificate"),
    ("bh_boundedness_counterexamples_2026_09_18_normalization.py", "rational_checks_plus_labeled_high_precision_diagnostics"),
    ("bh_boundedness_counterexamples_2026_09_18_inheritance.py", "exact_finite_field_rank_certificate"),
    ("bh_boundedness_positive_2026_09_18_blockunit.py", "exact_Gaussian_integer_certificate"),
    ("bh_boundedness_positive_2026_09_18_bilinear.py", "exact_Gaussian_integer_certificate"),
    ("bh_boundedness_barrier_2026_09_18_velocity_converse.py", "exact_symbolic_spectral_certificate"),
    ("bh_boundedness_barrier_2026_09_18_velocity_exact.py", "exact_algebraic_root_certificate"),
    ("bh_boundedness_director_2026_09_18_quadratic_phase.py", "numerical_diagnostics_only"),
    ("bh_boundedness_director_2026_09_18_flat_witness.py", "numerical_diagnostics_only"),
    ("bh_boundedness_counterexamples_2026_09_18_phase_l2_compiler.py", "numerical_diagnostics_only"),
    ("bh_boundedness_counterexamples_2026_09_18_checks.py", "finite_census_and_labeled_tensor_checks"),
]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent.parent
    scratch = root / "tmp" / "bh_boundedness_2026_09_18"
    scratch.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    for key in ("TMPDIR", "TEMP", "TMP"):
        env[key] = str(scratch)
    records = []
    started = time.time()
    for filename, status in REPLAYS:
        path = root / "computations" / filename
        tick = time.time()
        try:
            completed = subprocess.run([sys.executable, str(path)], cwd=str(root),
                                       env=env, capture_output=True, text=True,
                                       timeout=180)
            record = {"script": str(path.relative_to(root)), "evidence_type": status,
                      "returncode": completed.returncode,
                      "stdout": completed.stdout, "stderr": completed.stderr,
                      "elapsed_seconds": time.time()-tick}
        except subprocess.TimeoutExpired as exc:
            record = {"script": str(path.relative_to(root)), "evidence_type": status,
                      "returncode": None, "error": "TIMEOUT_NOT_MATHEMATICAL_EVIDENCE",
                      "elapsed_seconds": time.time()-tick}
        records.append(record)
        print(json.dumps({k: record[k] for k in ("script", "returncode", "elapsed_seconds")}), flush=True)
    sources = {}
    for path in sorted((root / "computations").glob("bh_boundedness_*2026_09_18*.py")):
        sources[str(path.relative_to(root))] = hashlib.sha256(path.read_bytes()).hexdigest()
    report = {"status": "PASS" if all(r["returncode"] == 0 for r in records) else "REVIEW_FAILURES",
              "python": sys.version, "started_utc": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(started)),
              "elapsed_seconds": time.time()-started, "records": records, "source_sha256": sources}
    args.output.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({"status": report["status"], "replays": len(records), "output": str(args.output)}), flush=True)
    if report["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
