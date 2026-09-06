"""Bounded timing and algebra audit of an attributed final catalogue level.

This does NOT replay the full level and does NOT establish catalogue completeness.
The pinned external module was read before execution; its main is not called.
"""
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "external_research/2026-09-06/robby955_mo_413935_61ac2268b6f9"


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    # These are standard BLAS runtime controls, set before importing NumPy.
    for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS"):
        os.environ[name] = "1"
    source = BASE / "verification/research_order15_certify.py"
    manifest = json.loads((ROOT / "research_archive/reviewed_external_dependencies.json").read_text())
    assert hashlib.sha256(source.read_bytes()).hexdigest() == manifest[str(source.relative_to(ROOT))]["sha256"]
    external = load_module("audited_external_order15", source)
    independent = load_module("independent_gray_replay", ROOT / "computations/transfer_adversary_external_order15_16_witness_verify_2026_09_06.py")
    catalogue = BASE / "verification/order15_level14_catalogue.g6.gz"
    payload = gzip.decompress(catalogue.read_bytes())
    assert hashlib.sha256(payload).hexdigest() == external.CATALOGUE_SHA256
    records = payload.splitlines()
    assert len(records) == external.CATALOGUE_RECORDS
    indices = sorted(set([0, len(records)-1] +
                         [j * (len(records)-1) // 255 for j in range(256)]))
    selected = [records[i] for i in indices]
    started = time.perf_counter()
    external._worker_init()
    # Warm-up is separate from the timed, fixed-size sample.
    external.replay_chunk(selected[:4])
    before = time.perf_counter()
    replay = external.replay_chunk(selected)
    elapsed = time.perf_counter() - before
    for raw in selected[::16]:
        matrix = external.graph6_signing(raw.decode("ascii"))
        cap = independent.scan(matrix)["cap"]
        energy = external.signing_energies(matrix, external._WORKER_SPINS)
        assert cap == int(abs(energy).max()) <= 25
    assert replay[0] == len(selected) and replay[1] == 0
    result = {
        "scope": "256-record final-level benchmark and 16 independent integer cap crosschecks; NOT full replay and NOT catalogue completeness",
        "external_commit": "61ac2268b6f9234c3ea268e7a7072ac72141f36b",
        "catalogue_compressed_sha256": hashlib.sha256(catalogue.read_bytes()).hexdigest(),
        "catalogue_uncompressed_sha256": hashlib.sha256(payload).hexdigest(),
        "catalogue_uncompressed_bytes": len(payload),
        "catalogue_records": len(records),
        "sample_indices": indices,
        "sample_result_processed_survivors_min_max": replay,
        "sample_elapsed_seconds_one_blas_thread": elapsed,
        "linear_extrapolation_seconds_one_worker": elapsed * len(records) / len(selected),
        "ideal_four_worker_extrapolation_seconds_not_a_guarantee": elapsed * len(records) / (4 * len(selected)),
        "integer_cap_crosschecks": len(selected[::16]),
        "blas_thread_controls": {name: os.environ[name] for name in
            ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS")},
        "total_elapsed_seconds": time.perf_counter() - started,
    }
    output = ROOT / "computations/results/transfer_adversary_external_order15_replay_assessment_2026_09_06.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
