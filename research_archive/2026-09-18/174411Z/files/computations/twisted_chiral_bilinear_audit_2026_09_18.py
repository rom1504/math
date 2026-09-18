#!/usr/bin/env python3
"""Exact finite Boolean bilinear norms and a deterministic full-sign obstruction."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def projective_spins(n):
    ids = np.arange(1 << (n - 1), dtype=np.uint64)
    tail = 2 * ((ids[:, None] >> np.arange(n - 1, dtype=np.uint64)) & 1).astype(np.int64) - 1
    return np.column_stack((np.ones(len(ids), dtype=np.int64), tail))


def signing_matrix(value):
    if not isinstance(value, list) or not value or not isinstance(value[0], list):
        return None
    n = len(value)
    if not 2 <= n <= 18 or not all(isinstance(row, list) and len(row) == n for row in value):
        return None
    try:
        a = np.asarray(value, dtype=np.int64)
    except (ValueError, TypeError):
        return None
    if not np.array_equal(a, a.T) or np.any(a.diagonal()):
        return None
    if not np.all(abs(a + np.eye(n, dtype=np.int64)) == 1):
        return None
    return a


def matrices_in(value, prefix=""):
    a = signing_matrix(value)
    if a is not None:
        yield prefix, a
    elif isinstance(value, dict):
        for key, child in value.items():
            yield from matrices_in(child, prefix + "/" + str(key))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            yield from matrices_in(child, prefix + "/" + str(index))


def profile(a, minimize_diagonal=False):
    n = len(a)
    x = projective_spins(n)
    fields = x @ a
    h = (fields * x).sum(axis=1) // 2
    row_values = abs(fields).sum(axis=1)
    k = int(np.argmax(row_values))
    beta = int(row_values[k])
    q = int(abs(h).max())
    other = np.where(fields[k] >= 0, 1, -1)
    assert int(other @ a @ x[k]) == beta
    record = {"n": n, "Q": q, "beta": beta, "beta_over_Q": beta / q,
              "twist_independent_bound_beta_minus_n": beta - n,
              "beta_witness_y": x[k].tolist(), "beta_witness_x": other.tolist(),
              "matrix": a.tolist(), "exact_spin_count": len(x),
              "beta_exceeds_2sqrt2_Q": beta * beta > 8 * q * q}
    if minimize_diagonal:
        # |a_i+d_i| = |a_i| + d_i sign(a_i) for nonzero integer a_i,
        # and equals 1 when a_i=0. Every summand is exact in float64.
        local = fields * x
        slopes = np.sign(local)
        offsets = row_values + (local == 0).sum(axis=1)
        diagonal_spins = np.concatenate((x, -x), axis=0)
        best = None
        for start in range(0, len(diagonal_spins), 256):
            ds = diagonal_spins[start:start + 256]
            values = ds.astype(float) @ slopes.T.astype(float) + offsets
            caps = values.max(axis=1)
            i = int(np.argmin(caps))
            if best is None or caps[i] < best[0]:
                best = int(caps[i]), ds[i].tolist()
        record["gamma_min_diagonal_bilinear"] = best[0]
        record["gamma_minimizing_diagonal"] = best[1]
        record["gamma_diagonal_count"] = len(diagonal_spins)
    return record


def hadamard(n):
    h = np.ones((1, 1), dtype=np.int64)
    while len(h) < n:
        h = np.block([[h, h], [h, -h]])
    assert len(h) == n
    return h


def obstruction(m=64):
    r = hadamard(m)
    k = np.ones((m, m), dtype=np.int64) - np.eye(m, dtype=np.int64)
    a = np.block([[k, r], [r.T, -k]])
    assert np.array_equal(r @ r.T, m * np.eye(m, dtype=np.int64))
    x = np.ones(2 * m, dtype=np.int64)
    y = np.r_[np.ones(m, dtype=np.int64), -np.ones(m, dtype=np.int64)]
    beta_witness = int(x @ a @ y)
    assert beta_witness == 2 * m * (m - 1)
    cap_upper = m * m / 2 + m ** 1.5
    child_lower = beta_witness - 2 * m
    return {"m": m, "n": 2 * m, "Q_upper_bound": cap_upper,
            "beta_lower_bound": beta_witness, "every_twist_child_lower_bound": child_lower,
            "target_2sqrt2_Q_upper_bound": 2 ** 1.5 * cap_upper,
            "matrix": a.tolist(), "x": x.tolist(), "y": y.tolist(),
            "full_sign_validated": bool(np.array_equal(a, a.T) and
            not np.any(a.diagonal()) and
            np.all(abs(a + np.eye(2*m, dtype=np.int64)) == 1))}


def padded_obstruction_certificate(k=12):
    """Succinct exact certificate: no allocation of the huge signing."""
    assert k % 2 == 0
    n, root, m = 4**k, 2**k, 2**(3*k//2)
    assert 2*m <= n and m % 2 == 0 and m*m == n*root
    q_upper = n*root//2 + m*m//2 + m*root + m
    beta_lower = n*root + 2*m*m - 6*m*root - 4*m
    child_lower = beta_lower-n
    strict_gap_certificate = child_lower**2 - 8*q_upper**2
    assert strict_gap_certificate > 0
    return {"k": k, "N": n, "m": m, "sqrt_N": root, "K": 1,
            "Q_upper_bound": q_upper, "beta_lower_bound": beta_lower,
            "every_twist_child_lower_bound": child_lower,
            "child_lower_squared_minus_8_Q_upper_squared": strict_gap_certificate,
            "matrix_definition_zero_based": [
                "A[i,i]=0",
                "A[i,j]=+1 if 0<=i,j<m and i!=j",
                "A[i,j]=-1 if m<=i,j<2m and i!=j",
                "otherwise A[i,j]=(-1)^popcount(i & j)"],
            "bent_vector_definition": "z_i=(-1)^(sum_j bit_(2j)(i)*bit_(2j+1)(i))",
            "matrix_explicitly_materialized": False,
            "certificate_type": "exact integer bounds proved in adversary artifact; no optimum claim"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--diagonal-max-n", type=int, default=0)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    paths = [ROOT / "computations/results" / ("exact_m%d.json" % n) for n in range(3, 11)]
    paths += sorted((ROOT / "computations/results").glob("m*_minimizer_orbits.json"))
    paths += [ROOT / p for p in [
        "artifacts/dependent_profile_recovery_witness.json",
        "artifacts/dependent_profile_recovery_witness(1).json",
        "computations/results/certified_m11_m12.json",
        "computations/results/certified_m13_m14.json",
        "computations/results/nested_10_in_11_cap17.json",
        "computations/results/extension_nested_m11_to_12.json",
        "computations/results/extension_nested_m12_to_13.json",
        "computations/results/bridge_6_7_sign1_cap20.json",
        "computations/results/conference_completion_m13.json",
        "computations/results/heuristic_m14_from_conference.json",
        "computations/results/heuristic_m15_cap25_search.json",
        "computations/results/heuristic_m16_from_class1_double.json",
        "computations/results/transfer_adversary_external_order15_16_witness_verify_2026_09_06.json",
    ]]
    found = ({p["sha256_int8"]: p for p in json.loads(args.output.read_text())["profiles"]}
             if args.resume and args.output.exists() else {})
    for path in paths:
        if not path.exists():
            continue
        value = json.loads(path.read_text())
        for location, a in matrices_in(value):
            digest = hashlib.sha256(a.astype(np.int8).tobytes()).hexdigest()
            source = {"path": str(path.relative_to(ROOT)), "json_location": location}
            if digest in found:
                if source not in found[digest]["sources"]:
                    found[digest]["sources"].append(source)
            else:
                p = profile(a, len(a) <= args.diagonal_max_n)
                p.update({"sources": [source], "sha256_int8": digest})
                found[digest] = p
                print(json.dumps({k: v for k, v in p.items() if k in [
                    "n", "Q", "beta", "beta_over_Q", "gamma_min_diagonal_bilinear"]}), flush=True)
    result = {"normalization": "H=sum_i<j", "profiles": list(found.values()),
              "hadamard_two_clique_obstruction": obstruction(),
              "padded_obstruction_certificate": padded_obstruction_certificate()}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
