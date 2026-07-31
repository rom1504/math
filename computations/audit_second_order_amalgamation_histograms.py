#!/usr/bin/env python3
"""Audit a histogram-level second-order bridge amalgamation criterion.

For fixed child signings A,D and bridge W, write

  eta(x,y) = M(A)+M(D)-|H_A(x)+H_D(y)|,
  zeta(u,v) = L(W)-|u^T W v|.

If P,Q are uniformly random diagonal switchings, then for every fixed child
pair (x,y), zeta(Px,Qy) has the uniform bridge-deficit histogram.  A union
bound therefore certifies a switching with eta+zeta >= g everywhere whenever

  sum_{e+z<g} N_eta(e) N_zeta(z) < 2^(m+n-2).

The resulting parent cap is at most M(A)+M(D)+L(W)-g.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bridge_block_cpsat import load_matrix, one_copy_energies


DEFAULT_INPUTS = (
    Path("computations/results/bridge_5_5_sign1.json"),
    Path("computations/results/bridge_6_6_sign1_cap18.json"),
    Path("computations/results/bridge_6_7_sign1_cap20.json"),
)


def histogram(values: np.ndarray) -> dict[str, int]:
    return {
        str(key): value
        for key, value in sorted(Counter(map(int, values.ravel())).items())
    }


def audit_arrays(
    source: str, a: np.ndarray, d: np.ndarray, w: np.ndarray
) -> dict[str, object]:
    spins_a, energies_a = one_copy_energies(a)
    spins_d, energies_d = one_copy_energies(d)
    spins_a = spins_a.astype(np.int64)
    spins_d = spins_d.astype(np.int64)
    m_a = int(np.max(np.abs(energies_a)))
    m_d = int(np.max(np.abs(energies_d)))
    internal = np.abs(energies_a[:, None] + energies_d[None, :])
    eta = m_a + m_d - internal
    cross = spins_a @ w @ spins_d.T
    bridge_cap = int(np.max(np.abs(cross)))
    zeta = bridge_cap - np.abs(cross)
    state_count = int(eta.size)
    if state_count != 1 << (len(a) + len(d) - 2):
        raise AssertionError("projective state count mismatch")

    count_eta = Counter(map(int, eta.ravel()))
    count_zeta = Counter(map(int, zeta.ravel()))
    certified_gain = 0
    expectations: list[dict[str, object]] = []
    for gain in range(m_a + m_d + bridge_cap + 2):
        bad_numerator = sum(
            count_e * count_z
            for e, count_e in count_eta.items()
            for z, count_z in count_zeta.items()
            if e + z < gain
        )
        if bad_numerator < state_count:
            certified_gain = gain
        expectations.append(
            {
                "gain": gain,
                "expected_bad_pair_numerator": bad_numerator,
                "expected_bad_pair_denominator": state_count,
            }
        )

    best_cap = m_a + m_d + bridge_cap
    best_switch: tuple[list[int], list[int]] | None = None
    for p in spins_a:
        left = spins_a * p
        for q in spins_d:
            right = spins_d * q
            switched_cross = np.abs(left @ w @ right.T)
            cap = int(np.max(internal + switched_cross))
            if cap < best_cap:
                best_cap = cap
                best_switch = (list(map(int, p)), list(map(int, q)))
    exact_switching_gain = m_a + m_d + bridge_cap - best_cap
    if exact_switching_gain < certified_gain:
        raise AssertionError((source, exact_switching_gain, certified_gain))

    return {
        "source": source,
        "orders": [len(a), len(d)],
        "child_caps": [m_a, m_d],
        "bridge_cap": bridge_cap,
        "independent_amalgamation_cap": m_a + m_d + bridge_cap,
        "projective_pair_count": state_count,
        "internal_sum_deficit_histogram": histogram(eta),
        "bridge_deficit_histogram": histogram(zeta),
        "histogram_first_moment_certified_gain": certified_gain,
        "histogram_certified_parent_cap": (
            m_a + m_d + bridge_cap - certified_gain
        ),
        "exhaustive_best_switching_gain": exact_switching_gain,
        "exhaustive_best_switching_parent_cap": best_cap,
        "best_switch": best_switch,
        "expectation_table": expectations,
    }


def audit(source: Path) -> dict[str, object]:
    payload = json.loads(source.read_text())
    if "bridge" not in payload:
        raise ValueError(f"{source} has no bridge witness")
    a = load_matrix(Path(payload["child_a"]), payload.get("child_a_class"))
    d = load_matrix(Path(payload["child_b"]), payload.get("child_b_class"))
    d = int(payload["sign_b"]) * d
    w = np.asarray(payload["bridge"], dtype=np.int64)
    return audit_arrays(str(source), a, d, w)


def sylvester(order: int) -> np.ndarray:
    matrix = np.ones((1, 1), dtype=np.int64)
    while len(matrix) < order:
        matrix = np.block([[matrix, matrix], [matrix, -matrix]])
    if len(matrix) != order:
        raise ValueError("Sylvester order must be a power of two")
    return matrix


def audit_sylvester(order: int) -> dict[str, object]:
    child = load_matrix(Path(f"computations/results/exact_m{order}.json"))
    return audit_arrays(
        f"structured Sylvester bridge with exact_m{order} children",
        child,
        child,
        sylvester(order),
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="*", type=Path, default=list(DEFAULT_INPUTS))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    records = [audit(source) for source in args.inputs]
    records.extend(audit_sylvester(order) for order in (4, 8))
    output = {
        "schema": "quadratic-signing-second-order-amalgamation-histogram-v1",
        "classification": (
            "proved finite histogram certificate plus exhaustive switching audit"
        ),
        "criterion": (
            "sum_{e+z<g} N_eta(e)N_zeta(z)<2^(m+n-2) implies a bridge "
            "switching with parent cap <= M_A+M_D+L_W-g"
        ),
        "records": records,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    for record in records:
        print(
            f"{record['orders']}: histogram_gain="
            f"{record['histogram_first_moment_certified_gain']} "
            f"best_switching_gain={record['exhaustive_best_switching_gain']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
