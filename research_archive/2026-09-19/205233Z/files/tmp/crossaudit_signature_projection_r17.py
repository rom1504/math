#!/usr/bin/env python3
"""Independent exhaustive audit of the Wave-17 signature projection memo."""

from __future__ import annotations

from itertools import combinations

import numpy as np

from audit_signature_partitions_r17 import (
    A8,
    A9,
    active_defect,
    mask_to_spin,
    qnorm,
    signature_partition,
    span,
    spin_to_mask,
    spins,
    split,
)
from check_response_dual_r16 import energies


def audit_matrix(A: np.ndarray, expected_groups: int) -> dict:
    n = len(A)
    Z = spins(n)
    ea = energies(A)
    qa = qnorm(A)
    grounds_a = sorted({spin_to_mask(z) for z, e in zip(Z, ea) if abs(int(e)) == qa})

    seen = set()
    rank_one_count = 0
    all_state_checks = 0
    all_cross_ground_checks = 0
    best = []
    for rank in (1, 2, 3):
        for generators in combinations(grounds_a, rank):
            group = span(generators)
            if len(group) != 1 << rank or group in seen:
                continue
            seen.add(group)
            group_spins = [mask_to_spin(g, n) for g in group]
            cells = signature_partition(group, n)
            C, D = split(A, cells)

            # Matrix-level Reynolds identity, stronger than checking one orbit.
            reynolds_num = sum((A * np.outer(g, g) for g in group_spins), np.zeros_like(A))
            assert np.array_equal(reynolds_num, len(group) * D)
            assert qnorm(D) <= qa

            ec = energies(C)
            ed = energies(D)
            qc = qnorm(C)
            # The matrix identity implies the quadratic identity for every
            # real vector, hence in particular for all projective spins.  We
            # still check the pointwise decomposition in one vectorized pass.
            assert np.array_equal(ea, ec + ed)
            all_state_checks += len(Z)

            if rank == 1:
                rank_one_count += 1
                g = group_spins[1]
                assert np.array_equal(D * np.outer(g, g), D)
                assert np.array_equal(C * np.outer(g, g), -C)
                qc2, response, gamma, _ = active_defect(C, D, 0)
                assert qc2 == qc
                abs_d = [abs(int(d)) for c, d in zip(ec, ed) if abs(int(c)) == qc]
                assert response == max(abs_d)
                assert gamma == 0
                assert qc <= qa

            if qc == 0:
                # The variance quotient is undefined, but the active defect is
                # still zero because every state is C-active.
                assert active_defect(C, D, 0)[2] == 0
                continue

            active_cache = {}
            for z, c0, d0 in zip(Z, ec, ed):
                c0 = int(c0)
                if abs(c0) != qc:
                    continue
                sigma = 1 if c0 > 0 else -1
                y = [sigma * int((z * g) @ C @ (z * g)) for g in group_spins]
                assert sum(y) == 0 and max(y) == qc
                delta = qc + min(y)
                assert 0 <= delta <= qc
                a = sigma * int(d0)
                if delta not in active_cache:
                    active_cache[delta] = active_defect(C, D, delta)
                qc2, response, gamma, _ = active_cache[delta]
                assert qc2 == qc
                assert response >= max(a, -a - delta)
                assert 2 * gamma <= delta

                # Mean-zero Bhatia--Davis, without division or floating point.
                v_num = sum(value * value for value in y)
                assert v_num <= len(group) * qc * (qc - delta)
                # Equivalently delta <= qc - v/qc.
                assert delta * qc * len(group) <= qc * qc * len(group) - v_num
                all_cross_ground_checks += 1

            cell_sizes = tuple(sorted(map(len, cells)))
            S = sum(qnorm(A[np.ix_(cell, cell)]) for cell in cells)
            qmin = {1: 0, 2: 2, 3: 6, 4: 8, 5: 8, 6: 10, 7: 18, 8: 20, 9: 24}
            X = S - sum(qmin[len(cell)] for cell in cells)
            best.append((X, rank, cell_sizes, qc, int(np.mean([g @ A @ g for g in group_spins]))))

    assert len(seen) == expected_groups
    return {
        "groups": len(seen),
        "rank_one": rank_one_count,
        "state_checks": all_state_checks,
        "cross_ground_checks": all_cross_ground_checks,
        "best": sorted(best, reverse=True)[:5],
    }


def main() -> None:
    print("A8", audit_matrix(A8, 79))
    print("A9", audit_matrix(A9, 1917))
    print("PASS: full-state Reynolds, rank-one, antipodal, and variance audits")


if __name__ == "__main__":
    main()
