#!/usr/bin/env python3
"""Exact finite checks for Wave 47B coarea/ground-alignment audit."""

from __future__ import annotations

import itertools
import sys
from fractions import Fraction

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from envelope_block_cover_r27 import A6, A8, A9, projective_spins
from hard_center_branching_r45_check import down_up_kernel
from high_replica_pressure_r44_explore import geometry


def exact_boundary_data(a: np.ndarray, m: int):
    z, selector_tuples, q, parent_e, row, deficits, h0 = geometry(a, m)
    selectors = list(map(frozenset, selector_tuples))
    f = deficits == 0
    nsel = len(selectors)
    degree = [Fraction(int(v), nsel) for v in f.sum(axis=1)]
    return z, selectors, q, parent_e, row, deficits.astype(int), h0, f, degree


def weighted_ratio(deficits, kernel, delta, weights, cap: int) -> Fraction:
    """(1/2) E|u(S)-u(T)| / ((1-lambda) E u); normalizers cancel."""
    u = np.maximum(cap - deficits, 0)
    slack_sum = sum(
        (w * int(u[iz, i]) for iz, w in enumerate(weights) for i in range(u.shape[1])),
        Fraction(0),
    )
    variation = Fraction(0)
    for iz, w in enumerate(weights):
        if not w:
            continue
        for i in range(u.shape[1]):
            for j in range(u.shape[1]):
                variation += w * kernel[i, j] * abs(int(u[iz, i]) - int(u[iz, j]))
    return variation / (2 * delta * slack_sum)


def hard_escape_ratio(f, kernel, delta, weights, degree) -> Fraction:
    numerator = Fraction(0)
    denominator = Fraction(0)
    nsel = f.shape[1]
    for iz, w in enumerate(weights):
        denominator += w * degree[iz]
        collision = sum(
            (kernel[i, j] for i in range(nsel) for j in range(nsel)
             if f[iz, i] and f[iz, j]),
            Fraction(0),
        ) / nsel
        numerator += w * (degree[iz] - collision)
    return numerator / (delta * denominator)


def raw_cancellation(a, z, selectors, deficits, kernel, weights) -> Fraction:
    abs_c = np.empty_like(deficits)
    q_s = np.empty(len(selectors), dtype=int)
    for j, ss in enumerate(selectors):
        ids = sorted(ss)
        child = a[np.ix_(ids, ids)]
        abs_c[:, j] = np.abs(np.einsum("bi,ij,bj->b", z[:, ids], child, z[:, ids]))
        q_s[j] = int(deficits[0, j] + abs_c[0, j])
        assert np.all(deficits[:, j] + abs_c[:, j] == q_s[j])
    combined = Fraction(0)
    separate = Fraction(0)
    for iz, w in enumerate(weights):
        if not w:
            continue
        for i in range(len(selectors)):
            for j in range(len(selectors)):
                mass = w * kernel[i, j]
                combined += mass * abs(int(deficits[iz, i]) - int(deficits[iz, j]))
                separate += mass * (
                    abs(int(q_s[i]) - int(q_s[j]))
                    + abs(int(abs_c[iz, i]) - int(abs_c[iz, j]))
                )
    assert combined <= separate
    return combined / separate if separate else Fraction(0)


def oriented_q(a: np.ndarray, ids: tuple[int, ...], tau: int) -> int:
    b = a[np.ix_(ids, ids)]
    return max(tau * int(x @ b @ x) for x in projective_spins(len(ids)))


def port_identity_audit(a: np.ndarray, m: int) -> int:
    """Verify Q-|E| = min_tau(orientation gap + core regret + port regret)."""
    n = len(a)
    z_all = list(projective_spins(n))
    checked = 0
    for u in itertools.combinations(range(n), m - 1):
        outside = [v for v in range(n) if v not in u]
        for v in outside:
            ids = tuple(sorted(u + (v,)))
            q_tau = {tau: oriented_q(a, ids, tau) for tau in (-1, 1)}
            q_abs = max(q_tau.values())
            au = a[np.ix_(u, u)]
            avu = a[v, list(u)]
            for z in z_all:
                zu = z[list(u)]
                core_e = int(zu @ au @ zu)
                field = int(avu @ zu)
                full_e = core_e + 2 * int(z[v]) * field
                pieces = {}
                for tau in (-1, 1):
                    orient_gap = q_abs - q_tau[tau]
                    core_response_regret = q_tau[tau] - (tau * core_e + 2 * abs(field))
                    port_regret = 2 * (abs(field) - tau * int(z[v]) * field)
                    assert orient_gap >= 0 and core_response_regret >= 0 and port_regret >= 0
                    lhs = q_abs - tau * full_e
                    assert lhs == orient_gap + core_response_regret + port_regret
                    pieces[tau] = lhs
                    checked += 1
                assert min(pieces.values()) == q_abs - abs(full_e)
    return checked


EXPECTED = {
    ("A6", 2): (Fraction(1, 3), Fraction(7, 24), Fraction(1, 6)),
    ("A6", 3): (Fraction(1, 3), Fraction(7, 24), Fraction(1, 6)),
    ("A6", 4): (Fraction(1, 3), Fraction(7, 24), Fraction(1, 6)),
    ("A8", 3): (Fraction(879, 1000), Fraction(1401, 1760), Fraction(73, 110)),
    ("A8", 4): (Fraction(239, 250), Fraction(1321, 1540), Fraction(467, 660)),
    ("A8", 5): (Fraction(111, 100), Fraction(1209, 1232), Fraction(35, 44)),
    ("A9", 4): (Fraction(8413, 9072), Fraction(50383, 59346), Fraction(4321, 5184)),
    ("A9", 5): (Fraction(13609, 13608), Fraction(80845, 89019), Fraction(6919, 7776)),
    ("A9", 6): (Fraction(433, 378), Fraction(10154, 9891), Fraction(433, 432)),
}


def finite_audit(name: str, a: np.ndarray, m: int) -> None:
    z, selectors, q, parent_e, row, deficits, h0, f, degree = exact_boundary_data(a, m)
    assert 0 < h0 < 4 and np.all(deficits % 4 == 0)
    # These finite diagnostic laws happen to lie entirely in the old C2 class.
    assert np.all(row <= 2 * len(a) * (len(a) - 1))
    laws = {
        "uniform": [Fraction(1) for _ in z],
        "ground-lift": degree,
        "parent-ground": [Fraction(int(abs(int(e)) == q)) for e in parent_e],
    }
    caps = list(range(1, int(deficits.max()) + 2, 4))
    print(f"{name}: actual H={h0:.9f} lies in (0,4); caps audited={caps}")
    for ell in range(max(0, m - 3), m):
        kernel = down_up_kernel(selectors, len(a), m, ell)
        lam = Fraction(ell * (len(a) - m), m * (len(a) - ell)) if ell else Fraction(0)
        delta = 1 - lam
        got = tuple(hard_escape_ratio(f, kernel, delta, laws[k], degree) for k in laws)
        assert got == EXPECTED[(name, ell)]
        # At cap one, coarea is exactly the hard escape ratio.
        assert tuple(weighted_ratio(deficits, kernel, delta, laws[k], 1) for k in laws) == got
        cap_rows = {
            k: [weighted_ratio(deficits, kernel, delta, w, cap) for cap in caps]
            for k, w in laws.items()
        }
        cancellation = {k: raw_cancellation(a, z, selectors, deficits, kernel, w)
                        for k, w in laws.items()}
        print(f"  ell={ell}, delta={delta}: H<4 ratios {dict(zip(laws, got))}")
        print(f"    cap-regime ratios={cap_rows}")
        print(f"    raw |Delta D| / separate triangle cost={cancellation}")
    checks = port_identity_audit(a, m)
    print(f"  adjacent-port identities checked: {checks}")


def main() -> None:
    for args in (("A6", A6, 5), ("A8", A8, 6), ("A9", A9, 7)):
        finite_audit(*args)
    print("PASS: exact coarea ratios, structured center laws, cancellation, and port decomposition")


if __name__ == "__main__":
    main()
