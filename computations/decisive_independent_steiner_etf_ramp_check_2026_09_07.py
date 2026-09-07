"""Exact cyclic Steiner design, real frame, full-sign and ramp-bound checks."""

from collections import Counter
from fractions import Fraction as Q
from itertools import combinations
import json
import numpy as np


def main():
    v, k, r = 61, 5, 15
    base = [(1, 9, 20, 34, 58), (4, 14, 19, 36, 49), (13, 15, 16, 22, 56)]
    differences = Counter((a-b) % v for block in base for a in block for b in block if a != b)
    assert differences == Counter({a: 1 for a in range(1, v)})
    blocks = [tuple(sorted((x+t) % v for x in block)) for block in base for t in range(v)]
    assert len(blocks) == len(set(blocks)) == v*r//k
    pairs = Counter(pair for block in blocks for pair in combinations(block, 2))
    assert pairs == Counter({pair: 1 for pair in combinations(range(v), 2)})
    incidence = [[i for i, block in enumerate(blocks) if point in block] for point in range(v)]
    assert all(len(row) == r for row in incidence)

    h = np.ones((1, 1), dtype=np.int64)
    for _ in range(4):
        h = np.block([[h, h], [h, -h]])
    assert np.array_equal(h@h.T, 16*np.eye(16, dtype=np.int64))
    assert np.all(h[0] == 1)
    d, n = len(blocks), v*(r+1)
    f = np.zeros((d, n), dtype=np.int64)
    for point, rows in enumerate(incidence):
        for local, row in enumerate(rows):
            f[row, point*16:(point+1)*16] = h[local+1]
    assert np.array_equal(f@f.T, 80*np.eye(d, dtype=np.int64))
    gram_numerator = f.T@f
    assert np.all(np.diag(gram_numerator) == r)
    a = gram_numerator-r*np.eye(n, dtype=np.int64)
    assert np.array_equal(a, a.T)
    assert np.all(np.diag(a) == 0)
    assert np.all(np.abs(a+np.eye(n, dtype=np.int64)) == 1)
    assert int(np.sum(a*a)) == n*(n-1)
    # These row-Gram and rank identities imply eigenvalues 65,-15 exactly.
    assert 80-r == 4*r+5 == 65
    assert (4*r+5)*d-r*(n-d) == 0
    assert (4*r+5)**2*d+r*r*(n-d) == n*(n-1)

    expected_energy_lower = Q(7*n*(n-1), 22*r)
    width_upper = Q(n*(5*r+5), 4)
    midpoint_upper = Q(n*(4*r+5), 4)
    covariance_norm_upper = 1+Q(2*(4*r+5), 3*(r-1))
    trace_lower = n*(expected_energy_lower-width_upper)/(midpoint_upper*covariance_norm_upper)
    assert expected_energy_lower == Q(222040, 11)
    assert width_upper == 19520
    assert covariance_norm_upper == Q(86, 21)
    assert trace_lower == Q(61488, 6149)
    # The asymptotic trace constant is decreasing in pi; pi<22/7 suffices.
    pi_upper = Q(22, 7)
    assert (16-5*pi_upper)/(4*(pi_upper+8)) == Q(1, 156)
    print(json.dumps({
        "status": "all exact checks passed",
        "steiner_points": v, "steiner_blocks": d,
        "pair_coverage": len(pairs), "frame_shape": [d, n],
        "full_sign_order": n, "exact_eigenvalues": [65, -15],
        "expected_energy_lower": str(expected_energy_lower),
        "width_upper": str(width_upper),
        "covariance_norm_upper": str(covariance_norm_upper),
        "finite_zero_error_PSD_ramp_trace_lower": str(trace_lower),
        "asymptotic_trace_fraction_strict_lower": "1/156",
    }, indent=2))


if __name__ == "__main__":
    main()
