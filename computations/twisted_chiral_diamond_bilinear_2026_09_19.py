#!/usr/bin/env python3
"""Exact rational 4/3-bilinear certificate for the coordinate-diamond norm.

The certificate was discovered by LP but is verified here by integer matrix
arithmetic; this script has no optimizer dependency. Types are (i,j,r),
where a=i,b=1-i,c=j*r,d=(1-j)*r after coordinate switching.
"""

import itertools
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
TERMS = [
    (2, [-1,-1,-1,-1,-1,-1,-1,1], [1,-1,1,-1,1,-1,-1,-1]),
    (2, [-1,-1,-1,-1,-1,-1,1,-1], [1,-1,-1,1,1,-1,-1,-1]),
    (4, [-1,-1,-1,-1,-1,1,1,1], [1,-1,1,1,1,-1,1,1]),
    (4, [-1,-1,-1,-1,1,-1,1,-1], [1,-1,-1,-1,1,-1,-1,1]),
    (2, [-1,-1,-1,-1,1,1,-1,1], [1,-1,1,-1,1,-1,1,1]),
    (2, [-1,-1,-1,-1,1,1,1,-1], [1,-1,-1,1,1,-1,1,1]),
    (1, [-1,-1,-1,1,-1,1,-1,1], [1,1,1,-1,1,-1,1,-1]),
    (2, [-1,-1,-1,1,1,-1,-1,1], [-1,-1,1,-1,1,-1,1,-1]),
    (4, [-1,-1,-1,1,1,1,-1,1], [1,-1,1,-1,1,1,1,-1]),
    (4, [-1,-1,1,-1,-1,-1,-1,-1], [1,-1,-1,1,-1,-1,-1,-1]),
    (1, [-1,-1,1,-1,-1,1,1,-1], [1,1,-1,1,1,-1,-1,1]),
    (2, [-1,1,-1,-1,-1,-1,-1,-1], [1,-1,1,1,-1,-1,-1,-1]),
    (1, [-1,1,1,1,-1,-1,-1,-1], [-1,1,1,1,-1,-1,-1,-1]),
    (1, [-1,1,1,1,1,1,1,1], [-1,1,1,1,1,1,1,1]),
]


def main():
    types = list(itertools.product([0, 1], [0, 1], [-1, 1]))
    a = np.array([i for i, j, r in types], dtype=np.int64)
    b = 1 - a
    c = np.array([j * r for i, j, r in types], dtype=np.int64)
    d = np.array([(1 - j) * r for i, j, r in types], dtype=np.int64)
    target = np.outer(a, a) - np.outer(c, c) + np.outer(b, d) + np.outer(d, b)
    reconstructed = np.zeros((8, 8), dtype=np.int64)
    for weight, x, y in TERMS:
        x, y = np.array(x, dtype=np.int64), np.array(y, dtype=np.int64)
        assert np.all(np.abs(x) == 1) and np.all(np.abs(y) == 1)
        reconstructed += weight * (np.outer(x, y) + np.outer(y, x))
    assert np.array_equal(reconstructed, 24 * target)
    assert sum(term[0] for term in TERMS) == 32

    rng = np.random.default_rng(2026091902)
    checks = 0
    for n in range(2, 21):
        for _ in range(100):
            upper = np.triu(rng.integers(-5, 6, (n, n)), 1)
            matrix = upper + upper.T
            labels = rng.integers(0, 8, n)
            aa, bb, cc, dd = [z[labels] for z in (a, b, c, d)]
            twice_energy = int(aa @ matrix @ aa - cc @ matrix @ cc + 2 * bb @ matrix @ dd)
            right = sum(weight * int(np.array(x)[labels] @ matrix @ np.array(y)[labels])
                        for weight, x, y in TERMS)
            assert 12 * twice_energy == right
            checks += 1

    sharp = np.array([
        [-3, -1, -1, -1, -2, 2],
        [-1, 3, 2, 2, -1, 1],
        [-1, 2, -2, 1, 0, 0],
        [-1, 2, 1, -2, 0, 0],
        [-2, -1, 0, 0, 2, 1],
        [2, 1, 0, 0, 1, 2],
    ], dtype=np.int64)
    spins = np.array(list(itertools.product([-1, 1], repeat=6)), dtype=np.int64)
    sharp_beta = int(np.max(np.abs(spins @ sharp @ spins.T)))
    sharp_quadratic = np.einsum("bi,ij,bj->b", spins, sharp, spins)
    aa, bb, cc, dd = [z[:6] for z in (a, b, c, d)]
    sharp_value = int((aa @ sharp @ aa - cc @ sharp @ cc) // 2 + bb @ sharp @ dd)
    assert sharp_beta == 18 and sharp_value == 24
    assert int(sharp_quadratic.max()) == 18 and int(sharp_quadratic.min()) == -18

    output = {
        "theorem": "Delta(A) <= (4/3) beta(A), A real hollow symmetric",
        "normalization": "H_A(x)=x^T A x/2; beta(A)=max_(x,y signs)|x^T A y|",
        "types": types,
        "coefficient_matrix_K": target.tolist(),
        "certificate": [{"weight_numerator": weight, "weight_denominator": 24,
                         "left": x, "right": y} for weight, x, y in TERMS],
        "weight_sum": "4/3",
        "integer_matrix_identity": "sum w_r (x_r y_r^T+y_r x_r^T) = 24 K",
        "matrix_identity_pass": True,
        "random_integer_identity_checks": checks,
        "sharp_weighted_witness": {
            "matrix": sharp.tolist(), "beta": sharp_beta,
            "real_cube_quadratic_cap": 9,
            "diamond_vectors": [z.tolist() for z in (aa, bb, cc, dd)],
            "diamond_value": sharp_value, "ratio_to_beta": "4/3",
            "caution": "matrix has a nonzero diagonal; hollow blowups give asymptotic sharpness",
            "hollow_r_blowup_attained_diamond_value": "24*r^2-10*r",
            "hollow_r_blowup_beta_error_bound": "abs(beta-18*r^2)<=14*r",
        },
    }
    path = ROOT / "computations/results/twisted_chiral_diamond_bilinear_2026_09_19.json"
    path.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"matrix_identity_pass": True, "weight_sum": "4/3", "integer_checks": checks}))


if __name__ == "__main__":
    main()
