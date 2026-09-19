#!/usr/bin/env python3
"""Independent finite checks for the final proof audit; no optimization solver."""

import itertools
import json
import random
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def multiply(a, x):
    return [sum(ai[j]*x[j] for j in range(len(x))) for ai in a]


def bilinear(a, x, y):
    return sum(xi*yi for xi, yi in zip(x, multiply(a, y)))


def main():
    rng = random.Random(202609192030)
    contraction_checks = 0
    for _ in range(10000):
        h = [rng.randrange(-10, 11) for _ in range(8)]
        z = [rng.randrange(-10, 11) for _ in range(8)]
        left = max(a+abs(b) for a, b in zip(h, z)) + max(a-abs(b) for a, b in zip(h, z))
        right = max(a+b for a, b in zip(h, z)) + max(a-b for a, b in zip(h, z))
        assert left <= right
        contraction_checks += 1

    cycle_checks = 0
    for n in range(4, 21):
        for _ in range(100):
            p = list(range(n))
            rng.shuffle(p)
            r = rng.randrange(1, n//4 + 1)
            removed = set(rng.sample(range(n), r))
            repaired = []
            for i in range(n):
                if i in removed:
                    repaired.append(i)
                else:
                    j = p[i]
                    while j in removed:
                        j = p[j]
                    repaired.append(j)
            assert sorted(repaired) == list(range(n))
            changed = {i for i in range(n) if p[i] != repaired[i]}
            union = removed | {p[i] for i in removed}
            assert changed <= removed | {i for i in range(n) if p[i] in removed}
            assert {p[i] for i in changed} <= union
            assert {repaired[i] for i in changed} <= union
            assert {repaired[i] for i in removed} == removed
            cycle_checks += 1

    insertion_checks = 0
    for n in range(1, 9):
        for r in range(1, 6):
            size = n+r
            for _ in range(20):
                p = list(range(n))
                rng.shuffle(p)
                signs = [rng.choice([-1, 1]) for _ in range(n)]
                error = [[0]*size for _ in range(size)]
                for i in range(size):
                    for j in range(i+1, size):
                        if j >= n:
                            error[i][j] = error[j][i] = rng.choice([-1, 1])
                x = [rng.choice([-1, 1]) for _ in range(size)]
                y = [rng.choice([-1, 1]) for _ in range(size)]
                # G e_i=s_i e_p(i), reconstructed without array operations.
                gx, gy = x[n:].copy(), y[n:].copy()
                gx = [0]*n + gx
                gy = [0]*n + gy
                for i in range(n):
                    gx[p[i]], gy[p[i]] = signs[i]*x[i], signs[i]*y[i]
                direct = (bilinear(error, x, x)-bilinear(error, y, y))//2 + bilinear(error, gx, gy)
                expanded = square_sum = 0
                for i in range(size):
                    for j in range(i+1, size):
                        if j < n:
                            continue
                        if i < n:
                            coefficient = x[j]*(x[i]+gy[i])+y[j]*(gx[i]-y[i])
                        else:
                            coefficient = x[i]*x[j]-y[i]*y[j]+x[i]*y[j]+x[j]*y[i]
                            assert abs(coefficient) == 2
                        expanded += error[i][j]*coefficient
                        square_sum += coefficient**2
                assert direct == expanded
                assert square_sum <= 8*n*r + 2*r*(r-1)
                insertion_checks += 1

    matrix = [
        [-3, -1, -1, -1, -2, 2],
        [-1, 3, 2, 2, -1, 1],
        [-1, 2, -2, 1, 0, 0],
        [-1, 2, 1, -2, 0, 0],
        [-2, -1, 0, 0, 2, 1],
        [2, 1, 0, 0, 1, 2],
    ]
    spins = list(itertools.product([-1, 1], repeat=6))
    beta = max(abs(bilinear(matrix, x, y)) for x in spins for y in spins)
    extrema = [bilinear(matrix, x, x)//2 for x in spins]
    a, b = [0, 0, 0, 0, 1, 1], [1, 1, 1, 1, 0, 0]
    c, d = [0, 0, -1, 1, 0, 0], [-1, 1, 0, 0, -1, 1]
    diamond = (bilinear(matrix, a, a)-bilinear(matrix, c, c))//2 + bilinear(matrix, b, d)
    diagonal_twice = sum(matrix[i][i]*(a[i]**2-c[i]**2+2*b[i]*d[i]) for i in range(6))
    diagonal_mass = sum(abs(matrix[i][i]) for i in range(6))
    assert (beta, min(extrema), max(extrema), diamond, diagonal_twice, diagonal_mass) == (18, -9, 9, 24, 20, 14)

    result = {
        "all_checks_pass": True,
        "factor_one_contraction_integer_checks": contraction_checks,
        "cycle_skipping_image_containment_checks": cycle_checks,
        "insertion_coefficient_expansion_checks": insertion_checks,
        "weighted_sharpness_plain_scalar_checks": {
            "bilinear_spin_pairs": 4096, "beta": beta,
            "quadratic_extrema": [min(extrema), max(extrema)],
            "diamond_witness_value": diamond,
            "hollow_blowup_witness_diagonal_loss_per_clone": diagonal_twice//2,
            "hollow_blowup_diagonal_mass_per_clone": diagonal_mass,
        },
    }
    (ROOT / "computations/results/twisted_chiral_uniform_final_audit_2026_09_19.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
