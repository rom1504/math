#!/usr/bin/env python3
"""Exact full-sign Hall profiles and scaled-permanent inequality checks."""

import json
import math
import random
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def matvec(a, x):
    return [sum(ai[j]*x[j] for j in range(len(x))) for ai in a]


def permanent(a):
    size = len(a)
    values = {0: Fraction(1)}
    for i in range(size):
        new = {}
        for mask, value in values.items():
            for j in range(size):
                if not mask & (1 << j):
                    target = mask | (1 << j)
                    new[target] = new.get(target, Fraction(0)) + value*a[i][j]
        values = new
    return values[(1 << size)-1]


def tensor_vector(x, y):
    return [a*b for a in x for b in y]


def tensor_matrix(x, y):
    return [[x[i][j]*y[a][b] for j in range(len(x)) for b in range(len(y))]
            for i in range(len(x)) for a in range(len(y))]


def main():
    records = []
    for m in [24, 28, 32, 40, 48, 64]:
        size = 2*m
        a = [[0 if i == j else (-1 if i >= m and j >= m else 1)
              for j in range(size)] for i in range(size)]
        x, y = [1]*size, [1]*m + [-1]*m
        ax, ay = matvec(a, x), matvec(a, y)
        alpha = [u*v for u, v in zip(x, ax)]
        gamma = [u*v for u, v in zip(y, ay)]
        cap = max(abs((u*u-v*v)//2+u*v)
                  for u in range(-m, m+1, 2) for v in range(-m, m+1, 2))
        assert cap == m*m
        for h in sorted({1, m//4, m//2, m-2}):
            z = [-1]*h + [1]*(size-h)
            w = z[:m] + [-1]*m
            az, aw = matvec(a, z), matvec(a, w)
            b = [u*v for u, v in zip(z, aw)]
            c = [u*v for u, v in zip(w, az)]
            margins_plus = [[min(alpha[i]+b[j], -gamma[i]+c[j])
                             for j in range(m)] for i in range(m)]
            margins_minus = [[min(alpha[i]+b[j], -gamma[i]+c[j])
                              for j in range(m, size)] for i in range(m, size)]
            assert all(row == [-2*m+2*h]*h + [2*m-2*h-2]*(m-h)
                       for row in margins_plus)
            assert all(row == [2*m-2*h]*m for row in margins_minus)
            assert 2*m-2*h-2 >= 2
            energy = (sum(u*v for u, v in zip(x, ax))-sum(u*v for u, v in zip(y, ay)))//2 + sum(b)
            assert energy == 4*m*m-4*m*h+4*h*h-2*m
            assert energy-2*m > 0 and (energy-2*m)**2 > 8*cap**2
            records.append({
                "half_order": m, "defect_count": h, "child_cap": cap,
                "base_energy": energy,
                "minimum_matching_energy": energy-2*m,
                "all_matching_energies_above_2sqrt2_cap_exact_square_check": True,
                "plus_row_degrees": [m-h]*m, "minus_row_degrees": [m]*m,
                "plus_sector_Hall_deficit": h,
                "exact_stable_conditional_probability": 0,
                "row_bound_symbolic": f"(({m-h}!)^(1/{m-h}))^{m}/{m}!",
                "log_row_bound_diagnostic": m*math.lgamma(m-h+1)/(m-h)-math.lgamma(m+1),
            })

    rng = random.Random(202609192041)
    scaled_checks = 0
    for size in range(1, 8):
        for k in range(1, size+1):
            for _ in range(10):
                matrices = []
                for __ in range(2):
                    rows, cols = list(range(size)), list(range(size))
                    rng.shuffle(rows); rng.shuffle(cols)
                    item = [[0]*size for i in range(size)]
                    for i in range(size):
                        for shift in range(k):
                            item[rows[i]][cols[(i+shift) % size]] = 1
                    matrices.append(item)
                b = [[Fraction(matrices[0][i][j]+matrices[1][i][j], 2*k)
                      for j in range(size)] for i in range(size)]
                assert all(sum(row) == 1 for row in b)
                assert all(sum(b[i][j] for i in range(size)) == 1 for j in range(size))
                assert all(v <= Fraction(1, k) for row in b for v in row)
                per_b = permanent(b)
                assert per_b**k * k**(size*k) <= math.factorial(k)**size
                aa = [rng.randrange(1, 5) for i in range(size)]
                bb = [rng.randrange(1, 5) for i in range(size)]
                original = [[b[i][j]/(aa[i]*bb[j]) for j in range(size)] for i in range(size)]
                factor = math.prod(aa)*math.prod(bb)
                assert permanent(original)*factor == per_b
                scaled_checks += 1
    h4 = [[1,1,1,1],[1,-1,1,-1],[1,1,-1,-1],[1,-1,-1,1]]
    x4, y4 = [1,1,1,-1], [-1,1,1,1]
    z4, w4 = [1,1,-1,1], [1,-1,1,1]
    hmat, u = [[1]], [1]
    sylvester = []
    for exponent in range(1, 4):
        hmat = tensor_matrix(h4, hmat)
        x, y, z, w = [tensor_vector(v, u) for v in [x4, y4, z4, w4]]
        size, root = len(hmat), 2**exponent
        delta = [hmat[i][i] for i in range(size)]
        assert sum(delta) == 0
        assert matvec(hmat, x) == [root*v for v in x]
        assert matvec(hmat, y) == [-root*v for v in y]
        assert matvec(hmat, z) == [root*v for v in w]
        assert matvec(hmat, w) == [root*v for v in z]
        a = [[hmat[i][j] if i != j else 0 for j in range(size)] for i in range(size)]
        t, relative = [x[i]*y[i] for i in range(size)], [z[i]*w[i] for i in range(size)]
        assert sum(t) == sum(relative) == 0
        ax, ay, az, aw = [matvec(a, v) for v in [x, y, z, w]]
        alpha = [x[i]*ax[i] for i in range(size)]
        gamma = [y[i]*ay[i] for i in range(size)]
        b, c = [z[i]*aw[i] for i in range(size)], [w[i]*az[i] for i in range(size)]
        assert alpha == [root-d for d in delta]
        assert [-v for v in gamma] == [root+d for d in delta]
        assert b == c == [root-delta[i]*relative[i] for i in range(size)]
        margin = min(min(alpha[i]+b[j], -gamma[i]+c[j])
                     for i in range(size) for j in range(size) if t[i] == relative[j])
        assert margin >= 2*root-2 >= 2
        base = (sum(x[i]*ax[i]-y[i]*ay[i] for i in range(size)))//2 + sum(b)
        assert base == 2*size*root-sum(delta[i]*relative[i] for i in range(size))
        p = [None]*size
        for sign in [-1, 1]:
            old = [i for i in range(size) if t[i] == sign]
            new = [j for j in range(size) if relative[j] == sign]
            for i, j in zip(old, new):
                p[i] = j
        switches = [z[p[i]]*x[i] for i in range(size)]
        bridge = [[switches[i]*switches[j]*a[p[i]][p[j]] for j in range(size)] for i in range(size)]
        # u_i=-1 is the worst possible matching energy; d_i=-t_i.
        for i in range(size):
            bridge[i][i] = -t[i]
        by, bx = matvec(bridge, y), matvec(bridge, x)
        fields = [x[i]*(ax[i]+by[i]) for i in range(size)] + [y[i]*(-ay[i]+bx[i]) for i in range(size)]
        assert min(fields) >= 1
        energy = sum(fields)//2
        assert energy == base-size
        cap = size*root//2
        if size >= 16:
            assert energy > 0 and energy**2 > 8*cap**2
        sylvester.append({"order": size, "quadratic_cap": cap,
                          "sector_sizes": [size//2, size//2],
                          "minimum_core_margin": margin, "base_energy": base,
                          "worst_matching_energy": energy,
                          "worst_matching_minimum_local_field": min(fields),
                          "exact_conditional_stability_probability": 1,
                          "exact_conditional_high_energy_probability_at_target": 1 if size >= 16 else None})
        u = tensor_vector(x4, u)

    output = {
        "all_checks_pass": True, "Hall_profile_checks": records,
        "scaled_permanent_rational_checks": scaled_checks,
        "low_cap_complete_support_calibration": sylvester,
        "scope": "conditional row-only false positives; no lower bound for the averaged selector",
    }
    (ROOT / "computations/results/twisted_chiral_stability_row_bound_audit_2026_09_19.json").write_text(json.dumps(output, indent=2)+"\n")
    print(json.dumps({"all_checks_pass": True, "Hall_profiles": len(records), "scaled_permanent_checks": scaled_checks,
                      "low_cap_calibration": sylvester}))


if __name__ == "__main__":
    main()
