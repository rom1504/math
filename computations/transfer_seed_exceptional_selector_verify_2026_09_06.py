"""Finite covariance clipping and exceptional-selector defect checks.

No input or output files. Rational constants are exact; matrix checks
are finite floating-point diagnostics, not the proofs of the inequalities.
"""
from fractions import Fraction as F
import itertools
import json
import math

import numpy as np


def multiply(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def semicircle(k):
    return 0 if k % 2 else math.comb(k, k // 2) // (k // 2 + 1)


def exact_constants():
    p = [F(-3, 5), F(-2, 5), F(8, 5), F(1)]
    r = multiply(p, p)
    rr = multiply(r, r)
    v = sum(c * semicircle(k) for k, c in enumerate(r))
    e = sum(c * semicircle(k + 1) for k, c in enumerate(r))
    s = sum(c * semicircle(k) for k, c in enumerate(rr))
    assert (v, e, s) == (F(178, 25), F(288, 25), F(314747, 625))
    coeff = sum(map(abs, r))
    coeff2 = sum(map(abs, rr))
    assert (coeff, coeff2) == (F(276, 25), F(53896, 625))
    eta, b = F(1, 100000), F(1, 10)
    vp, vm, em, sp = v + coeff * eta, v - coeff * eta, e - coeff * eta, s + coeff2 * eta
    assert sp < 23 ** 2 and em - 2 * coeff * eta * 23 / b > 0
    finite = F(7, 22) * ((em - 2 * coeff * eta * 23 / b) / (vp + b)
                          - sp / ((vm + b) ** 2 * 1000))
    assert finite == F(145071245497770289701, 288781684033417350314)
    assert finite > F(50235, 100000) > F(1, 2)
    defect = F(58, 175) ** 2 / (8 * F(11, 7) * 23)
    assert defect == F(841, 2213750) > F(1, 3000)
    return {'P_squared_coefficients': list(map(str, r)),
            'P_fourth_even_coefficients': {str(k): str(c) for k, c in enumerate(rr) if k % 2 == 0},
            'semicircle_v_e_s': list(map(str, (v, e, s))),
            'coefficient_absolute_sums': list(map(str, (coeff, coeff2))),
            'finite_cap_lower': str(finite), 'finite_cap_lower_float': float(finite),
            'limiting_local_defect_lower': str(defect)}


def exact_cap(a):
    n = len(a)
    spins = np.array([(1,) + x for x in itertools.product((-1, 1), repeat=n - 1)])
    return float(np.max(np.abs(np.einsum('bi,ij,bj->b', spins, a, spins))) / (2 * n ** 1.5))


def check_clipping(a, r, b):
    n = len(a)
    ell = a / math.sqrt(n)
    v, e = np.trace(r) / n, np.sum(ell * r) / n
    s = np.sum(r * r) / n
    delta = np.std(np.diag(r))
    kept = np.flatnonzero(np.diag(r) <= v + b)
    removed = n - len(kept)
    assert removed / n <= delta ** 2 / b ** 2 + 1e-10
    ri, li = r[np.ix_(kept, kept)], ell[np.ix_(kept, kept)]
    linear = np.sum(li * ri) / n
    assert abs(e - linear) <= 2 * delta * math.sqrt(s) / b + 1e-9
    cov = ri / (v + b)
    if len(kept):
        np.fill_diagonal(cov, 1)
        assert np.linalg.eigvalsh(cov)[0] >= -1e-9
    gaussian = np.sum(li * np.arcsin(np.clip(cov, -1, 1))) / (math.pi * n)
    bound = ((e - 2 * delta * math.sqrt(s) / b) / (v + b)
             - s / ((v + b) ** 2 * math.sqrt(n))) / math.pi
    assert gaussian >= bound - 1e-9
    assert exact_cap(a) >= gaussian - 1e-9


def check_cubic_defect(a):
    n = len(a)
    ell = a / math.sqrt(n)
    powers = [np.eye(n)]
    for _ in range(12):
        powers.append(powers[-1] @ ell)
    mu = [np.trace(x) / n for x in powers]
    rs = []
    for sign in (1, -1):
        poly = sign * powers[3] + 1.6 * powers[2] - 0.4 * sign * ell - 0.6 * np.eye(n)
        rs.append(poly @ poly)
    vb = (np.trace(rs[0]) + np.trace(rs[1])) / (2 * n)
    eb = np.sum(ell * (rs[0] - rs[1])) / (2 * n)
    assert abs(vb - (9 / 25 - 44 / 25 * mu[2] + 44 / 25 * mu[4] + mu[6])) < 1e-8
    assert abs(eb - (12 / 25 * mu[2] - 62 / 25 * mu[4] + 16 / 5 * mu[6])) < 1e-8
    s = sum(np.sum(r * r) / n for r in rs) / 2
    delta = math.sqrt(sum(np.var(np.diag(r)) for r in rs) / 2)
    rooted = np.var(np.diag(44 / 25 * powers[4] + powers[6])) + np.var(np.diag(-62 / 25 * powers[3] + 16 / 5 * powers[5]))
    assert abs(delta ** 2 - rooted) < 1e-7
    c = exact_cap(a)
    gap = max(0, eb - math.pi * c * vb)
    necessary = max(0, gap ** 2 / (8 * math.pi * c * math.sqrt(s)) - math.sqrt(s) / (2 * math.sqrt(n)))
    assert delta >= necessary - 1e-9


def main():
    rng = np.random.default_rng(20260906)
    clipping = cubic = 0
    for n in range(2, 9):
        for trial in range(40):
            upper = np.triu(rng.choice((-1, 1), size=(n, n)), 1)
            a = upper + upper.T
            rank = 1 + trial % n
            factor = rng.normal(size=(n, rank)) * np.exp(rng.uniform(-2, 2, size=(n, 1)))
            if trial % 3 == 0:
                factor[0] = 0  # Singular covariance with an exact zero diagonal.
            r = factor @ factor.T
            for b in (0.05, 0.25, 1.0, 10.0):
                check_clipping(a, r, b)
                clipping += 1
            check_cubic_defect(a)
            cubic += 1
    print(json.dumps({'exact': exact_constants(), 'seed': 20260906,
                      'finite_clipping_checks': clipping,
                      'finite_two_sided_cubic_checks': cubic}, indent=2))


if __name__ == '__main__':
    main()
