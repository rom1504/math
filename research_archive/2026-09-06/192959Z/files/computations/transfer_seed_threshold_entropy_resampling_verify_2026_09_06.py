"""Exact rational threshold-entropy certificate and finite diagnostics.

No input or output files. Matrix diagnostics use floats and are not the proof.
"""
from fractions import Fraction as F
import itertools
import json
import math

import numpy as np


def rational_certificate():
    delta, c, s, r = F(1, 2 ** 24), F(33, 32), F(1, 32), F(1, 25)
    assert 5 / (F(78) * F(11, 4) ** 13) > delta
    assert 5 / (F(104) * F(8, 3) ** 14) < delta
    assert F(225, 29) / (1 - delta) ** 2 < 8
    expo = F(29, 25)
    degree = 12
    exp_upper = sum(expo ** k / math.factorial(k) for k in range(degree + 1))
    exp_upper += (expo ** (degree + 1) / math.factorial(degree + 1)
                  / (1 - expo / (degree + 2)))
    assert exp_upper < F(16, 5)
    assert F(55, 29) ** 2 * F(625, 624) < F(19, 10) ** 2
    lam = s * s / (2 * (r * r - s * s * c * c))
    assert lam == F(320000, 367951) < 1
    kap_upper = 32 * delta * r * F(19, 10)
    assert kap_upper < F(1, 10 ** 6)
    tiny = F(1, 10 ** 6)
    assert 8 * tiny / (1 - tiny) < 9 * tiny
    fb_upper = 1 / ((1 - tiny) ** 2 * (1 - 36 * tiny))
    assert fb_upper < F(1001, 1000)
    factor = 2 * F(19, 10) ** 2 * F(1001, 1000) / (1 - (s * c / r) ** 2)
    assert factor == F(23682154496, 1149846875) < 22
    xlo, xhi = F(25, 4), F(8)
    assert F(1, 8) - F(22, 2048) * xhi > 0
    new_credit = xlo / 8 - F(11, 2048) * xlo * xlo
    old_credit = F(33, 4) ** 2 / 128
    assert new_credit == F(18725, 32768)
    assert old_credit == F(1089, 2048)
    gap = new_credit - old_credit
    assert gap == F(1301, 32768) > 0
    p = F(31, 32)
    assert c * c * p > 1 and F(64) / p < F(33, 4) ** 2
    return {'delta': str(delta), 'C': str(c), 's': str(s), 'r': str(r),
            'determinant_lambda': str(lam), 'kappa_upper': str(kap_upper),
            'KL_factor_upper': str(factor), 'KL_factor_upper_float': float(factor),
            'new_credit_over_v_squared_lower': str(new_credit),
            'old_credit_over_v_squared_upper': str(old_credit),
            'credit_gap_over_v_squared_lower': str(gap),
            'credit_gap_lower_float_diagnostic': float(gap * (4 * delta * (1 - delta)) ** 2)}


def resampling_expectation(a, s, r, kappa):
    n = len(a)
    base = a / math.sqrt(n)
    c = float(np.linalg.norm(base, 2))
    lam = s * s / (2 * (r * r - s * s * c * c))
    assert abs(s) * c < r < 1
    determinant_mean = count_mean = 0.0
    for bits in itertools.product((0, 1), repeat=n):
        kept = np.flatnonzero(bits)
        k = len(kept)
        mass = kappa ** k * (1 - kappa) ** (n - k)
        restricted = base[np.ix_(kept, kept)]
        eigen = np.linalg.eigvalsh(restricted)
        value = math.exp(-0.5 * float(np.log1p(-(s / r) ** 2 * eigen ** 2).sum()))
        bound = math.exp(lam * k * (k - 1) / n)
        assert value <= bound + 1e-11
        determinant_mean += mass * value
        count_mean += mass * bound
    fb = math.exp(4 * lam * kappa * math.exp(2 * lam) / (1 - kappa)) / (1 - kappa) ** 2
    assert math.log(count_mean) <= n * lam * kappa * kappa * fb + math.log(n + 1) + 1e-11
    return determinant_mean, count_mean


def finite_checks():
    rng = np.random.default_rng(20260906)
    r, s = 0.5, 0.07
    kappa = 2 / math.pi * math.asin(r)  # Exact threshold formula at t=0.
    small_threshold = general = 0
    for n in (2, 3):
        edges = list(itertools.combinations(range(n), 2))
        for signs in itertools.product((-1, 1), repeat=len(edges)):
            a = np.zeros((n, n))
            for (i, j), value in zip(edges, signs):
                a[i, j] = a[j, i] = value
            determinant, count = resampling_expectation(a, s, r, kappa)
            # At n<=3 and t=0, only the exact pair sign moments survive.
            chi2 = sum((2 / math.pi * math.asin(s * a[i, j] / math.sqrt(n))) ** 2
                       for i, j in edges)
            assert 1 + chi2 <= determinant + 1e-11 <= count + 2e-11
            small_threshold += 1
    for n in range(2, 9):
        for _ in range(8):
            upper = np.triu(rng.choice((-1, 1), size=(n, n)), 1)
            a = upper + upper.T
            resampling_expectation(a, s, r, kappa)
            general += 1
    return {'exact_formula_threshold_instances_float_checked': small_threshold,
            'general_flat_determinant_instances_float_checked': general,
            'seed': 20260906}


if __name__ == '__main__':
    print(json.dumps({'rational_certificate': rational_certificate(),
                      'finite_diagnostics': finite_checks()}, indent=2))
