"""Exact finite spectral-witness constants and optional compression diagnostics.

No input or output files. Use --experiment for finite sampled-Walsh checks.
"""
import argparse
from fractions import Fraction as F
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
    q = multiply(p, p)
    q[0] += F(1, 10)
    mean = sum((x * semicircle(k) for k, x in enumerate(q)), F(0))
    energy = sum((x * semicircle(k + 1) for k, x in enumerate(q)), F(0))
    second = sum((x * semicircle(k) for k, x in enumerate(multiply(q, q))), F(0))
    assert mean == F(361, 50)
    assert energy == F(288, 25)
    lower = F(7, 22) * energy / mean
    assert lower == F(2016, 3971) and lower > F(507, 1000)
    return {'polynomial_low_to_high': list(map(str, p)), 'semicircle_R_mean': str(mean),
            'semicircle_LR_mean': str(energy), 'semicircle_R_squared_mean': str(second),
            'rational_cap_lower': str(lower), 'rational_cap_lower_display': float(lower)}


def sampled_walsh(n, depth, rng):
    assert n <= 1 << depth
    words = rng.choice(1 << depth, size=n, replace=False).astype(np.uint64)
    parity = np.bitwise_and(words[:, None], words[None, :])
    for shift in (32, 16, 8, 4, 2, 1):
        parity ^= parity >> np.uint64(shift)
    a = 1 - 2 * (parity & np.uint64(1)).astype(np.int64)
    np.fill_diagonal(a, 0)
    return a


def diagnostics(a):
    n = len(a)
    ell = a / math.sqrt(n)
    identity = np.eye(n)
    p = ell @ ell @ ell + 1.6 * (ell @ ell) - 0.4 * ell - 0.6 * identity
    q = 0.1 * identity + p @ p
    diagonal = np.diag(q)
    v = float(diagonal.mean())
    delta = float(np.sqrt(np.mean((diagonal - v) ** 2)))
    s = float(np.sum(q * q) / n)
    a_lower = float(diagonal.min())
    normalized = q / np.sqrt(diagonal[:, None] * diagonal[None, :])
    normalized = np.clip(normalized, -1, 1)
    linear_ratio = float(np.sum(ell * q) / (n * v))
    actual_linear = float(np.sum(ell * normalized) / n)
    gaussian = float(np.sum(ell * np.arcsin(normalized)) / (math.pi * n))
    bound = (linear_ratio - delta * math.sqrt(s) / a_lower ** 2
             - s / (a_lower ** 2 * math.sqrt(n))) / math.pi
    assert abs(actual_linear - linear_ratio) <= delta * math.sqrt(s) / a_lower ** 2 + 1e-10
    assert abs(gaussian - actual_linear / math.pi) <= s / (math.pi * a_lower ** 2 * math.sqrt(n)) + 1e-10
    assert gaussian >= bound - 1e-10
    return {'n': n, 'R_mean': v, 'R_diagonal_standard_deviation': delta,
            'R_second_trace_mean': s, 'R_minimum_diagonal': a_lower,
            'linear_trace_ratio': linear_ratio, 'finite_formula_lower': bound,
            'exact_Gaussian_expectation_float': gaussian}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--experiment', action='store_true')
    parser.add_argument('--depth', type=int, default=18)
    parser.add_argument('--orders', default='64,128,256,512')
    parser.add_argument('--seed', type=int, default=20260906)
    args = parser.parse_args()
    out = {'exact_constants': exact_constants()}
    if args.experiment:
        rng = np.random.default_rng(args.seed)
        out.update(seed=args.seed, depth=args.depth)
        out['finite_checks'] = [diagnostics(sampled_walsh(n, args.depth, rng))
                                for n in map(int, args.orders.split(','))]
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
