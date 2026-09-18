#!/usr/bin/env python3
"""Finite replays for scoped Boolean BH counterexample barriers.

Integer granularity/support checks are exact. Floating comparisons of
complex norms verify the implementation, not the analytic theorems.
"""
import itertools
import json
import math
import numpy as np


def fwht(a):
    a = np.array(a, copy=True)
    width = 1
    while width < len(a):
        for start in range(0, len(a), 2 * width):
            left = a[start:start + width].copy()
            right = a[start + width:start + 2 * width].copy()
            a[start:start + width] = left + right
            a[start + width:start + 2 * width] = left - right
        width *= 2
    return a


def weight(mask):
    return bin(int(mask)).count("1")


def degree(a, tol=1e-11):
    return max([weight(i) for i, v in enumerate(a) if abs(v) > tol] + [0])


def ratio(a, m=None):
    if m is None:
        m = degree(a)
    if m == 0:
        return 1.0
    q = 2.0 * m / (m + 1)
    return float(np.sum(np.abs(a) ** q) ** (1.0 / q) /
                 np.max(np.abs(fwht(a))))


def rank_gf2(rows):
    basis = {}
    for row in rows:
        row = int(row)
        while row:
            pivot = row.bit_length() - 1
            if pivot in basis:
                row ^= basis[pivot]
            else:
                basis[pivot] = row
                break
    return len(basis)


def mul_char(a, mask, c, s):
    return c * a + 1j * s * a[np.arange(len(a)) ^ mask]


def main():
    rng = np.random.RandomState(20260918)
    counts = {}
    # Exact sixth-root extremizer in the Eisenstein basis (1,zeta).
    root_pairs = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
    masks6 = [0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 16, 17, 18, 20, 24]
    phases6 = [0, 2, 0, 5, 1, 4, 2, 4, 3, 1, 4, 1, 2, 2, 3, 2]
    assert masks6 == [s for s in range(32) if weight(s) <= 2]
    sixth_histogram = {}
    for x in range(32):
        ar = sum((-1) ** weight(x & s) * root_pairs[k][0]
                 for s, k in zip(masks6, phases6))
        ai = sum((-1) ** weight(x & s) * root_pairs[k][1]
                 for s, k in zip(masks6, phases6))
        assert ar % 4 == 0 and ai % 4 == 0
        assert (ar // 4, ai // 4) in root_pairs
        assert ar * ar + ar * ai + ai * ai == 16
        key = str((ar, ai))
        sixth_histogram[key] = sixth_histogram.get(key, 0) + 1
    assert len(masks6) ** 3 == 2 ** 4 * 4 ** 4  # fourth power of R_2=2
    counts["sixth_root_extremizer_exact"] = 32
    # Real Boolean granularity at every actual degree: complete n=4 census.
    for pattern in range(1 << 16):
        values = np.array([1 if (pattern >> i) & 1 else -1
                           for i in range(16)], dtype=np.int64)
        raw = fwht(values)
        m = degree(raw, 0)
        if m:
            assert np.all(raw % (2 ** (5 - m)) == 0)
            assert np.count_nonzero(raw) <= 4 ** (m - 1)
    counts["boolean_n4_exact"] = 1 << 16

    # Fourth-root-valued functions: complete n=3 census, integer pairs.
    phases = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    for digits in itertools.product(range(4), repeat=8):
        re = fwht(np.array([phases[j][0] for j in digits], dtype=np.int64))
        im = fwht(np.array([phases[j][1] for j in digits], dtype=np.int64))
        supp = np.nonzero((re != 0) | (im != 0))[0]
        m = max([weight(j) for j in supp] + [0])
        assert np.all(re % (2 ** (3 - m)) == 0)
        assert np.all(im % (2 ** (3 - m)) == 0)
        assert len(supp) <= 4 ** m
    counts["four_phase_n3_exact"] = 4 ** 8

    # Generic coded phase products: all reachable characters survive.
    worst_generic = 0.0
    for n in range(1, 9):
        for trial in range(30):
            masks = rng.randint(1, 1 << n, size=n + 3)
            angles = rng.uniform(.12, .63, size=len(masks))
            a = np.zeros(1 << n, dtype=complex)
            a[0] = 1
            reachable = {0}
            for mask, theta in zip(masks, angles):
                a = mul_char(a, int(mask), math.cos(theta), math.sin(theta))
                reachable |= {v ^ int(mask) for v in list(reachable)}
            assert set(np.nonzero(np.abs(a) > 1e-11)[0]) == reachable
            m = degree(a)
            r = rank_gf2(masks)
            assert r <= m
            assert np.max(np.abs(np.abs(fwht(a)) - 1)) < 1e-10
            rr = ratio(a, m)
            assert rr <= math.sqrt(2) + 1e-10
            worst_generic = max(worst_generic, rr)
    counts["generic_coded_products"] = 240

    # Changing-seed tensors and arbitrary shared-variable selector branches.
    worst_selector = 0.0
    for trial in range(200):
        factors = []
        budgets = []
        for j in range(2 + trial % 2):
            n = 2 + j % 2
            d = 1 + (trial + j) % n
            a = rng.normal(size=1 << n) + 1j * rng.normal(size=1 << n)
            a[[i for i in range(1 << n) if weight(i) > d]] = 0
            a /= np.max(np.abs(fwht(a)))
            factors.append(a)
            budgets.append(d)
        tensor = factors[0]
        for a in factors[1:]:
            tensor = np.kron(a, tensor)
        D = sum(budgets)
        geometric = math.exp(sum(d * math.log(ratio(a, d))
                                for a, d in zip(factors, budgets)) / D)
        assert ratio(tensor, D) <= geometric + 2e-10
        g, h = factors[:2]
        size = max(len(g), len(h))
        g = np.pad(g, (0, size - len(g)))
        h = np.pad(h, (0, size - len(h)))
        selected = np.concatenate(((g + h) / 2, (g - h) / 2))
        Dsel = 1 + max(budgets[:2])
        C = max(1, ratio(g, budgets[0]), ratio(h, budgets[1]))
        bound = 2 ** (1.0 / Dsel) * C ** ((Dsel - 1.0) / Dsel)
        observed = ratio(selected, Dsel)
        assert observed <= bound + 2e-10
        worst_selector = max(worst_selector, observed)
    counts["changing_seed_tensors"] = 200
    counts["shared_variable_selectors"] = 200

    # Exact entropy identity, including a non-unimodular seed's decay.
    tensor_limits = []
    for values in [np.array([1, 1j, -1, -1j]),
                   np.array([1, .2j, -.4, .7j])]:
        a = fwht(values) / 4
        d = degree(a)
        v = float(np.sum(np.abs(a) ** 2))
        pi = np.abs(a[np.abs(a) > 1e-13]) ** 2 / v
        rows = []
        for k in [1, 2, 3, 4]:
            aa = a.copy()
            for unused in range(k - 1):
                aa = np.kron(aa, a)
            alpha = (k * d) / (k * d + 1.0)
            H = math.log(float(np.sum(pi ** alpha))) / (1 - alpha)
            exact = v ** (k / 2.0) * math.exp(H / (2 * d))
            assert abs(ratio(aa, k * d) - exact) < 1e-9
            rows.append([k, exact])
        tensor_limits.append(rows)
    counts["tensor_entropy_identities"] = 8
    print(json.dumps({"status": "PASS", "counts": counts,
                      "sixth_root_histogram": sixth_histogram,
                      "worst_generic_ratio": worst_generic,
                      "worst_selector_ratio": worst_selector,
                      "tensor_examples": tensor_limits}, indent=2))


if __name__ == "__main__":
    main()
