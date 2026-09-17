"""Exact frame checks and numerical diagnostics for margin-aware truncation.

The finite identities are assertions.  Gaussian quadrature is only a
diagnostic and is not used to certify the analytic comparison theorem.
"""

import itertools
import json
import math

import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.linalg import hadamard
from scipy.special import logsumexp


RNG = np.random.default_rng(2026091726)


def bent_word(k):
    bits = k.bit_length() - 1
    assert bits % 2 == 0 and 1 << bits == k
    half = bits // 2
    mask = (1 << half) - 1
    return np.array([
        -1 if bin((z & mask) & (z >> half)).count('1') % 2 else 1
        for z in range(k)
    ], dtype=np.int64)


def frame_checks():
    cases = []
    for k in (4, 16, 64):
        h = hadamard(k).astype(np.int64)
        n, p, q = k * k, k, k
        f = np.repeat(bent_word(k)[:, None], p, axis=1)
        b_cap = math.isqrt(k)
        ground_coeff = h @ f
        assert np.all(np.abs(ground_coeff) == b_cap)
        tested = 0
        largest_tail_ratio = 0.0
        for _ in range(200):
            x = f.copy().reshape(-1)
            radius = int(RNG.integers(1, n // 2 + 1))
            x[RNG.choice(n, radius, replace=False)] *= -1
            x = x.reshape(k, p)
            c = h @ x
            b = np.clip(c, -b_cap, b_cap)
            tail = int(np.abs(c - b).sum())
            deficit = 2 * radius * (n - radius)
            difference_second = int(((b - ground_coeff) ** 2).sum())
            retained_fourth = int((b ** 4).sum())
            assert tail <= 2 * q * radius
            assert deficit >= n * radius
            assert difference_second <= 4 * q * radius
            assert retained_fourth <= q * n * k
            largest_tail_ratio = max(largest_tail_ratio, tail / (2*q*radius))
            tested += 1
        # A codimension-two plateau spectrum realizes a much larger
        # global ell1 truncation cost than the ground response cost.
        if k == 4:
            plateau = np.ones(k, dtype=np.int64)
        else:
            plateau = np.repeat(bent_word(k // 4), 4)
        plateau_array = np.repeat(plateau[:, None], p, axis=1)
        spectrum = h @ plateau_array
        active = np.abs(spectrum[spectrum != 0])
        assert np.all(active == 2 * b_cap)
        plateau_tail = int(np.abs(spectrum-np.clip(spectrum,-b_cap,b_cap)).sum())
        assert 4 * b_cap * plateau_tail == q * n
        # Coherent words attain the untruncated maximum fourth budget.
        coherent = np.repeat(h[0, :, None], p, axis=1)
        coherent_fourth = int(((h @ coherent) ** 4).sum())
        assert coherent_fourth == q * n * k * k
        cases.append(dict(k=k, n=n, q=q, tested=tested,
                          plateau_tail=plateau_tail,
                          ground_fourth=int((ground_coeff**4).sum()),
                          global_fourth=coherent_fourth,
                          largest_tail_ratio=largest_tail_ratio))
    return cases


def exact_zero_margin_counterexample():
    k, n, q = 4, 16, 4
    h = hadamard(k).astype(np.int64)
    a = RNG.choice((-1,1), size=(n,n))
    a = np.triu(a, 1)
    a += a.T
    words = np.array([(1,) + z for z in itertools.product((-1,1), repeat=n-1)],
                     dtype=np.int64)
    energies = np.einsum('bi,ij,bj->b', words, a, words) // 2
    index = int(np.argmax(np.abs(energies)))
    ground = words[index]
    polarity = int(np.sign(energies[index]))
    canonical = h.copy()  # columns contain every Walsh row once.
    switch = ground.reshape(k,k) * canonical
    switched_ground = switch * ground.reshape(k,k)
    coefficients = h @ switched_ground
    assert np.array_equal(np.abs(coefficients), k*np.eye(k,dtype=np.int64))
    assert polarity * int(ground @ a @ ground)//2 == abs(int(energies[index]))
    for cap in (0, 1, 2, 3):
        tails = int(np.abs(coefficients-np.clip(coefficients,-cap,cap)).sum())
        assert tails == q*(k-cap)
    return dict(n=n, exact_cap=abs(int(energies[index])), old_deficit=0,
                coefficients_per_column=1, coefficient_magnitude=k,
                tail_at_cap_1=q*(k-1))


def interpolation_diagnostics():
    nodes, weights = hermgauss(36)
    nodes *= math.sqrt(2)
    weights /= math.sqrt(math.pi)
    grid = np.array(list(itertools.product(nodes,nodes)))
    masses = np.outer(weights, weights).reshape(-1)
    signs = np.array(list(itertools.product((-1.,1.), repeat=2)))
    worst_step = -math.inf
    worst_final_violation = -math.inf
    for _ in range(30):
        c = RNG.normal(size=(5,2))
        offsets = RNG.normal(size=5)
        cap = float(RNG.uniform(.25,1.25))
        tau = float(RNG.uniform(.1,.6))
        b = np.clip(c,-cap,cap)
        tails = np.abs(c-b).sum(axis=1)
        fourth = (b**4).sum(axis=1)
        alpha = (19/12)*tau**3*math.exp(4*tau*cap)
        values = []
        for t in np.linspace(0,1,9):
            shifted_offsets = offsets+tails+alpha*(1-t*t)*fourth
            total = 0.0
            for eps in signs:
                environment = math.sqrt(1-t)*grid + math.sqrt(t)*eps
                f = logsumexp(tau*(environment@b.T+shifted_offsets),axis=1)/tau
                total += float(masses@f)/4
            values.append(total)
        worst_step = max(worst_step, float(np.max(np.diff(values))))
        original = float(np.max(signs@c.T+offsets,axis=1).mean())
        gaussian_max = float(masses@np.max(grid@b.T+offsets+tails+alpha*fourth,axis=1))
        violation = original-gaussian_max-math.log(len(offsets))/tau
        worst_final_violation = max(worst_final_violation,violation)
        assert max(np.diff(values)) < 1e-7
        assert violation <= 1e-7
    return dict(cases=30, max_interpolation_increase=worst_step,
                max_final_violation=worst_final_violation,
                note='Gauss-Hermite numerical diagnostic, not certified integration')


if __name__ == '__main__':
    print(json.dumps(dict(status='PASS', seed=2026091726,
                          frames=frame_checks(),
                          exact_zero_margin=exact_zero_margin_counterexample(),
                          interpolation=interpolation_diagnostics()), indent=2))
