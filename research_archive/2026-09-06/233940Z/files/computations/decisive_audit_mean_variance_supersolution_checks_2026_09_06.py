"""Finite-channel checks of the precision Schur supersolution proof.

No optimizer or asymptotic assertion is inferred from these diagnostics.
Run from the repository with .venv/bin/python and this path.
"""

from fractions import Fraction
import math
import numpy as np


def entropy(p):
    a = np.asarray(p, dtype=float).ravel()
    a = a[a > 0]
    return float(-np.dot(a, np.log(a)))


def marginal(p, axes):
    discarded = tuple(i for i in range(p.ndim) if i not in axes)
    return p.sum(axis=discarded) if discarded else p


def mutual(p, x_axes, y_axes):
    return (entropy(marginal(p, x_axes)) + entropy(marginal(p, y_axes))
            - entropy(marginal(p, tuple(sorted(set(x_axes + y_axes))))))


def variance_cost(p, values, value_axis, label_axes):
    kept = tuple(sorted((value_axis,) + label_axes))
    joint = marginal(p, kept)
    axis = kept.index(value_axis)
    shape = [1] * joint.ndim
    shape[axis] = len(values)
    vals = np.asarray(values).reshape(shape)
    mass = joint.sum(axis=axis)
    first = (joint * vals).sum(axis=axis)
    second = (joint * vals * vals).sum(axis=axis)
    ratio = np.divide(first * first, mass, out=np.zeros_like(mass), where=mass > 0)
    return float((second - ratio).sum())


def main():
    rng = np.random.default_rng(2026090623)
    exact_cases = 0
    for _ in range(400):
        l1, l2 = (Fraction(int(x), 17) for x in rng.integers(1, 100, size=2))
        ar = (l1 + l2) / 2
        hr = l1 * l2 / ar
        dr = (l1 - l2) / (l1 + l2)
        x, y = (Fraction(int(x), 7) for x in rng.integers(-30, 31, size=2))
        # Rotation square identity, without introducing irrational sqrt(2).
        assert (l1 * (x + y) ** 2 + l2 * (x - y) ** 2) / 2 == ar * (x + dr * y) ** 2 + hr * y * y
        exact_cases += 1

    minima = [math.inf, math.inf]
    max_info_error = 0.0
    for case in range(500):
        na, nb = (int(x) for x in rng.integers(1, 6, size=2))
        va = np.arange(na, dtype=float) - na / 2
        vb = np.arange(nb, dtype=float) - nb / 2
        base = rng.random((na, nb))
        if case % 7 == 0:
            base *= rng.random((na, nb)) > 0.5
            if base.sum() == 0:
                base[0, 0] = 1
        base /= base.sum()
        u = (va[:, None] + vb[None, :]) / math.sqrt(2)
        v = (va[:, None] - vb[None, :]) / math.sqrt(2)
        us, ui = np.unique(u, return_inverse=True)
        vs, vi = np.unique(v, return_inverse=True)
        ui, vi = ui.reshape(na, nb), vi.reshape(na, nb)
        cm = rng.random((len(us), 3))
        cn = rng.random((len(vs), 4))
        cm /= cm.sum(axis=1, keepdims=True)
        cn /= cn.sum(axis=1, keepdims=True)
        p = base[:, :, None, None] * cm[ui, :, None] * cn[vi, None, :]

        pu = np.zeros((len(us), 3))
        pv = np.zeros((len(vs), 4))
        for i in range(na):
            for j in range(nb):
                pu[ui[i, j]] += base[i, j] * cm[ui[i, j]]
                pv[vi[i, j]] += base[i, j] * cn[vi[i, j]]
        ium, ivn = mutual(pu, (0,), (1,)), mutual(pv, (0,), (1,))
        lhs_info = mutual(p, (0,), (1, 2, 3)) + mutual(p, (1,), (2, 3))
        rhs_info = mutual(p, (0,), (1,)) + ium + ivn - mutual(p, (2,), (3,))
        max_info_error = max(max_info_error, abs(lhs_info - rhs_info))
        assert abs(lhs_info - rhs_info) < 1e-11

        t = float(np.exp(rng.uniform(-2, 2)))
        l1, l2 = t * np.exp(rng.uniform(-10, 0, size=2))
        ar, hr = (l1 + l2) / 2, 2 * l1 * l2 / (l1 + l2)
        parent_cost = ar * variance_cost(p, va, 0, (1, 2, 3)) + hr * variance_cost(p, vb, 1, (2, 3))
        child_cost = l1 * variance_cost(pu, us, 0, (1,)) + l2 * variance_cost(pv, vs, 0, (1,))
        slack = child_cost - parent_cost
        assert slack >= -1e-10
        minima[0] = min(minima[0], slack)
        c = lambda z: math.log(z * (2 * t - z) / (t * t)) / 4
        curvature = c(ar) + c(hr) - c(l1) - c(l2)
        assert curvature >= -1e-12
        minima[1] = min(minima[1], curvature)
    print({"status": "PASS", "exact_rational_factorizations": exact_cases,
           "finite_channel_cases": 500, "minimum_cost_slack": minima[0],
           "minimum_precision_slack": minima[1], "maximum_information_error": max_info_error})


if __name__ == "__main__":
    main()
