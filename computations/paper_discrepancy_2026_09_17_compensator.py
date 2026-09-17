"""Exact rational counterexample to averaging the GS phase compensator.

Every branching probability, projector, and exponential-moment coefficient
is computed over the rationals. The final falsification uses the single
(-1,-1,-1) leaf and a rational lower bound on exp(8), not floating-point
evaluation of an exponential. Floating point is printed only for context.
"""

from __future__ import annotations

import json

import sympy as sp


def projection(v, active):
    if not active:
        return sp.zeros(v.rows)
    columns = v[:, active]
    return columns * (columns.T * columns).inv() * columns.T


def exact_walk(v, initial):
    leaves = []

    def visit(z, probability, gamma, pivot=None, phase_start=None):
        active = [i for i in range(z.rows) if abs(z[i]) < 1]
        if pivot is not None and pivot not in active:
            difference = phase_start - projection(v, active)
            residual = difference * v[:, pivot]
            gamma += (residual.T * residual)[0] * difference
            pivot, phase_start = None, None
        if not active:
            leaves.append((probability, z, gamma))
            return
        if pivot is None:
            pivot = active[-1]
            phase_start = projection(v, active)
        others = [i for i in active if i != pivot]
        direction = sp.zeros(z.rows, 1)
        direction[pivot] = 1
        if others:
            columns = v[:, others]
            coefficients = -(columns.T * columns).inv() * columns.T * v[:, pivot]
            for i, coefficient in zip(others, coefficients):
                direction[i] = coefficient
        lower, upper = -sp.oo, sp.oo
        for i in active:
            if direction[i]:
                endpoints = [(-1-z[i])/direction[i], (1-z[i])/direction[i]]
                lower = max(lower, min(endpoints))
                upper = min(upper, max(endpoints))
        assert lower < 0 < upper
        for delta, chance in [(lower, upper/(upper-lower)),
                              (upper, -lower/(upper-lower))]:
            visit(z+delta*direction, probability*chance, gamma, pivot, phase_start)

    visit(initial, sp.Integer(1), sp.zeros(v.rows))
    return leaves


def main():
    r = sp.Rational
    v = sp.Matrix([[-1, r(1, 3), 0], [0, r(2, 3), r(3, 5)],
                   [0, -r(2, 3), -r(4, 5)]])
    initial = sp.Matrix([r(99, 100), r(1, 2), -r(4, 5)])
    theta = sp.Matrix([8, -3, 3])
    assert all((v[:, i].T*v[:, i])[0] == 1 for i in range(3))
    leaves = exact_walk(v, initial)
    assert sum(p for p, _, _ in leaves) == 1
    assert sum((p*z for p, z, _ in leaves), sp.zeros(3, 1)) == initial
    gamma_mean = sum((p*g for p, _, g in leaves), sp.zeros(3))
    proxy = (theta.T*gamma_mean*theta)[0]/2
    counter = [(p, z, g) for p, z, g in leaves if z == -sp.ones(3, 1)]
    assert len(counter) == 1
    rare_probability, rare_sign, _ = counter[0]
    rare_linear = (theta.T*v*(rare_sign-initial))[0]
    assert rare_probability == r(13, 14000)
    assert rare_linear == r(469, 25)
    assert rare_linear-proxy > 8
    exp8_lower = sum(r(8**j, sp.factorial(j)) for j in range(12))
    assert exp8_lower > 14000/r(13)
    # Therefore p*exp(rare_linear-proxy)>p*exp(8)>1, rigorously.
    mgf = sum(p*sp.exp((theta.T*v*(z-initial))[0]) for p, z, _ in leaves)
    compensated_mgf = sum(p*sp.exp((theta.T*v*(z-initial))[0]
                                  -(theta.T*g*theta)[0]/2)
                          for p, z, g in leaves)
    print(json.dumps({
        "status": "exact rational averaged-compensator falsification passed",
        "leaf_count": len(leaves),
        "rare_probability": str(rare_probability),
        "rare_linear": str(rare_linear),
        "half_mean_compensator": str(proxy),
        "rare_linear_minus_proxy": str(rare_linear-proxy),
        "exp8_rational_lower": str(exp8_lower),
        "averaged_compensator": [[str(x) for x in row] for row in gamma_mean.tolist()],
        "log_mgf": float(sp.log(mgf)),
        "log_mgf_minus_half_mean_compensator": float(sp.log(mgf)-proxy),
        "log_random_compensated_mgf": float(sp.log(compensated_mgf)),
    }, indent=2))


if __name__ == "__main__":
    main()
