#!/usr/bin/env python3
"""Numerical stress test of the explicit flat-diagonal quadratic witness.

All mathematical inequalities have analytic proofs. This script checks
the adaptive Bloch signs, tensor ordering, uniform diagonal, and finite
cap/L2 commutator estimate. It is not an SDP optimality certificate.
"""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np


def helpers():
    source = Path(__file__).with_name("bh_boundedness_director_2026_09_18_quadratic_phase.py")
    spec = importlib.util.spec_from_file_location("phase_helpers", source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def check(theta, linear, helper):
    n = len(theta)
    size = 1 << n
    verts = np.arange(size)
    spins = 1 - 2*((verts[:, None] >> np.arange(n)) & 1)
    g = np.exp(1j*(np.einsum("bi,ij,bj->b", spins, theta, spins)/2 + spins@linear))
    J = np.zeros((size, size), complex)
    for i in range(n):
        J[verts, verts ^ (1 << i)] = 0.5
    A = g.conj()[:, None]*J*g[None, :] - J
    theta0 = (theta+np.pi/4) % (np.pi/2) - np.pi/4
    np.fill_diagonal(theta0, 0)
    alpha0 = (linear+np.pi/4) % (np.pi/2) - np.pi/4
    gauge = np.rint((theta-theta0)/(np.pi/2)).astype(int).sum(axis=1)
    gauge += np.rint((linear-alpha0)/(np.pi/2)).astype(int)
    sigma = 1 - 2*(gauge % 2)
    rho = np.zeros((size, size), complex)
    row_witness_sum = 0.0
    for z in spins:
        M = np.imag(np.exp(2j*alpha0) * np.prod(
            np.cos(2*theta0)+0.5j*z[None, :]*np.sin(2*theta0), axis=1))
        Y = sigma*np.sign(M)*np.sqrt(3)/2
        # Either zero Y choice is allowed when M=0; use the pure positive one.
        Y[M == 0] = np.sqrt(3)/2
        state = np.array([[1.0+0j]])
        for i in reversed(range(n)):
            local = np.array([[1+z[i]/2, -1j*Y[i]],
                              [1j*Y[i], 1-z[i]/2]])/2
            state = np.kron(state, local)
        rho += state/size
        row_witness_sum += np.sqrt(3)/4*np.abs(M).sum()/size
    assert np.max(np.abs(np.diag(rho)-1/size)) < 1e-12
    assert np.min(np.linalg.eigvalsh(rho)) > -1e-11
    witness = float(np.trace(rho@A).real)
    assert abs(abs(witness)-row_witness_sum) < 1e-10
    coef = helper.fwht(g)
    prob = np.abs(coef)**2
    degree = np.array([bin(int(v)).count("1") for v in verts])
    influence = float(prob@degree)
    H = helper.entropy(prob)
    certificate = max(influence, abs(witness))
    assert H <= 64*certificate+1e-10
    finite = []
    for D in range(n+1):
        truncated = coef.copy()
        truncated[degree>D] = 0
        Q = helper.fwht(truncated)*size
        cap = float(np.max(np.abs(Q)))
        delta = float(np.linalg.norm(Q-g)/np.sqrt(size))
        assert D*cap+n*delta+1e-10 >= certificate
        finite.append({"degree_budget": D, "cap": cap, "L2_error": delta,
                       "finite_rhs": D*cap+n*delta})
    return {"n": n, "theta": theta.tolist(), "linear": linear.tolist(),
            "entropy": H, "influence": influence,
            "flat_product_witness": witness,
            "finite_projection_checks": finite}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    helper = helpers()
    rng = np.random.default_rng(18462026)
    rows = []
    for n in range(2, 7):
        for scale in (1e-4, 0.05, 0.5, 3):
            theta = np.triu(rng.normal(size=(n, n))*scale, 1)
            theta += theta.T
            rows.append(check(theta, rng.normal(size=n)*scale, helper))
    result = {"status": "NUMERICAL_DIAGNOSTIC_ONLY", "seed": 18462026,
              "cases": len(rows), "rows": rows}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": "PASS_DIAGNOSTICS", "cases": len(rows),
                      "finite_checks": sum(len(r["finite_projection_checks"]) for r in rows)}))


if __name__ == "__main__":
    main()
