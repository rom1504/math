"""Diagnostic optimization of an actual two-field response mask.

Floating-point Gaussian quadrature is exploratory, NOT a certificate.
The mask is H(g,z)=1{|g+q*z| <= max(0,r(g))}; its rooted first mask
is h(g)=(E_z H(g,z)-p)/c, with p=E H and c^2=Var(E_z H).
Only the Gaussian integral is evaluated; no signing surrogates are used.
"""
import argparse
import json
import numpy as np
from scipy.special import ndtr
from scipy.optimize import minimize, differential_evolution
from numpy.polynomial.legendre import leggauss


def phi(x):
    return np.exp(-np.asarray(x)**2/2)/np.sqrt(2*np.pi)


def evaluator(order=640):
    nodes, weights = leggauss(order)
    g = 8*nodes
    w = 8*weights*phi(g)

    def run(theta, detail=False):
        q, a, b, d = theta[:4]
        if len(theta) == 4:
            radius = a + b*np.exp(-d*g*g)
        else:
            # Polynomial on a bounded coordinate gives a flexible even mask.
            coord = g*g/(1+g*g)
            radius = sum(t*coord**i for i, t in enumerate(theta[1:]))
        radius = np.maximum(radius, 0)
        lo, hi = (-radius-g)/q, (radius-g)/q
        k = ndtr(hi)-ndtr(lo)
        p = float(w@k)
        c = float(np.sqrt(max(w@(k*k)-p*p, 0)))
        if c < 1e-12 or p < 1e-12:
            return {} if detail else -1.0
        mu = p*g
        z0 = -mu/c
        full = 2*c*phi(mu/c)+mu*(2*ndtr(mu/c)-1)
        def signed(l, u):
            return mu*(ndtr(u)-ndtr(l))+c*(phi(l)-phi(u))
        split = np.minimum(np.maximum(z0, lo), hi)
        inside = signed(split, hi)-signed(lo, split)
        value = float(w@(full-inside))
        if detail:
            return {"parameters": list(map(float, theta)), "p": p, "c": c,
                    "q_actual": c/p, "value": value,
                    "Gaussian_quadrature_order": order,
                    "evidence": "floating-point diagnostic, not certified"}
        return value
    return run


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--order", type=int, default=640)
    parser.add_argument("--iterations", type=int, default=160)
    args = parser.parse_args()
    evaluate = evaluator(args.order)
    result = differential_evolution(lambda x: -evaluate(x),
        [(0.05, 2), (0.01, 3), (-2, 3), (.01, 5)],
        seed=20260905, maxiter=args.iterations, popsize=18, tol=1e-9,
        polish=True)
    best = minimize(lambda x: -evaluate(x), result.x, method="Nelder-Mead",
        options={"maxiter": 1200, "xatol": 1e-10, "fatol": 1e-12})
    report = evaluate(best.x, True)
    report["checks"] = [evaluator(k)(best.x, True) for k in (320, 960, 1600)]
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
