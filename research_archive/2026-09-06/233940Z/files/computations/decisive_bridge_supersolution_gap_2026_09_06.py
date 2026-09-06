"""Numerical finite-grid Gaussian-kernel dual for T and one Bellman step.

No continuum upper certificate: every reported value is a floating diagnostic.
"""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize


def kernel_grid(x, t, points=121):
    xmax = max(x)
    z = 4 * t * xmax * xmax
    rho = z / (1 + np.sqrt(1 + z * z))
    lower = t * (1 - rho)
    lam, y = np.meshgrid(np.geomspace(lower, t, points), np.linspace(0, xmax, points))
    lam, y = lam.ravel(), y.ravel()
    factor = (lam * (2 * t - lam) / (t * t)) ** 0.25
    return np.array([factor * (np.exp(-lam * (a - y) ** 2) + np.exp(-lam * (a + y) ** 2)) / 2 for a in x]).T


def envelope(p, kernels, details=False):
    p = np.asarray(p)
    keep = p > 1e-12
    p = p[keep]
    k = kernels[:, keep]
    b = p.copy()
    active = list(set(np.argmax(k, axis=0).tolist()))
    for iteration in range(100):
        constraints = k[active]
        result = minimize(lambda b: -np.dot(p, np.log(b / p)), b,
                          jac=lambda b: -p / b,
                          bounds=[(1e-12, 1e4)] * len(p),
                          constraints={"type": "ineq", "fun": lambda b: 1 - constraints @ b,
                                       "jac": lambda b: -constraints},
                          method="SLSQP", options={"ftol": 1e-12, "maxiter": 200})
        b = result.x
        values = np.sum(k * b, axis=1)
        peak = int(np.argmax(values))
        maximum = values[peak]
        if maximum <= 1 + 1e-9 or peak in active:
            # Exact feasible rescaling for the finite grid, modulo floats.
            upper = float(-np.dot(p, np.log(b / p)) + np.log(maximum))
            if not details:
                return upper
            selected = k[active].T
            weights = np.full(len(active), 1 / len(active))
            primal = minimize(lambda q: -np.dot(p, np.log(selected @ q)), weights,
                              jac=lambda q: -selected.T @ (p / (selected @ q)),
                              bounds=[(0, 1)] * len(active),
                              constraints={"type": "eq", "fun": lambda q: q.sum() - 1,
                                           "jac": lambda q: np.ones_like(q)},
                              method="SLSQP", options={"ftol": 1e-13, "maxiter": 500})
            return {"upper_grid": upper, "lower_feasible": float(-primal.fun),
                    "active": active, "weights": primal.x.tolist()}
        active.append(peak)
    raise RuntimeError("Cutting-plane dual did not settle")


def ent(p):
    p = np.asarray(p)
    return float(-np.sum(p[p > 0] * np.log(p[p > 0])))


def main():
    p, t = 31 / 32, 4.0
    a = 1 / np.sqrt(p)
    root = envelope([1 - p, p], kernel_grid([0, a], t, 201))
    kernels = kernel_grid([0, a / np.sqrt(2), a * np.sqrt(2)], t, 151)
    root_entropy = ent([1 - p, p / 2, p / 2])
    records = []
    for cross in np.linspace(0, 2 * (1 - p), 25):
        both = p - cross / 2
        zero = 1 - p - cross / 2
        for fraction in np.linspace(0.5, 1, 31):
            same, opposite = both * fraction, both * (1 - fraction)
            plus = envelope([zero + opposite, cross, same], kernels)
            minus = envelope([zero + same, cross, opposite], kernels)
            pair_entropy = ent([zero] + [cross / 4] * 4 + [same / 2] * 2 + [opposite / 2] * 2)
            information = 2 * root_entropy - pair_entropy
            value = (plus + minus - information) / 2
            records.append({"cross": float(cross), "same_fraction": float(fraction),
                            "value": value, "plus": plus, "minus": minus,
                            "information": information})
    best = max(records, key=lambda r: r["value"])
    output = {"status": "finite grids and floating optimization ONLY; no continuum Bellman upper bound",
              "root_T_grid": root, "best_policy_grid": best,
              "apparent_gap": root - best["value"], "records": records}
    Path("computations/decisive_bridge_supersolution_gap_2026_09_06.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({k: v for k, v in output.items() if k != "records"}, indent=2))


if __name__ == "__main__":
    main()
