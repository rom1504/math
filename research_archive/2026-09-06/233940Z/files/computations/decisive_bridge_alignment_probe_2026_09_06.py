"""Search explicit common-precision child channels at a fixed ternary pair."""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
from decisive_bridge_supersolution_gap_2026_09_06 import envelope, ent


def common_precision(p, x, t, lam, details=False):
    y = np.linspace(0, max(x), 501)
    c = 0.25 * np.log(lam * (2 * t - lam) / (t * t))
    kernels = np.array([np.exp(c) * (np.exp(-lam * (a - y) ** 2) + np.exp(-lam * (a + y) ** 2)) / 2 for a in x]).T
    value = envelope(p, kernels, details=details)
    if details:
        value["reproductions"] = [float(y[i]) for i in value["active"]]
        value["lambda"] = float(lam)
    return value


def optimize(p, x, t):
    grid = np.geomspace(0.02, t, 90)
    values = [common_precision(p, x, t, lam) for lam in grid]
    candidates = [(v, l) for v, l in zip(values, grid)]
    for i in range(1, len(grid) - 1):
        if values[i] >= max(values[i-1], values[i+1]):
            result = minimize_scalar(lambda lam: -common_precision(p, x, t, lam),
                                     bounds=(grid[i-1], grid[i+1]), method="bounded",
                                     options={"xatol": 1e-10})
            candidates.append((-result.fun, result.x))
    _, lam = max(candidates)
    return common_precision(p, x, t, lam, details=True)


def main():
    p, t, cross = 31 / 32, 4., 1 / 24
    a = 1 / np.sqrt(p)
    both, zero = p - cross / 2, 1 - p - cross / 2
    x = [0, a / np.sqrt(2), a * np.sqrt(2)]
    plus = optimize([zero, cross, both], x, t)
    minus = optimize([zero + both, cross, 0], x, t)
    root = optimize([1-p, p], [0, a], t)
    pair_entropy = ent([zero] + [cross / 4] * 4 + [both / 2] * 2)
    information = 2 * ent([1-p,p/2,p/2]) - pair_entropy
    lower = (plus["lower_feasible"] + minus["lower_feasible"] - information) / 2
    result = {"status": "floating feasible child channels; root global upper not certified",
              "p":p,"t":t,"cross":cross,"plus":plus,"minus":minus,"root":root,
              "pair_information":information,"policy_lower":lower,
              "gap_against_root_diagnostic":lower-root["upper_grid"]}
    Path("computations/decisive_bridge_alignment_probe_2026_09_06.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__ == "__main__":
    main()
