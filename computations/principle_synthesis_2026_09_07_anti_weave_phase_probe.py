"""Numerical diagnostic only: paired-row anti-weave root phase.

Uses the archived exact ternary scalar variational expression. No number
from this exploratory floating-point search is a proof certificate.
"""
import argparse
import json
import math

import numpy as np
from scipy.optimize import minimize


def ternary_energy(t, density, amplitude_sq):
    if density <= 1e-12:
        return 0.0

    def objective(point):
        precision, z = point
        if precision <= 0 or precision > t or z < 0 or z > 1:
            return 1e6
        c = 0.25 * math.log(precision * (2 * t - precision) / (t * t))
        a = math.expm1(precision * amplitude_sq * z * z)
        b = math.cosh(2 * precision * amplitude_sq * z) - 1
        if a * b <= 1e-20:
            gain = 0.0
        else:
            weight = max(0.0, min(1.0, (density * b - a) /
                                  ((1 - density) * a * b)))
            gain = density * math.log1p(weight * b) - math.log1p(weight * a)
        return -(c - precision * density * amplitude_sq + gain)

    starts = [(t * s, z) for s in (.03, .1, .25, .5, .8, .99)
              for z in (.0, .25, .5, .75, 1.)]
    results = [minimize(objective, point, method="Nelder-Mead",
                        bounds=((1e-9, t), (0, 1)),
                        options={"xatol": 1e-9, "fatol": 1e-10,
                                 "maxiter": 500}) for point in starts]
    best = min(results, key=lambda result: result.fun)
    return float(-best.fun)


def entropy(theta):
    if theta <= 0 or theta >= 1:
        return 0.0
    return -theta * math.log(theta) - (1 - theta) * math.log1p(-theta)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=float, default=.96)
    parser.add_argument("--t", type=float, default=4.85)
    parser.add_argument("--points", type=int, default=21)
    args = parser.parse_args()
    rows = []
    for theta in np.linspace(0, .5, args.points):
        eplus = ternary_energy(args.t, args.p * theta, 2 / args.p)
        eminus = ternary_energy(args.t, args.p * (1 - theta), 2 / args.p)
        row = args.p / 2 * (math.log(2) + entropy(theta)) + (eplus + eminus) / 2
        cap = (args.t + row) / (2 * args.t * math.sqrt(args.p))
        rows.append({"theta": float(theta), "Eplus": eplus,
                     "Eminus": eminus, "row_exponent": row,
                     "display_cap": cap})
    print(json.dumps({"status": "FLOATING_DIAGNOSTIC_NOT_CERTIFICATE",
                      "p": args.p, "t": args.t, "rows": rows,
                      "largest_grid_cap": max(x["display_cap"] for x in rows)},
                     indent=2))


if __name__ == "__main__":
    main()
