"""Floating lower witnesses only for a conditional-variance envelope."""
import argparse
import json
import numpy as np


def entropy(x):
    with np.errstate(divide="ignore", invalid="ignore"):
        return -np.where(x > 0, x * np.log(x), 0) - np.where(
            x < 1, (1 - x) * np.log1p(-x), 0)


def run(mesh):
    p, t = 31 / 32, 4
    z = np.linspace(0, 1, mesh)
    s = np.linspace(0, 1, mesh)
    hull = []
    for a in z:
        variance = (a - a * a * s * s) / p
        rho = 4 * t * variance / (np.sqrt(1 + 16 * t * t * variance**2) + 1)
        gaussian = -t * variance * (1 - rho) + .25 * np.log1p(-rho * rho)
        value = float(np.max(entropy(a) + a * entropy((1 + s) / 2) + gaussian))
        while len(hull) > 1 and (
            (hull[-1][1] - hull[-2][1]) * (a - hull[-1][0])
            <= (value - hull[-1][1]) * (hull[-1][0] - hull[-2][0])
        ):
            hull.pop()
        hull.append((float(a), value))
    for left, right in zip(hull, hull[1:]):
        if left[0] <= p <= right[0]:
            value = (left[1] * (right[0] - p) + right[1] * (p - left[0])) / (
                right[0] - left[0])
            return {"mesh": mesh, "p": p, "t": t,
                    "offset": float(value - entropy(p) + t * (1 - np.sqrt(p))),
                    "posterior_z_support": [left, right],
                    "status": "floating grid LOWER witness; not an upper certificate"}
    raise AssertionError("The mesh hull must cover the source mass")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh", type=int, default=2001)
    args = parser.parse_args()
    if args.mesh < 2:
        parser.error("mesh must be at least two")
    print(json.dumps(run(args.mesh), indent=2))
