#!/usr/bin/env python3
"""Exploratory, complete face enumeration of a two-type weighted limit.

This finite-dimensional calculation is floating-point discovery, NOT a
certificate for full-sign matrices or asymptotic minimizers. Zero-curvature
stationary faces can be moved to a boundary face without changing the value;
thus strictly negative definite free Hessians and vertices suffice in exact
arithmetic. Rational reconstruction/proof is a separate obligation.
"""
import itertools
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def quadratic_box_max(h):
    n = len(h)
    answer = (-float("inf"), None, None)
    for mask in range(1 << n):
        free = [i for i in range(n) if mask >> i & 1]
        fixed = [i for i in range(n) if not (mask >> i & 1)]
        ff = h[np.ix_(free, free)]
        if free and np.linalg.eigvalsh(ff).max() >= -1e-10:
            continue
        signs = np.array(list(itertools.product([-1.0, 1.0], repeat=len(fixed))))
        x = np.zeros((len(signs), n))
        x[:, fixed] = signs
        if free:
            x[:, free] = np.linalg.solve(ff, -h[np.ix_(free, fixed)] @ signs.T).T
            x = x[np.abs(x[:, free]).max(axis=1) <= 1 + 1e-10]
        if len(x):
            vals = np.einsum("bi,ij,bj->b", x, h, x) / 2
            best = int(vals.argmax())
            if vals[best] > answer[0]:
                answer = (float(vals[best]), x[best].tolist(), free)
    return answer


def main():
    records = []
    for t in [0, .1, .25, 2 ** .5 - 1, .5, .75, 1]:
        macro = np.array([[1, t], [t, -1]])
        block = np.array([0, 0, 1, 1])
        sign = np.array([1, -1, 1, -1])
        a = macro[np.ix_(block, block)] / 16
        b = sign[:, None] * a * sign[None, :]
        h = np.block([[a, b], [b, -a]])
        mx, witness, free = quadratic_box_max(h)
        seed = (1 + t*t) / 8
        records.append({"cross_mean": t, "seed_cap": seed,
            "balanced_twist_cap_candidate": mx, "ratio": mx / seed,
            "witness": witness, "free_coordinates": free,
            "untwisted_cap_candidate": quadratic_box_max(np.block([[a,a],[a,-a]]))[0]})
    result = {"evidence": "floating-point exhaustive stationary-face discovery, not an exact certificate",
        "model": "two equal constant-sign blocks with means +1,-1 and cross mean t; balanced within-block switching",
        "records": records}
    out = ROOT / "computations/results/twisted_chiral_two_type_twist_2026_09_19.json"
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(records, indent=2))


if __name__ == "__main__":
    main()
