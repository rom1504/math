"""Finite LP diagnostic for fractional ternary-support energies."""
import itertools
import json
from pathlib import Path
from fractions import Fraction
import numpy as np
from scipy.optimize import linprog


def compute(a):
    n = len(a)
    vectors = np.array(list(itertools.product([-1, 0, 1], repeat=n)), dtype=np.int64)
    energy = np.einsum("bi,ij,bj->b", vectors, a, vectors) // 2
    support = ((vectors != 0) * (1 << np.arange(n))).sum(axis=1)
    incidence = ((np.arange(1, 1 << n)[:, None] >> np.arange(n)) & 1).T
    out = []
    for sigma in [1, -1]:
        values = np.full(1 << n, -n*n, dtype=np.int64)
        np.maximum.at(values, support, sigma * energy)
        lp = linprog(-values[1:].astype(float), A_ub=incidence,
                     b_ub=np.ones(n), bounds=(0, None), method="highs")
        assert lp.success
        primal = [Fraction(float(w)).limit_denominator(1000000) for w in lp.x]
        dual = [Fraction(float(-w)).limit_denominator(1000000) for w in lp.ineqlin.marginals]
        assert all(w >= 0 for w in primal + dual)
        assert all(sum(primal[s] for s in range((1 << n)-1) if incidence[i, s]) <= 1 for i in range(n))
        assert all(sum(dual[i] for i in range(n) if incidence[i, s]) >= int(values[s+1]) for s in range((1 << n)-1))
        exact = sum(primal[s] * int(values[s+1]) for s in range((1 << n)-1))
        assert exact == sum(dual)
        out.append(dict(polarity=sigma, value=float(-lp.fun),
                        exact_value=str(exact), exact_primal_dual_checked=True,
                        weights=[dict(support=int(s+1), weight=float(w),
                                      support_energy=int(values[s+1]))
                                 for s, w in enumerate(lp.x) if w > 1e-8],
                        dual=(-lp.ineqlin.marginals).tolist()))
    return out


def main():
    output = []
    for n in range(3, 11):
        source = Path(f"computations/results/exact_m{n}.json")
        if source.exists():
            a = np.array(json.loads(source.read_text())["matrix"], dtype=np.int64)
            output.append(dict(n=n, source=str(source), packing=compute(a)))
    Path("computations/results/flatify_construct_2026_09_07_partial_packing.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps([dict(n=row["n"], values=[p["exact_value"] for p in row["packing"]]) for row in output], indent=2))


if __name__ == "__main__":
    main()
