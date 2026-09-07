"""Exact all-spin/edge-law regression for the physical stability correction.

At t=(k/2)log(2), tilted weights are integer powers of two after a common
shift. Finner is checked after squaring, using Python integers only. This
checks finite normalization/diagonal/factorization, not an asymptotic bound.
"""
import itertools
import json
from pathlib import Path

import numpy as np


def main():
    m, k = 4, 3
    h2 = np.array([[1, 1], [1, -1]], dtype=np.int64)
    h = np.kron(h2, h2)
    selectors = [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
    # Deterministic different row/column choices; no ignored data dependency.
    bases = [h[selectors[i]][:, np.roll(np.arange(m), i)] for i in range(m)]
    edges = list(itertools.combinations(range(m), 2))
    models = []
    for signs in itertools.product([-1, 1], repeat=len(edges)):
        s = np.eye(m, dtype=np.int64)
        for (i, j), z in zip(edges, signs):
            s[i, j] = s[j, i] = z
        c = np.empty((m*k, m*k), dtype=np.int64)
        for i in range(m):
            for j in range(m):
                c[i*k:(i+1)*k, j*k:(j+1)*k] = s[i, j]*np.outer(bases[i][:, j], bases[j][:, i])
        np.fill_diagonal(c, 0)
        assert np.array_equal(c, c.T)
        models.append(c)
    models = np.array(models)
    caps = np.zeros(len(models), dtype=np.int64)
    stable_caps = np.zeros_like(caps)
    checks = strict = 0
    smallest_num = smallest_den = None
    witness = None
    for bits in itertools.product([-1, 1], repeat=m*k-1):
        x = np.array((1,)+bits, dtype=np.int64)
        fields = np.einsum("sij,j->si", models, x)
        energies2 = fields @ x
        assert not np.any(energies2 % 2)
        energies = energies2//2
        caps = np.maximum(caps, np.abs(energies))
        for sigma in [-1, 1]:
            z = sigma*energies
            shifted = z-int(z.min())
            weights = [1 << int(power) for power in shifted]
            total = sum(weights)
            fibre_stable = (sigma*fields*x).reshape(len(models), m, k).min(axis=2)>=0
            all_stable = fibre_stable.all(axis=1)
            stable_caps = np.maximum(stable_caps, np.where(all_stable, z, 0))
            joint = sum(w for w, valid in zip(weights, all_stable) if valid)
            marginals = [sum(w for w, valid in zip(weights, fibre_stable[:, i]) if valid)
                         for i in range(m)]
            lhs = joint**2*total**(m-2)
            rhs = 1
            for q in marginals:
                rhs *= q
            assert lhs <= rhs
            checks += 1
            strict += lhs < rhs
            if joint and (smallest_num is None or joint*smallest_den<smallest_num*total):
                smallest_num, smallest_den = joint, total
                witness = {"spin": x.tolist(), "sigma": sigma,
                           "joint_weight": str(joint), "total_weight": str(total),
                           "fibre_weights": list(map(str, marginals))}
    assert np.array_equal(caps, stable_caps)
    report = {"status": "EXACT_INTEGER_PHYSICAL_STABILITY_FINNER_PASS",
              "m": m, "k": k, "actual_sign_models": len(models),
              "projective_spins": 2**(m*k-1), "oriented_checks": checks,
              "strict_finner_checks": strict, "caps": caps.tolist(),
              "max_equals_stable_max_in_every_model": True,
              "smallest_nonzero_stability_probability_witness": witness,
              "scope": "Finite exact regression, not an asymptotic cap certificate."}
    output = Path(__file__).resolve().parent/"results/principle_director_2026_09_07_stable_weave_exact.json"
    output.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({k: v for k, v in report.items() if k not in {"caps", "smallest_nonzero_stability_probability_witness"}}, indent=2))


if __name__ == "__main__":
    main()
