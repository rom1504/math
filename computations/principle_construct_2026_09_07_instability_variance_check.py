"""Bounded exact-cube checks of the variance proof, not asymptotic evidence."""
import json
from pathlib import Path
import numpy as np

root = Path(__file__).resolve().parents[1]
records = []
rng = np.random.default_rng(20260907)
for n in range(3, 15):
    x = 1 - 2 * ((np.arange(1 << n)[:, None] >> np.arange(n)) & 1)
    matrices = [("all_positive", np.ones((n, n), dtype=np.int64)-np.eye(n, dtype=np.int64))]
    random_a = np.triu(rng.choice([-1, 1], size=(n, n)), 1)
    matrices.append(("random", random_a + random_a.T))
    exact_path = root / f"computations/results/exact_m{n}.json"
    if exact_path.exists():
        raw = json.loads(exact_path.read_text())
        if "matrix" in raw:
            matrices.append(("stored_exact", np.asarray(raw["matrix"], dtype=np.int64)))
    for name, a in matrices:
        fields = x @ a
        f = np.abs(fields).sum(axis=1)
        quadratic = (x * fields).sum(axis=1)
        instability = np.maximum(-x * fields, 0).sum(axis=1)
        assert np.array_equal(2 * instability, f - quadratic)
        signs = np.sign(fields)
        gradient = signs @ a
        rhs = 2 * np.square(gradient).sum(axis=1).mean()
        var_f = f.var()
        assert var_f <= rhs + 1e-8
        assert abs(quadratic.var() - 2*n*(n-1)) < 1e-8
        exact_mean = n/2 * np.abs(x[:, :n-1].sum(axis=1)).mean()
        assert instability.mean() == exact_mean
        a2 = a @ a
        tr4 = np.square(a2).sum()
        records.append(dict(n=n, family=name, mean=float(exact_mean),
                            variance=float(instability.var()),
                            variance_over_bound_scale=float(instability.var()/(n*n+tr4/n)),
                            convex_poincare_lhs=float(var_f),
                            convex_poincare_rhs=float(rhs), tr4=int(tr4)))
result = dict(status="PASS", records=records,
              scope="Finite exact-cube identities and convex variance bound only; universal constants are proved analytically.")
target = root / "computations/results/principle_construct_2026_09_07_instability_variance_check.json"
target.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(dict(status="PASS", cases=len(records),
                     max_observed_ratio=max(r["variance_over_bound_scale"] for r in records))))
