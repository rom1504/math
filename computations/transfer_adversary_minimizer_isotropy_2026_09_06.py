"""Finite LP screen with exact rational replay; no asymptotic inference.

The stored minimizer labels are imported from their stated source certificates.
Caps, proposed isotropic laws, and separating functionals are checked afresh.
"""
import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.optimize import linprog
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]


def source(path):
    payload = (ROOT / path).read_bytes()
    return json.loads(payload), hashlib.sha256(payload).hexdigest()


def cases():
    for n in range(3, 9):
        path = f"computations/results/m{n}_minimizer_orbits.json"
        data, digest = source(path)
        for item in data["classes"]:
            yield f"n{n}_class{item['class']}", item["representative_matrix"], {
                "path": path, "file_sha256": digest,
                "json_pointer": f"/classes/{item['class']}/representative_matrix",
                "minimum_cap_imported": data["target_cap"],
                "minimum_certificate_status": data["classification"],
            }
    extras = [
        (9, "computations/results/exact_m9.json", "matrix", 12),
        (10, "computations/results/exact_m10.json", "matrix", 13),
        (11, "computations/results/nested_10_in_11_cap17.json", "matrix", 17),
        (12, "computations/results/extension_nested_m11_to_12.json", "parent_matrix", 18),
        (13, "computations/results/bridge_6_7_sign1_cap20.json", "parent_matrix", 20),
        (14, "computations/results/conference_completion_m13.json", "conference_matrix", 21),
    ]
    for n, path, key, cap in extras:
        data, digest = source(path)
        yield f"n{n}_stored", data[key], {
            "path": path, "file_sha256": digest, "json_pointer": "/" + key,
            "minimum_cap_imported": cap,
            "minimum_certificate_status": (
                "stored exact solver certificate; global lower bound not rerun here"
                if n <= 10 else "certified_m11_m12.json / certified_m13_m14.json; solver-dependent lower bounds not rerun here"
            ),
        }


def exact_weights(result, features):
    indices = np.flatnonzero(result.x > 1e-9)
    for bound in (1000, 1000000, 1000000000):
        weights = [F(float(result.x[i])).limit_denominator(bound) for i in indices]
        if sum(weights) != 1 or any(w < 0 for w in weights):
            continue
        if all(sum(w * int(features[i, e]) for i, w in zip(indices, weights)) == 0
               for e in range(features.shape[1])):
            return [(int(i), w) for i, w in zip(indices, weights)]
    rows = [[1] * len(indices)] + [[int(features[i, e]) for i in indices]
                                   for e in range(features.shape[1])]
    values = exact_linear_solution(rows, [1] + [0] * features.shape[1], result.x[indices])
    if values is None or any(w < 0 for w in values):
        return None
    return [(int(i), w) for i, w in zip(indices, values)]


def exact_linear_solution(rows, rhs, approximate):
    """Recover an exact point on the numerical solution's affine face."""
    columns = len(approximate)
    augmented = sp.Matrix([row + [int(b)] for row, b in zip(rows, rhs)])
    reduced, pivots = augmented.rref()
    if columns in pivots:
        return None
    free = [j for j in range(columns) if j not in pivots]
    values = [F(0)] * columns
    for j in free:
        values[j] = F(float(approximate[j])).limit_denominator(1000000)
    for i, j in enumerate(pivots):
        values[j] = F(reduced[i, columns]) - sum(F(reduced[i, k]) * values[k] for k in free)
    assert all(sum(F(a) * x for a, x in zip(row, values)) == b for row, b in zip(rows, rhs))
    return values


def separator(features):
    result = linprog(np.zeros(features.shape[1]), A_ub=-features.astype(float),
                     b_ub=-np.ones(len(features)), bounds=(None, None), method="highs")
    assert result.success, result.message
    weights = [F(float(y)).limit_denominator(1000000) for y in result.x]
    values = [sum(int(x) * w for x, w in zip(row, weights)) for row in features]
    assert min(values) > 0
    return {
        "coefficients_in_lexicographic_edge_order": [str(w) for w in weights],
        "minimum_on_selected_states": str(min(values)),
        "status": "exact rational separator, uniform-cube expectation zero",
    }


def check_case(name, matrix, provenance):
    a = np.asarray(matrix, dtype=np.int64)
    n = len(a)
    assert a.shape == (n, n) and np.array_equal(a, a.T)
    assert np.all(np.diag(a) == 0)
    edges = list(itertools.combinations(range(n), 2))
    assert all(abs(int(a[i, j])) == 1 for i, j in edges)
    spins = np.asarray([(1,) + x for x in itertools.product((-1, 1), repeat=n-1)], dtype=np.int64)
    features = np.asarray([spins[:, i] * spins[:, j] for i, j in edges], dtype=np.int64).T
    energies = features @ np.asarray([a[i, j] for i, j in edges])
    cap = int(np.max(np.abs(energies)))
    assert cap == provenance["minimum_cap_imported"], (name, cap)
    slack = cap - np.abs(energies)
    result = {"case": name, "n": n, "source": provenance, "cap_verified_exact": cap,
              "positive_ground_count_projective": int(np.sum(energies == cap)),
              "negative_ground_count_projective": int(np.sum(energies == -cap)),
              "ground_matrix_rank_numerical": int(np.linalg.matrix_rank(spins[slack == 0])),
              "slack_screen": []}
    for threshold in sorted(set(int(v) for v in slack)):
        keep = np.flatnonzero(slack <= threshold)
        f = features[keep]
        eq = np.vstack((np.ones(len(keep)), f.T))
        rhs = np.zeros(1 + len(edges)); rhs[0] = 1
        lp = linprog(np.zeros(len(keep)), A_eq=eq, b_eq=rhs, bounds=(0, None), method="highs")
        if lp.success:
            law = exact_weights(lp, f)
            assert law is not None, (name, threshold, "rational reconstruction failed")
            cert = [{"projective_index": int(keep[i]), "spin": spins[keep[i]].tolist(),
                     "weight": str(w), "slack": int(slack[keep[i]])} for i, w in law]
            result["slack_screen"].append({"maximum_slack": threshold,
                                           "status": "FEASIBLE, exactly replayed rational law",
                                           "law": cert})
            result["minimum_uniform_slack_for_isotropy"] = threshold
            break
        assert lp.status == 2, lp.message
        result["slack_screen"].append({"maximum_slack": threshold,
                                       "status": "INFEASIBLE, exact rational separator",
                                       "separator": separator(f)})
    else:
        raise AssertionError("The uniform full-cube law must be isotropic")

    eq = np.vstack((np.ones(len(spins)), features.T))
    rhs = np.zeros(1 + len(edges)); rhs[0] = 1
    lp = linprog(slack.astype(float), A_eq=eq, b_eq=rhs, bounds=(0, None), method="highs")
    assert lp.success
    law = exact_weights(lp, features)
    assert law is not None
    value = sum(w * int(slack[i]) for i, w in law)
    dual = [F(float(z)).limit_denominator(1000000) for z in lp.eqlin.marginals]
    valid = lambda y: y is not None and y[0] == value and all(
        y[0] + sum(z * int(f) for z, f in zip(y[1:], row)) <= int(g)
        for row, g in zip(features, slack))
    if not valid(dual):
        tight = np.flatnonzero(np.abs(eq.T @ lp.eqlin.marginals - slack) < 1e-7)
        rows = [[1] + [int(z) for z in features[i]] for i in tight]
        dual = exact_linear_solution(rows, [int(slack[i]) for i in tight], lp.eqlin.marginals)
    assert valid(dual), (name, "dual rational recovery failed")
    assert all(dual[0] + sum(z * int(f) for z, f in zip(dual[1:], row)) <= int(g)
               for row, g in zip(features, slack))
    assert dual[0] == value
    result["minimum_expected_slack_isotropic"] = {
        "exact_value": str(value),
        "dual_constant": str(dual[0]), "dual_edge_coefficients": [str(z) for z in dual[1:]],
        "law": [{"projective_index": i, "spin": spins[i].tolist(), "weight": str(w),
                 "slack": int(slack[i])} for i, w in law],
        "status": "exact primal law and pointwise dual inequality agree",
    }
    print(json.dumps({"case": name, "cap": cap, "grounds": int(np.sum(slack == 0)),
                      "minimum_support_slack": result["minimum_uniform_slack_for_isotropy"],
                      "minimum_expected_slack": str(value)}), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    results = [check_case(*case) for case in cases()]
    a4 = np.ones((4, 4), dtype=np.int64) - np.eye(4, dtype=np.int64)
    a4[2, 3] = a4[3, 2] = -1
    x4 = np.asarray([(1,) + z for z in itertools.product((-1, 1), repeat=3)], dtype=np.int64)
    h4 = np.einsum("bi,ij,bj->b", x4, a4, x4) // 2
    assert np.array_equal(np.abs(h4), 2 + x4[:, 0] * x4[:, 1] + x4[:, 2] * x4[:, 3])
    iso4 = np.asarray([[1, 1, 1, -1], [1, 1, -1, 1],
                       [1, -1, 1, 1], [1, -1, -1, -1]], dtype=np.int64)
    assert np.array_equal(iso4.T @ iso4, 4 * np.eye(4, dtype=np.int64))
    assert np.all(np.abs(np.einsum("bi,ij,bj->b", iso4, a4, iso4) // 2) == 2)
    minimizing_count = 0
    edges4 = list(itertools.combinations(range(4), 2))
    f4 = np.asarray([x4[:, i] * x4[:, j] for i, j in edges4]).T
    for signs in itertools.product((-1, 1), repeat=6):
        energy = f4 @ np.asarray(signs)
        cap = int(np.max(np.abs(energy)))
        assert cap >= 4
        if cap == 4:
            minimizing_count += 1
            assert int(np.sum(np.abs(energy) == 4)) == 2
    assert minimizing_count == 48
    report = {"status": "finite exact primal/dual certificates; no asymptotic inference",
              "spin_convention": "first coordinate +1; lexicographic tail (-1,+1)",
              "edge_convention": "lexicographic unordered pairs, zero based",
              "elementary_n4_counterexample": {
                  "matrix": a4.tolist(), "all_signings_checked": 64,
                  "global_minimizers": minimizing_count,
                  "ground_projective_count_each_minimizer": 2,
                  "absolute_energy_identity": "|H(x)|=2+x0*x1+x2*x3",
                  "every_isotropic_law_mean_slack": 2,
                  "attaining_uniform_four_atom_law": iso4.tolist(),
              },
              "cases": results}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"completed_cases": len(results), "output": str(args.output)}))


if __name__ == "__main__":
    main()
