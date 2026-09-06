"""Finite non-Gaussian sequential-label tests; floating diagnostics only."""
import json
import math
from pathlib import Path

import numpy as np


def entropy(p):
    p = np.asarray(p)
    return float(-np.sum(p[p > 0] * np.log(p[p > 0])))


def aggregate(keys, probs):
    out = {}
    for key, p in zip(keys, probs):
        out[key] = out.get(key, 0.0) + p
    return out


def information(xkeys, lkeys, probs):
    xp = aggregate(xkeys, probs)
    lp = aggregate(lkeys, probs)
    joint = aggregate(list(zip(xkeys, lkeys)), probs)
    return entropy(list(xp.values())) + entropy(list(lp.values())) - entropy(list(joint.values()))


def reward(values, lkeys, probs, g):
    groups = {}
    for x, key, p in zip(values, lkeys, probs):
        if key not in groups:
            groups[key] = [0.0, 0.0, 0.0]
        groups[key][0] += p
        groups[key][1] += p * x
        groups[key][2] += p * x * x
    answer = 0.0
    for mass, first, second in groups.values():
        if mass > 0:
            answer += mass * g(max(0.0, second / mass - (first / mass) ** 2))
    return answer


def gaussian_g(t):
    def g(v):
        z = 4.0 * t * v
        rho = z / (math.sqrt(1.0 + z * z) + 1.0)
        return -t * v * (1.0 - rho) + 0.25 * math.log1p(-rho * rho)
    return g


def test_case(rng, d, case):
    count = 9 if d == 2 else 12
    atoms = rng.integers(-3, 4, (count, d)).astype(float)
    if case % 7 == 0:
        atoms[:, 1] = atoms[:, 0] ** 2  # nonlinear deterministic conditioning
    if case % 11 == 0:
        atoms[:, 1] = atoms[:, 0]  # singular covariance
    if d == 2 and case % 2 == 0:
        orthogonal = np.array([[1, 1], [1, -1]]) / math.sqrt(2)
    else:
        orthogonal = np.linalg.qr(rng.normal(size=(d, d)))[0]
    children = atoms @ orthogonal.T
    probs = rng.dirichlet(np.full(count, 0.2 if case % 3 == 0 else 1.0))
    labels_per_child = 2
    channels = []
    for i in range(d):
        unique, inverse = np.unique(np.round(children[:, i], 11), return_inverse=True)
        kernel = rng.dirichlet(np.full(labels_per_child, 0.15 if case % 5 == 0 else 1.0), size=len(unique))
        channels.append(kernel[inverse])
    combinations = list(np.ndindex(*(labels_per_child,) * d))
    atom_indices, label_rows, joint = [], [], []
    for j in range(count):
        for labels in combinations:
            atom_indices.append(j)
            label_rows.append(labels)
            joint.append(probs[j] * math.prod(channels[i][j, labels[i]] for i in range(d)))
    joint = np.array(joint)
    expanded_a = atoms[atom_indices]
    expanded_u = children[atom_indices]
    g = gaussian_g(float(10 ** rng.uniform(-2, 2)))
    parent_reward = parent_info = child_reward = child_info = 0.0
    for i in range(d):
        parent_labels = [(label_rows[j], tuple(expanded_a[j, :i])) for j in range(len(joint))]
        child_labels = [labels[i] for labels in label_rows]
        parent_reward += reward(expanded_a[:, i], parent_labels, joint, g)
        parent_info += information(list(expanded_a[:, i]), parent_labels, joint)
        child_reward += reward(expanded_u[:, i], child_labels, joint, g)
        child_info += information(list(np.round(expanded_u[:, i], 11)), child_labels, joint)
    total_correlation = sum(entropy(list(aggregate(list(atoms[:, i]), probs).values())) for i in range(d))
    total_correlation -= entropy(list(aggregate([tuple(row) for row in atoms], probs).values()))
    return {
        "dimension": d,
        "case": case,
        "reward_slack": parent_reward - child_reward,
        "information_slack": child_info + total_correlation - parent_info,
        "total_slack": parent_reward - parent_info - child_reward + child_info + total_correlation,
    }


def main():
    rng = np.random.default_rng(202609067)
    results = [test_case(rng, d, case) for d in (2, 3, 4) for case in range(160)]
    summary = {
        "status": "floating diagnostics, not a proof",
        "seed": 202609067,
        "cases": len(results),
        "minimum_reward_slack": min(results, key=lambda r: r["reward_slack"]),
        "minimum_information_slack": min(results, key=lambda r: r["information_slack"]),
        "minimum_total_slack": min(results, key=lambda r: r["total_slack"]),
        "results": results,
    }
    assert min(r["reward_slack"] for r in results) > -1e-8
    assert min(r["information_slack"] for r in results) > -1e-8
    assert min(r["total_slack"] for r in results) > -1e-8
    output = Path("computations/decisive_bridge_supersolution_audit_2026_09_06.json")
    output.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "results"}, indent=2))


if __name__ == "__main__":
    main()
