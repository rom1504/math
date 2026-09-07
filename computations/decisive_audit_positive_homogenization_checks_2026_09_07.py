"""Replay relative-entropy homogenization and exact row-repair coupling."""
from collections import Counter
from itertools import permutations, product
import math
import numpy as np


def entropy(p):
    p = np.asarray(p)
    return -float(np.sum(p[p > 0] * np.log(p[p > 0])))


def relative(p, q):
    use = p > 0
    assert np.all(q[use] > 0)
    return float(np.sum(p[use] * np.log(p[use] / q[use])))


def compositions(n, d):
    if d == 1:
        yield (n,)
    else:
        for a in range(n + 1):
            for rest in compositions(n - a, d - 1):
                yield (a,) + rest


def main():
    rng = np.random.default_rng(20260907)
    information_checks = 0
    for _ in range(600):
        r, d = int(rng.integers(1, 6)), int(rng.integers(2, 5))
        pi = rng.dirichlet(np.ones(r))
        gamma = np.empty((r, r, d, d))
        for a in range(r):
            for b in range(a, r):
                g = rng.dirichlet(np.ones(d * d)).reshape(d, d)
                if a == b:
                    g = (g + g.T) / 2
                gamma[a, b], gamma[b, a] = g, g.T
        nu = np.einsum("b,abcd->ac", pi, gamma)
        alpha = pi[:, None] * nu
        joint = np.einsum("a,b,abcd->acbd", pi, pi, gamma)
        baseline = np.einsum("ac,bd->acbd", alpha, alpha)
        collapsed = joint.sum(axis=(0, 2))
        mean = alpha.sum(axis=0)
        full_kl = relative(joint, baseline)
        reduced_kl = relative(collapsed, mean[:, None] * mean[None, :])
        assert full_kl + 3e-13 >= reduced_kl
        conditional_entropy = sum(pi[a] * pi[b] * entropy(gamma[a, b]) for a in range(r) for b in range(r))
        expected_kl = 2 * sum(pi[a] * entropy(nu[a]) for a in range(r)) - conditional_entropy
        assert abs(full_kl - expected_kl) < 3e-13
        information_checks += 1

    repair_checks = 0
    counts = list(compositions(4, 3))
    for first in counts:
        for second in counts:
            pairs = []
            surplus, deficit = [], []
            for c in range(3):
                pairs.extend([(c, c)] * min(first[c], second[c]))
                surplus.extend([c] * max(first[c] - second[c], 0))
                deficit.extend([c] * max(second[c] - first[c], 0))
            pairs.extend(zip(surplus, deficit))
            left, right = Counter(), Counter()
            differences = sum(abs(first[c] - second[c]) for c in range(3)) // 2
            for ordering in permutations(range(4)):
                x = tuple(pairs[j][0] for j in ordering)
                y = tuple(pairs[j][1] for j in ordering)
                left[x] += 1
                right[y] += 1
                assert sum(a != b for a, b in zip(x, y)) == differences
            assert len(set(left.values())) == len(set(right.values())) == 1
            repair_checks += 1

    finite_checks = 0
    for _ in range(30):
        n, d = 4, 3
        row_slots = [[j for j in range(n) if j != i] for i in range(n)]
        rows = []
        nu = []
        for i in range(n):
            colors = tuple(map(int, rng.integers(0, d, size=n-1)))
            rows.append(list(set(permutations(colors))))
            nu.append(np.bincount(colors, minlength=d) / (n-1))
        arrays = list(product(*rows))
        probability = rng.dirichlet(np.ones(len(arrays)))
        edge_laws = {(i, j): np.zeros((d, d)) for i in range(n) for j in range(i+1, n)}
        for array, weight in zip(arrays, probability):
            for i, j in edge_laws:
                a = array[i][row_slots[i].index(j)]
                b = array[j][row_slots[j].index(i)]
                edge_laws[i, j][a, b] += weight
        nu = np.asarray(nu)
        alpha = nu / n
        joint = np.zeros((n, d, n, d))
        for (i, j), g in edge_laws.items():
            joint[i, :, j, :] = g / (n * (n-1))
            joint[j, :, i, :] = g.T / (n * (n-1))
        baseline = np.einsum("ic,jd->icjd", alpha, alpha)
        kl = relative(joint, baseline)
        average_edge_h = np.mean([entropy(g) for g in edge_laws.values()])
        identity = math.log(n / (n-1)) + 2 * np.mean([entropy(v) for v in nu]) - average_edge_h
        assert abs(kl - identity) < 4e-13
        assert entropy(probability) <= sum(entropy(g) for g in edge_laws.values()) + 4e-13
        finite_checks += 1
    print(f"PASS: {information_checks} heterogeneous data-processing checks; {repair_checks} exact repair couplings; {finite_checks} literal-array entropy identities")


if __name__ == "__main__":
    main()
