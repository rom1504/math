"""Finite exact-state replay for the designer Gibbs nonprojectivity note."""
from itertools import combinations
import math
import numpy as np


def catalogue(n):
    edges = list(combinations(range(n), 2))
    masks = np.arange(1 << len(edges), dtype=np.uint64)
    signs = 1 - 2 * ((masks[:, None] >> np.arange(len(edges), dtype=np.uint64)) & 1).astype(np.int8)
    spin_masks = np.arange(1 << (n - 1), dtype=np.uint64)
    spins = np.ones((len(spin_masks), n), dtype=np.int8)
    spins[:, 1:] = 1 - 2 * ((spin_masks[:, None] >> np.arange(n - 1, dtype=np.uint64)) & 1).astype(np.int8)
    tests = np.array([spins[:, i] * spins[:, j] for i, j in edges], dtype=np.int16)
    caps = np.max(np.abs(signs.astype(np.int16) @ tests), axis=1)
    return edges, masks, caps


def law(caps, n, beta):
    logs = -beta * math.sqrt(n) * caps
    maximum = float(logs.max())
    partition = maximum + math.log(float(np.exp(logs - maximum).sum()))
    return np.exp(logs - partition), partition


def main():
    data = {n: catalogue(n) for n in range(2, 7)}
    checks = 0
    worst = 0.0
    for n in range(3, 7):
        edges, masks, caps = data[n]
        for beta in (0.1, 1.0, 4.0):
            parent, _ = law(caps, n, beta)
            for k in range(2, n):
                child_edges, _, child_caps = data[k]
                positions = [edges.index(edge) for edge in child_edges]
                code = np.zeros(len(masks), dtype=np.uint64)
                for j, pos in enumerate(positions):
                    code |= ((masks >> pos) & 1) << j
                marginal = np.bincount(code.astype(np.int64), weights=parent, minlength=len(child_caps))
                reference, log_partition = law(child_caps, k, beta)
                active = marginal > 0
                entropy = -float(np.dot(marginal[active], np.log(marginal[active])))
                mean_cap = float(np.dot(marginal, child_caps))
                # Use exact log-density expression, avoiding underflow in reference.
                log_reference = -beta * math.sqrt(k) * child_caps - log_partition
                divergence = float(np.dot(marginal[active], np.log(marginal[active]) - log_reference[active]))
                value = entropy / k**2 - beta * mean_cap / k**1.5
                discrepancy = abs(log_partition / k**2 - value - divergence / k**2)
                worst = max(worst, discrepancy)
                assert discrepancy < 2e-12
                assert divergence >= -2e-12
                for threshold in sorted(set(child_caps)):
                    probability = float(reference[child_caps >= threshold].sum())
                    exponent = len(child_edges) * math.log(2) - beta * math.sqrt(k) * (threshold - child_caps.min())
                    assert probability <= math.exp(min(0.0, exponent)) + 2e-12
                    checks += 1
    print(f"PASS: {checks} one-atom tail checks; all marginal KL identities; maximum error {worst:.3g}")
    print(f"Concrete asymptotic gap (1-log(2))/2 = {(1-math.log(2))/2:.15f}")


if __name__ == "__main__":
    main()
