"""Finite diagnostic only: actual archived children and one H2 lift.

Run with the repository virtualenv. No claim of asymptotic minimality or
of a universal tensor inequality is made by these finite computations.
"""
import json
from pathlib import Path
import numpy as np


def cube(n):
    bits = (np.arange(1 << (n - 1), dtype=np.uint64)[:, None]
            >> np.arange(n - 1, dtype=np.uint64)) & 1
    return np.column_stack((np.ones(len(bits), dtype=np.int64),
                            1 - 2 * bits.astype(np.int64)))


def main():
    results = []
    for n in range(3, 11):
        source = Path(f"computations/results/exact_m{n}.json")
        if not source.exists():
            continue
        a = np.array(json.loads(source.read_text())["matrix"], dtype=np.int64)
        x = cube(n)
        h = np.einsum("bi,ij,bj->b", x, a, x) // 2
        beta = int(np.abs(x @ a).sum(axis=1).max())
        # Pair energies use both y and -y, hence absolute bridge exactly.
        cross = x @ a @ x.T
        lifted = int(np.max(np.abs(h[:, None] - h[None, :]) + np.abs(cross)))
        q = int(np.abs(h).max())
        results.append(dict(n=n, source=str(source), Q=q, beta=beta,
                            hollow_tensor_cap=lifted,
                            tensor_ratio=lifted / q,
                            desired_ratio=2 * np.sqrt(2),
                            additive_excess=lifted - 2 * np.sqrt(2) * q))
    Path("computations/results/flatify_construct_2026_09_07_tensor_diagnostic.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
