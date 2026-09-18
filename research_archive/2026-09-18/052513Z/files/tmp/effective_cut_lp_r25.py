#!/usr/bin/env python3
"""Exact finite audit of the all-cut effective-loss LP from Wave 25."""

from collections import Counter
from fractions import Fraction
import importlib.util
import itertools
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "selector_codebook_r22", ROOT / "tmp" / "selector_codebook_r22.py"
)
SC = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(SC)


def all_cut_types(A, threshold=Fraction(0)):
    """Return exact (effective-loss coverage, row-square, deficit) types."""
    n = len(A)
    p2 = Fraction(n - 2, n)  # one deletion: (n-1)_2/(n)_2
    parent_q = 2 * SC.PC.matrix_norm(A)[0]
    child_q = [
        2 * SC.PC.matrix_norm(SC.PC.principal(A, (i,))[0])[0]
        for i in range(n)
    ]
    types = []
    detailed = []
    for orientation in (-1, 1):
        for tail in itertools.product((-1, 1), repeat=n - 1):
            x = (1,) + tail
            parent_energy = 2 * SC.score(A, orientation, x)
            deficit = parent_q - parent_energy
            assert deficit >= 0
            rows = tuple(
                orientation
                * x[i]
                * sum(A[i][j] * x[j] for j in range(n) if j != i)
                for i in range(n)
            )
            assert sum(rows) == parent_energy
            r2 = sum(r * r for r in rows)
            effective = []
            for i in range(n):
                child_energy = 2 * SC.score(A, orientation, x, (i,))
                loss = child_q[i] - child_energy
                assert loss >= 0
                effective.append(Fraction(loss) - p2 * deficit)
            coverage = Fraction(sum(v <= threshold for v in effective), n)
            types.append((coverage, r2))
            detailed.append((coverage, r2, deficit, tuple(effective)))
    return Counter(types), detailed


def ground_types(A):
    n = len(A)
    parent_q = 2 * SC.PC.matrix_norm(A)[0]
    out = []
    for orientation in (-1, 1):
        for tail in itertools.product((-1, 1), repeat=n - 1):
            x = (1,) + tail
            parent_energy = 2 * SC.score(A, orientation, x)
            if parent_energy != parent_q:
                continue
            rows = tuple(
                orientation
                * x[i]
                * sum(A[i][j] * x[j] for j in range(n) if j != i)
                for i in range(n)
            )
            child_q = [
                2 * SC.PC.matrix_norm(SC.PC.principal(A, (i,))[0])[0]
                for i in range(n)
            ]
            coverage = Fraction(
                sum(2 * SC.score(A, orientation, x, (i,)) == child_q[i] for i in range(n)),
                n,
            )
            out.append((coverage, sum(r * r for r in rows)))
    return Counter(out)


def dual_value(types, budget):
    """Exact value min_theta>=0 max(0,u[1+theta(C-c)])."""
    actions = list(types) + [(Fraction(0), Fraction(0))]
    lines = [(u, u * (budget - cost)) for u, cost in actions]
    candidates = {Fraction(0)}
    for a, (ia, sa) in enumerate(lines):
        for ib, sb in lines[a + 1 :]:
            if sa == sb:
                continue
            theta = (ib - ia) / (sa - sb)
            if theta >= 0:
                candidates.add(theta)
    return min(max(i + theta * s for i, s in lines) for theta in candidates)


def nondominated(types):
    points = sorted(set(types))
    return [
        (u, c)
        for u, c in points
        if u > 0
        and not any(u2 >= u and c2 <= c and (u2, c2) != (u, c) for u2, c2 in points)
    ]


def main():
    expected_fronts = {
        "A6": [(Fraction(5, 6), 30)],
        "A8": [(Fraction(3, 8), 32), (Fraction(1, 2), 64)],
        "A9": [
            (Fraction(1, 9), 24),
            (Fraction(2, 9), 32),
            (Fraction(1, 3), 48),
            (Fraction(4, 9), 112),
        ],
    }
    for name, data in SC.CASES.items():
        A = data[0]
        types, detailed = all_cut_types(A)
        grounds = ground_types(A)
        assert sum(types.values()) == 2 ** len(A)
        assert all(types[p] >= multiplicity for p, multiplicity in grounds.items())
        front = nondominated(types)
        assert front == expected_fronts[name]
        costs = sorted({c for u, c in front if u > 0})
        probes = sorted(
            set(
                costs
                + [Fraction(costs[i] + costs[i + 1], 2) for i in range(len(costs) - 1)]
            )
        )
        strict = False
        rows = []
        for budget in probes:
            za = dual_value(types, budget)
            zg = dual_value(grounds, budget)
            assert za >= zg
            strict = strict or za > zg
            rows.append((budget, za, zg))
        print(name, "all-cut nondominated", front)
        representatives = []
        for point in front:
            representatives.append(next((deficit, eff) for u, c, deficit, eff in detailed if (u, c) == point))
        print(name, "representative (deficit,effective losses)", representatives)
        print(name, "LP probes", rows)
        print(name, "strict improvement", strict)
        # The average effective loss is independent of the full cut.
        means = {sum(eff, Fraction(0)) / len(A) for _, _, _, eff in detailed}
        assert len(means) == 1
        q = 2 * SC.PC.matrix_norm(A)[0]
        mean_child_q = Fraction(
            sum(
                2 * SC.PC.matrix_norm(SC.PC.principal(A, (i,))[0])[0]
                for i in range(len(A))
            ),
            len(A),
        )
        assert means == {mean_child_q - Fraction(len(A) - 2, len(A)) * q}
        print(name, "cut-independent mean effective loss", next(iter(means)))
    print("PASS: all-cut effective-loss identities and LP domination")


if __name__ == "__main__":
    main()
