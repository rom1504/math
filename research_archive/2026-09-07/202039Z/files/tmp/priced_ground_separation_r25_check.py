#!/usr/bin/env python3
"""Exact checks for the Wave 25 priced-ground separation memo."""

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


def grounds(A):
    n = len(A)
    parent_m = SC.PC.matrix_norm(A)[0]
    child_m = [
        SC.PC.matrix_norm(SC.PC.principal(A, (i,))[0])[0]
        for i in range(n)
    ]
    out = []
    for orientation in (-1, 1):
        for tail in itertools.product((-1, 1), repeat=n - 1):
            x = (1,) + tail
            if SC.score(A, orientation, x) != parent_m:
                continue
            W = tuple(
                tuple(
                    0 if i == j else orientation * A[i][j] * x[i] * x[j]
                    for j in range(n)
                )
                for i in range(n)
            )
            rows = tuple(sum(W[i]) for i in range(n))
            assert min(rows) >= 0
            assert sum(rows) == 2 * parent_m
            coverage = Fraction(
                sum(
                    SC.score(A, orientation, x, (i,)) == child_m[i]
                    for i in range(n)
                ),
                n,
            )
            row_square = sum(r * r for r in rows)

            # R2 = ||Ax||^2, independently of orientation.
            Ax = tuple(sum(A[i][j] * x[j] for j in range(n)) for i in range(n))
            assert row_square == sum(value * value for value in Ax)
            out.append((orientation, x, W, rows, coverage, row_square))
    return parent_m, child_m, out


EXPECTED = {
    "A6": Counter({(Fraction(5, 6), 30): 12}),
    "A8": Counter({(Fraction(3, 8), 64): 8}),
    "A9": Counter(
        {
            (Fraction(1, 9), 80): 2,
            (Fraction(2, 9), 96): 12,
            (Fraction(1, 3), 112): 4,
            (Fraction(2, 9), 112): 3,
            (Fraction(1, 3), 128): 4,
        }
    ),
}


def verify_exchanges(A, parent_m, child_m, states):
    n = len(A)
    for g in states:
        sigma, x, W, rows_g, _, r2_g = g
        for h in states:
            tau, y, _, rows_h, _, r2_h = h
            epsilon = sigma * tau
            z = tuple(x[i] * y[i] for i in range(n))
            cross_degree = tuple(
                sum(W[i][j] for j in range(n) if z[j] != z[i])
                for i in range(n)
            )
            shore_weight = sum(
                W[i][j]
                for i in range(n)
                for j in range(i + 1, n)
                if z[i] != z[j]
            )
            assert shore_weight == (0 if epsilon == 1 else parent_m)
            predicted_rows = tuple(
                epsilon * (rows_g[i] - 2 * cross_degree[i])
                for i in range(n)
            )
            assert predicted_rows == rows_h
            assert r2_h - r2_g == 4 * sum(
                a_i * a_i - r_i * a_i
                for a_i, r_i in zip(cross_degree, rows_g)
            )

            for omitted in range(n):
                selected = [i for i in range(n) if i != omitted]
                energy_g = sum(
                    W[i][j]
                    for ii, i in enumerate(selected)
                    for j in selected[ii + 1 :]
                )
                restricted_shore = sum(
                    W[i][j]
                    for ii, i in enumerate(selected)
                    for j in selected[ii + 1 :]
                    if z[i] != z[j]
                )
                payoff_g = 2 * SC.score(A, sigma, x, (omitted,))
                payoff_h = 2 * SC.score(A, tau, y, (omitted,))
                assert payoff_g == 2 * energy_g
                loss_difference = payoff_g - payoff_h
                predicted = (
                    2 * (1 - epsilon) * energy_g
                    + 4 * epsilon * restricted_shore
                )
                assert loss_difference == predicted


def positive_part(value):
    return max(Fraction(0), value)


def dual_value(types, budget):
    actions = list(types) + [(Fraction(0), Fraction(0))]
    lines = [(u, u * (budget - cost)) for u, cost in actions]
    candidates = {Fraction(0)}
    for index, (a0, a1) in enumerate(lines):
        for b0, b1 in lines[index + 1 :]:
            if a1 == b1:
                continue
            theta = (b0 - a0) / (a1 - b1)
            if theta >= 0:
                candidates.add(theta)
    return min(max(a0 + theta * a1 for a0, a1 in lines) for theta in candidates)


def pair_value(types, budget):
    answer = max((u for u, cost in types if cost <= budget), default=Fraction(0))
    for u_a, c_a in types:
        for u_b, c_b in types:
            if not (c_a < budget < c_b) or u_a == 0 or u_b == 0:
                continue
            inverse = (
                Fraction(c_b - budget, c_b - c_a) / u_a
                + Fraction(budget - c_a, c_b - c_a) / u_b
            )
            answer = max(answer, 1 / inverse)
    return answer


def verify_price_edge_cases():
    """Audit dummy, positive part, endpoints, and zero-coverage actions."""
    zero = Fraction(0)
    one_third = Fraction(1, 3)
    one_half = Fraction(1, 2)

    # A zero-coverage action below budget cannot subsidize an action above it.
    types = Counter({(zero, zero): 1, (one_half, Fraction(10)): 1})
    assert dual_value(types, Fraction(5)) == zero
    assert pair_value(types, Fraction(5)) == zero

    # At the lowest positive cost, only columns exactly at that cost survive.
    types = Counter(
        {
            (zero, Fraction(-10)): 1,
            (one_third, Fraction(10)): 1,
            (one_half, Fraction(20)): 1,
        }
    )
    assert dual_value(types, Fraction(9)) == zero
    assert pair_value(types, Fraction(9)) == zero
    assert dual_value(types, Fraction(10)) == one_third
    assert pair_value(types, Fraction(10)) == one_third
    assert dual_value(types, Fraction(20)) == one_half
    assert pair_value(types, Fraction(20)) == one_half

    # Strict crossing interpolation and its two singleton endpoints agree.
    crossing = Counter(
        {(one_third, Fraction(10)): 1, (one_half, Fraction(20)): 1}
    )
    expected_midpoint = Fraction(2, 5)
    assert dual_value(crossing, Fraction(15)) == expected_midpoint
    assert pair_value(crossing, Fraction(15)) == expected_midpoint

    # Including the dummy is exactly taking the positive part of real lines.
    price = Fraction(1, 2)
    budget = Fraction(5)
    real_lines = [u * (1 + price * (budget - cost)) for u, cost in crossing]
    with_dummy = max([zero] + real_lines)
    positive_parts = max(positive_part(value) for value in real_lines)
    assert with_dummy == positive_parts


def main():
    verify_price_edge_cases()
    print("dummy, positive-part, endpoint, and zero-coverage audits pass")

    matrices = {name: data[0] for name, data in SC.CASES.items()}
    all_types = {}
    for name, A in matrices.items():
        parent_m, child_m, states = grounds(A)
        types = Counter((state[4], state[5]) for state in states)
        assert types == EXPECTED[name]
        verify_exchanges(A, parent_m, child_m, states)
        all_types[name] = types
        print(f"{name}: {len(states)} grounds; all pair exchanges pass")

    assert dual_value(all_types["A6"], Fraction(29)) == 0
    assert dual_value(all_types["A6"], Fraction(30)) == Fraction(5, 6)
    assert dual_value(all_types["A8"], Fraction(63)) == 0
    assert dual_value(all_types["A8"], Fraction(64)) == Fraction(3, 8)

    for budget, expected, price in (
        (Fraction(80), Fraction(1, 9), None),
        (Fraction(88), Fraction(4, 27), Fraction(1, 24)),
        (Fraction(96), Fraction(2, 9), None),
        (Fraction(104), Fraction(4, 15), Fraction(1, 40)),
        (Fraction(112), Fraction(1, 3), Fraction(0)),
    ):
        dual = dual_value(all_types["A9"], budget)
        primal = pair_value(all_types["A9"], budget)
        assert dual == primal == expected
        if price is not None:
            priced = max(
                positive_part(u * (1 + price * (budget - cost)))
                for u, cost in all_types["A9"]
            )
            assert priced == expected
        print(f"A9 C={budget}: Z*={expected}")

    print("priced exact-ground separation certificates passed")


if __name__ == "__main__":
    main()
