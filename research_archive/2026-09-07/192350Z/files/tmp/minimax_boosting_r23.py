#!/usr/bin/env python3
"""Exact A9 minimax/cover certificates for the Wave 23 boosting memo."""

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


def tolerance_patterns(one_copy_tolerance):
    """Patterns use M=Q/2, so ledger tolerance is twice the argument."""
    A = SC.PC.A9
    n = len(A)
    child_q = [
        SC.PC.matrix_norm(SC.PC.principal(A, (i,))[0])[0]
        for i in range(n)
    ]
    patterns = set()
    for orientation in (-1, 1):
        for tail in itertools.product((-1, 1), repeat=n - 1):
            x = (1,) + tail
            pattern = tuple(
                i
                for i in range(n)
                if child_q[i] - SC.score(A, orientation, x, (i,))
                <= one_copy_tolerance
            )
            if pattern:
                patterns.add(pattern)
    return patterns


def regular_parent_patterns(ledger_tolerance, row_cap):
    """Coverage patterns of A9 parent grounds with the stated row cap."""
    A = SC.PC.A9
    n = len(A)
    parent_q = SC.PC.matrix_norm(A)[0]
    child_q = [
        SC.PC.matrix_norm(SC.PC.principal(A, (i,))[0])[0]
        for i in range(n)
    ]
    patterns = set()
    cap_counts = []
    for orientation in (-1, 1):
        for tail in itertools.product((-1, 1), repeat=n - 1):
            x = (1,) + tail
            if SC.score(A, orientation, x) != parent_q:
                continue
            rows = tuple(
                orientation
                * x[i]
                * sum(A[i][j] * x[j] for j in range(n) if j != i)
                for i in range(n)
            )
            assert min(rows) >= 0
            assert sum(rows) == 2 * parent_q
            cap_counts.append(max(rows))
            if max(rows) > row_cap:
                continue
            pattern = tuple(
                i
                for i in range(n)
                if 2 * (child_q[i] - SC.score(A, orientation, x, (i,)))
                <= ledger_tolerance
            )
            if pattern:
                patterns.add(pattern)
    assert len(cap_counts) == 25
    return patterns, cap_counts


CASES = {
    0: {
        "count": 22,
        "alpha": Fraction(4, 13),
        "dual": tuple(Fraction(x, 13) for x in (1, 2, 1, 1, 2, 2, 2, 1, 1)),
        "primal": {
            (3, 6, 8): Fraction(1, 26),
            (4, 7, 8): Fraction(1, 26),
            (4, 6): Fraction(3, 26),
            (1, 5): Fraction(3, 26),
            (0, 4, 8): Fraction(4, 26),
            (1, 3, 8): Fraction(2, 26),
            (2, 3, 5): Fraction(5, 26),
            (1, 2, 7): Fraction(3, 26),
            (0, 6, 7): Fraction(4, 26),
        },
        "cover": ((4, 7, 8), (2, 3, 5), (1, 2, 7), (0, 6, 7)),
    },
    2: {
        "count": 46,
        "alpha": Fraction(5, 8),
        "dual": tuple(Fraction(x, 8) for x in (1, 1, 1, 0, 2, 1, 1, 1, 0)),
        "primal": {
            (3, 4, 5, 6, 7, 8): Fraction(1, 8),
            (1, 3, 4, 5, 7, 8): Fraction(1, 8),
            (0, 3, 4, 5, 6, 8): Fraction(1, 8),
            (0, 2, 3, 4, 5, 8): Fraction(1, 8),
            (1, 2, 3, 4, 5): Fraction(1, 8),
            (0, 1, 2, 6, 7, 8): Fraction(3, 8),
        },
        "cover": ((1, 3, 5, 6, 7), (0, 2, 3, 4, 5, 8)),
    },
    4: {
        "count": 78,
        "alpha": Fraction(1),
        "dual": (Fraction(1),) + (Fraction(0),) * 8,
        "primal": {tuple(range(9)): Fraction(1)},
        "cover": (tuple(range(9)),),
    },
}


REGULAR_CASES = {
    (4, 0): {"alpha": Fraction(0)},
    (4, 4): {"alpha": Fraction(0)},
    (4, 8): {
        "alpha": Fraction(1),
        "dual": (Fraction(1),) + (Fraction(0),) * 8,
        "primal": {tuple(range(9)): Fraction(1)},
        "cover": (tuple(range(9)),),
    },
    (6, 0): {
        "alpha": Fraction(2, 7),
        "dual": tuple(Fraction(x, 7) for x in (1, 1, 1, 0, 1, 1, 1, 1, 0)),
        "primal": {
            (0, 4, 8): Fraction(2, 7),
            (1, 6): Fraction(1, 7),
            (1, 7): Fraction(1, 7),
            (2, 3, 5): Fraction(2, 7),
            (6, 7): Fraction(1, 7),
        },
        "cover": ((0, 4, 8), (1, 6), (1, 7), (2, 3, 5)),
    },
    (6, 4): {
        "alpha": Fraction(1, 2),
        "dual": tuple(Fraction(x, 6) for x in (1, 0, 1, 0, 1, 1, 0, 2, 0)),
        "primal": {
            (0, 3, 6, 7, 8): Fraction(1, 2),
            (1, 2, 3, 4, 5): Fraction(1, 2),
        },
        "cover": ((0, 1, 6, 7, 8), (1, 2, 3, 4, 5)),
    },
    (6, 8): {
        "alpha": Fraction(1),
        "dual": (Fraction(1),) + (Fraction(0),) * 8,
        "primal": {tuple(range(9)): Fraction(1)},
        "cover": (tuple(range(9)),),
    },
    (8, 0): CASES[0],
    (8, 4): CASES[2],
    (8, 8): CASES[4],
}


def covers(patterns):
    return set().union(*(set(pattern) for pattern in patterns)) == set(range(9))


def verify_case(tolerance, data):
    patterns = tolerance_patterns(tolerance)
    assert len(patterns) == data["count"]
    alpha = data["alpha"]
    dual = data["dual"]
    primal = data["primal"]

    # Dual selector weighting: every hypothesis covers at most alpha.
    assert sum(dual) == 1
    assert all(sum(dual[i] for i in pattern) <= alpha for pattern in patterns)

    # Primal common prior: every selector is covered with probability >= alpha.
    assert sum(primal.values()) == 1
    assert set(primal).issubset(patterns)
    selector_mass = tuple(
        sum(weight for pattern, weight in primal.items() if i in pattern)
        for i in range(9)
    )
    assert min(selector_mass) == alpha

    # Exact integral cover and exhaustive exclusion of every smaller cover.
    cover = data["cover"]
    assert set(cover).issubset(patterns)
    assert covers(cover)
    pattern_list = sorted(patterns)
    for size in range(1, len(cover)):
        assert not any(covers(choice) for choice in itertools.combinations(pattern_list, size))

    print(
        f"ledger_t={2*tolerance}: patterns={len(patterns)}, "
        f"alpha={alpha}, minimum_cover={len(cover)}"
    )


def verify_game(patterns, data):
    alpha = data["alpha"]
    if alpha == 0:
        assert set().union(*(set(pattern) for pattern in patterns)) != set(range(9))
        return None

    dual = data["dual"]
    primal = data["primal"]
    assert sum(dual) == 1
    assert all(sum(dual[i] for i in pattern) <= alpha for pattern in patterns)
    assert sum(primal.values()) == 1
    assert set(primal).issubset(patterns)
    selector_mass = tuple(
        sum(weight for pattern, weight in primal.items() if i in pattern)
        for i in range(9)
    )
    assert min(selector_mass) == alpha

    cover = data["cover"]
    assert set(cover).issubset(patterns)
    assert covers(cover)
    pattern_list = sorted(patterns)
    for size in range(1, len(cover)):
        assert not any(covers(choice) for choice in itertools.combinations(pattern_list, size))
    return len(cover)


def verify_regular_cases():
    all_caps = None
    for (row_cap, tolerance), data in REGULAR_CASES.items():
        patterns, cap_counts = regular_parent_patterns(tolerance, row_cap)
        all_caps = cap_counts
        cover_size = verify_game(patterns, data)
        print(
            f"regular_B={row_cap}, ledger_t={tolerance}: "
            f"patterns={len(patterns)}, alpha={data['alpha']}, "
            f"minimum_cover={cover_size}"
        )
    assert all_caps is not None
    assert {cap: all_caps.count(cap) for cap in set(all_caps)} == {4: 2, 6: 16, 8: 7}


def main():
    for tolerance, data in CASES.items():
        verify_case(tolerance, data)
    verify_regular_cases()
    print("A9 minimax boosting certificates passed")


if __name__ == "__main__":
    main()
