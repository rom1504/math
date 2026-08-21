#!/usr/bin/env python3
"""Exact finite checks for the 2026-08-21 blank-slate attack.

The script uses only integer arithmetic.  It verifies

* the two order-eight radial-moment certificates;
* the order-eight radius-two local-minimum counterexample; and
* small instances of the scalable one-edge-local family; and
* exhaustive tiny cases of the degree-two rectangular inequality.

The first spin is fixed to +1 because a global spin flip leaves every
quadratic energy unchanged.
"""

from __future__ import annotations

import argparse
import itertools
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path


def edges(n: int) -> list[tuple[int, int]]:
    return list(itertools.combinations(range(n), 2))


def gauge_matrix_from_mask(n: int, mask: int) -> dict[tuple[int, int], int]:
    """First row positive; mask bits encode negative internal edges."""
    ans = {edge: 1 for edge in edges(n)}
    internal = list(itertools.combinations(range(1, n), 2))
    for bit, edge in enumerate(internal):
        if (mask >> bit) & 1:
            ans[edge] = -1
    return ans


def energies(n: int, signing: dict[tuple[int, int], int]) -> list[int]:
    ans: list[int] = []
    edge_list = edges(n)
    for tail in itertools.product((-1, 1), repeat=n - 1):
        spin = (1,) + tail
        ans.append(
            sum(signing[i, j] * spin[i] * spin[j] for i, j in edge_list)
        )
    return ans


def cap(n: int, signing: dict[tuple[int, int], int]) -> int:
    return max(map(abs, energies(n, signing)))


def flip(
    signing: dict[tuple[int, int], int],
    changed: tuple[tuple[int, int], ...] | list[tuple[int, int]],
) -> dict[tuple[int, int], int]:
    ans = signing.copy()
    for edge in changed:
        ans[edge] *= -1
    return ans


def even_moments(values: list[int], highest: int = 8) -> dict[str, int]:
    return {
        f"m{degree}": sum(value**degree for value in values) // len(values)
        for degree in range(2, highest + 1, 2)
    }


def eval_z_polynomial(moments: dict[str, int], coefficients: dict[int, int]) -> int:
    """Evaluate E sum_power coefficients[power] (q^2)^power."""
    total = coefficients.get(0, 0)
    for power, coefficient in coefficients.items():
        if power:
            total += coefficient * moments[f"m{2 * power}"]
    return total


def radial_checks() -> dict[str, object]:
    certificates = {
        "central_cap_8": {1: 64, 2: -1},
        "central_cap_10": {2: 100, 3: -1},
        "central_even_lattice_cap_12": {1: 9216, 2: -2944, 3: 164, 4: -1},
    }
    records: dict[str, object] = {}
    for mask in (6875, 6887):
        vals = energies(8, gauge_matrix_from_mask(8, mask))
        moments = even_moments(vals)
        records[str(mask)] = {
            "cap": max(map(abs, vals)),
            "histogram": dict(sorted(Counter(vals).items())),
            "moments": moments,
            "certificate_expectations": {
                name: eval_z_polynomial(moments, polynomial)
                for name, polynomial in certificates.items()
            },
        }

    assert records["6875"]["cap"] == 14
    assert records["6887"]["cap"] == 12
    assert records["6875"]["moments"] == {
        "m2": 28,
        "m4": 2152,
        "m6": 263008,
        "m8": 37459072,
    }
    assert records["6887"]["moments"] == {
        "m2": 28,
        "m4": 2152,
        "m6": 228448,
        "m8": 27782272,
    }
    assert records["6875"]["certificate_expectations"] == {
        "central_cap_8": -360,
        "central_cap_10": -47808,
        "central_even_lattice_cap_12": -403200,
    }
    assert records["6887"]["certificate_expectations"] == {
        "central_cap_8": -360,
        "central_cap_10": -13248,
        "central_even_lattice_cap_12": 3605760,
    }
    return records


def radius_two_counterexample() -> dict[str, object]:
    n = 8
    negative_one_based = {
        (1, 6),
        (2, 6),
        (2, 7),
        (3, 5),
        (3, 7),
        (4, 5),
        (4, 8),
        (6, 8),
    }
    negative = {(i - 1, j - 1) for i, j in negative_one_based}
    signing = {edge: (-1 if edge in negative else 1) for edge in edges(n)}
    base_values = energies(n, signing)
    neighbor_caps: Counter[int] = Counter()
    for radius in (1, 2):
        for changed in itertools.combinations(edges(n), radius):
            neighbor_caps[cap(n, flip(signing, changed))] += 1

    improving_one_based = ((1, 3), (2, 5), (3, 7))
    improving = tuple((i - 1, j - 1) for i, j in improving_one_based)
    improved_cap = cap(n, flip(signing, improving))

    assert max(map(abs, base_values)) == 12
    assert neighbor_caps == Counter({14: 197, 16: 107, 12: 102})
    assert improved_cap == 10
    return {
        "negative_edges_one_based": sorted(negative_one_based),
        "histogram": dict(sorted(Counter(base_values).items())),
        "single_and_double_flip_cap_counts": dict(sorted(neighbor_caps.items())),
        "improving_triple_one_based": improving_one_based,
        "improving_triple_cap": improved_cap,
    }


def scalable_local_family() -> dict[str, object]:
    records = []
    for n in (6, 8, 10, 12):
        m = n // 2
        right = set(range(m, n))
        signing = {
            edge: (-1 if edge[0] in right and edge[1] in right else 1)
            for edge in edges(n)
        }
        q = cap(n, signing)
        single_caps = Counter(
            cap(n, flip(signing, (edge,))) for edge in edges(n)
        )
        assert q == m * m
        assert min(single_caps) >= q
        records.append(
            {
                "n": n,
                "cap": q,
                "single_flip_cap_counts": dict(sorted(single_caps.items())),
            }
        )
    return {"formula": "Q(A_2m)=m^2", "finite_checks": records}


def rectangular_fourier_checks() -> dict[str, object]:
    """Exhaustively verify the degree-two rectangular inequality at tiny sizes."""
    records = []
    for row_count, column_count in ((1, 2), (2, 2), (2, 3), (3, 3), (2, 4)):
        vectors = list(itertools.product((-1, 1), repeat=column_count))
        mu = Fraction(
            sum(abs(sum(vector)) for vector in vectors), len(vectors)
        )
        eta = Fraction(
            _binomial(column_count - 2, (column_count - 2) // 2),
            1 << (column_count - 2),
        )
        checked = 0
        minimum_slack: Fraction | None = None
        for flat in itertools.product((-1, 1), repeat=row_count * column_count):
            rows = tuple(
                flat[i * column_count : (i + 1) * column_count]
                for i in range(row_count)
            )
            norm = max(
                sum(abs(sum(row[j] * vector[j] for j in range(column_count))) for row in rows)
                for vector in vectors
            )
            correlation_square_sum = sum(
                sum(row[p] * row[q] for row in rows) ** 2
                for p, q in itertools.combinations(range(column_count), 2)
            )
            rhs = row_count * mu + (
                eta * eta * correlation_square_sum / (row_count * mu)
            )
            slack = Fraction(norm) - rhs
            assert slack >= 0
            minimum_slack = slack if minimum_slack is None else min(minimum_slack, slack)
            checked += 1
        records.append(
            {
                "rows": row_count,
                "columns": column_count,
                "matrices_checked": checked,
                "minimum_slack": str(minimum_slack),
            }
        )
    return {"finite_exhaustive_checks": records}


def augmented_cut_code_masks(n: int) -> list[int]:
    edge_list = edges(n)
    all_ones = (1 << len(edge_list)) - 1
    cuts = set()
    for tail in itertools.product((-1, 1), repeat=n - 1):
        spin = (1,) + tail
        mask = 0
        for bit, (i, j) in enumerate(edge_list):
            if spin[i] * spin[j] == -1:
                mask |= 1 << bit
        cuts.add(mask)
    augmented = cuts | {mask ^ all_ones for mask in cuts}
    assert len(augmented) == 1 << n
    return sorted(augmented)


def bonferroni_checks() -> dict[str, object]:
    n = 5
    edge_count = len(edges(n))
    code = augmented_cut_code_masks(n)
    records: dict[str, object] = {}
    for radius in (2, 3):
        multiplicities = []
        for word in range(1 << edge_count):
            multiplicities.append(
                sum(_popcount(word ^ center) <= radius for center in code)
            )
        histogram = Counter(multiplicities)
        covered = sum(count for z, count in histogram.items() if z)
        highest = max(histogram)
        partials = {}
        running = 0
        for order in range(1, highest + 1):
            intersection_sum = sum(
                count * _binomial(z, order) for z, count in histogram.items()
            )
            running += (1 if order % 2 else -1) * intersection_sum
            partials[order] = running
        records[str(radius)] = {
            "multiplicity_histogram": dict(sorted(histogram.items())),
            "covered_words": covered,
            "bonferroni_partial_sums": partials,
        }

    assert records["2"] == {
        "multiplicity_histogram": {0: 192, 1: 352, 3: 480},
        "covered_words": 832,
        "bonferroni_partial_sums": {1: 1792, 2: 352, 3: 832},
    }
    assert records["3"]["covered_words"] == 1024
    assert records["3"]["multiplicity_histogram"] == {1: 32, 3: 480, 7: 320, 10: 192}
    assert records["3"]["bonferroni_partial_sums"] == {
        1: 5632,
        2: -11168,
        3: 23552,
        4: -27968,
        5: 27136,
        6: -15424,
        7: 7936,
        8: -704,
        9: 1216,
        10: 1024,
    }
    return records


def _binomial(n: int, k: int) -> int:
    if k < 0 or k > n:
        return 0
    answer = 1
    for j in range(1, k + 1):
        answer = answer * (n - j + 1) // j
    return answer


def _popcount(value: int) -> int:
    return bin(value).count("1")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = {
        "radial_certificates": radial_checks(),
        "radius_two_counterexample": radius_two_counterexample(),
        "scalable_one_edge_local_family": scalable_local_family(),
        "rectangular_degree_two_inequality": rectangular_fourier_checks(),
        "bonferroni_union_checks": bonferroni_checks(),
    }
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
