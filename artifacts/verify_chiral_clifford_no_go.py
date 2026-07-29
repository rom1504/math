#!/usr/bin/env python3
"""Exact finite audit of the inherited-chiral 4-lift ansatz.

The script:
  * enumerates every compatible (P,R,K,E) Clifford mask;
  * proves 912/1008 have a uniform-fibre lower bound above 320;
  * quotients the 96 survivors into four coordinate-permutation orbits;
  * finds one two-pattern Boolean witness for each orbit;
  * uses exact subset-sum dynamic programming to prove that no compatible
    choice of the remaining 24 diagonal signs can lower that witness to 320.

No optimization solver is used in the final certificate.
"""

import json
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations, product
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
TARGET = 320


def exact_rank(rows):
    """Rank over Q, used only on matrices with at most 16 columns."""
    matrix = [
        [Fraction(int(entry)) for entry in row]
        for row in rows
    ]
    if not matrix:
        return 0
    row = 0
    for column in range(len(matrix[0])):
        pivot = next(
            (
                index
                for index in range(row, len(matrix))
                if matrix[index][column]
            ),
            None,
        )
        if pivot is None:
            continue
        matrix[row], matrix[pivot] = matrix[pivot], matrix[row]
        pivot_value = matrix[row][column]
        matrix[row] = [entry / pivot_value for entry in matrix[row]]
        for index in range(len(matrix)):
            if index == row or not matrix[index][column]:
                continue
            multiplier = matrix[index][column]
            matrix[index] = [
                old - multiplier * new
                for old, new in zip(matrix[index], matrix[row])
            ]
        row += 1
    return row


def certify_full_fibre_algebra(base):
    """Prove that two cross-fibre blocks generate M_4(R)."""
    first = base[0:4, 4:8]
    second = base[0:4, 8:12]
    identity = np.eye(4, dtype=np.int64)
    basis = []
    queue = [identity]
    while queue and len(basis) < 16:
        candidate = queue.pop(0)
        trial = basis + [candidate]
        if exact_rank([matrix.reshape(-1) for matrix in trial]) == len(basis):
            continue
        basis.append(candidate)
        queue.extend((candidate @ first, candidate @ second))
    assert len(basis) == 16


def exact_characteristic_polynomial(matrix):
    """Faddeev--LeVerrier over Z."""
    size = len(matrix)
    identity = np.eye(size, dtype=object)
    matrix = matrix.astype(object)
    auxiliary = identity.copy()
    coefficients = [1]
    for degree in range(1, size + 1):
        product_matrix = matrix @ auxiliary
        trace = sum(product_matrix[index, index] for index in range(size))
        assert trace % degree == 0
        coefficient = -trace // degree
        coefficients.append(int(coefficient))
        auxiliary = product_matrix + coefficient * identity
    return coefficients


def polynomial_remainder(dividend, divisor):
    """Remainder over Q for leading-coefficient-first polynomials."""
    dividend = [Fraction(coefficient) for coefficient in dividend]
    divisor = [Fraction(coefficient) for coefficient in divisor]
    while dividend and dividend[0] == 0:
        dividend.pop(0)
    while divisor and divisor[0] == 0:
        divisor.pop(0)
    while len(dividend) >= len(divisor):
        multiplier = dividend[0] / divisor[0]
        for index, coefficient in enumerate(divisor):
            dividend[index] -= multiplier * coefficient
        while dividend and dividend[0] == 0:
            dividend.pop(0)
    return dividend


def certify_squarefree_characteristic_polynomial(base):
    coefficients = exact_characteristic_polynomial(base)
    assert coefficients == [
        1,
        0,
        -66,
        0,
        1627,
        0,
        -18604,
        0,
        98663,
        0,
        -196850,
        0,
        16925,
    ]
    derivative = [
        (len(coefficients) - index - 1) * coefficient
        for index, coefficient in enumerate(coefficients[:-1])
    ]
    first = [Fraction(coefficient) for coefficient in coefficients]
    second = [Fraction(coefficient) for coefficient in derivative]
    while second:
        first, second = second, polynomial_remainder(first, second)
    assert len(first) == 1


def base_data():
    data = json.loads(
        (HERE / "dependent_profile_recovery_witness.json").read_text()
    )
    base = np.asarray(data["matrix"], dtype=np.int64)
    permutation = data["signed_antiautomorphism"]["permutation_image"]
    signs = data["signed_antiautomorphism"]["column_signs"]
    chiral = np.zeros_like(base)
    for source, (target, sign) in enumerate(zip(permutation, signs)):
        chiral[target, source] = sign
    companion = chiral @ base
    assert np.array_equal(chiral @ chiral, -np.eye(12, dtype=np.int64))
    assert np.array_equal(companion, companion.T)
    return base, chiral, companion, permutation


def signed_involutions_four():
    answer = []
    identity = np.eye(4, dtype=np.int64)
    for permutation in permutations(range(4)):
        if any(permutation[permutation[i]] != i for i in range(4)):
            continue
        for signs in product((-1, 1), repeat=4):
            if any(signs[permutation[i]] != signs[i] for i in range(4)):
                continue
            matrix = np.zeros((4, 4), dtype=np.int64)
            for i in range(4):
                matrix[permutation[i], i] = signs[i]
            if not np.array_equal(matrix, matrix.T):
                continue
            if not np.array_equal(matrix @ matrix, identity):
                continue
            if not any(np.array_equal(matrix, old) for old in answer):
                answer.append(matrix)
    return answer


def joint_action_profile(base, companion):
    answer = set()
    for spin in product((-1, 1), repeat=12):
        x = np.asarray(spin, dtype=np.int64)
        answer.add((int(x @ base @ x), int(x @ companion @ x)))
    return answer


def enumerate_templates(base, companion):
    pairs = list(combinations(range(4), 2))
    fibre_spins = [
        np.asarray(spin, dtype=np.int64)
        for spin in product((-1, 1), repeat=4)
    ]
    joint_profile = joint_action_profile(base, companion)
    templates = []
    lower_histogram = Counter()

    for involution in signed_involutions_four():
        # A lower bound of 320 is possible only when K has two supported
        # undirected pairs, one +1 and one -1.  We enumerate that class
        # directly; the remaining support sizes are enumerated afterward
        # only for the total histogram assertion.
        for k_size in (2, 4, 6):
            for k_indices in combinations(range(6), k_size):
                k_support = set(k_indices)
                for positive_positions in combinations(
                    range(k_size), k_size // 2
                ):
                    k_signs = [-1] * k_size
                    for position in positive_positions:
                        k_signs[position] = 1
                    k_matrix = np.zeros((4, 4), dtype=np.int64)
                    for index, sign in zip(k_indices, k_signs):
                        a, b = pairs[index]
                        k_matrix[a, b] = k_matrix[b, a] = sign
                    if not np.array_equal(
                        involution @ k_matrix, k_matrix @ involution
                    ):
                        continue

                    r_support = [(i, i) for i in range(4)] + [
                        pairs[index]
                        for index in range(6)
                        if index not in k_support
                    ]
                    for r_signs in product((-1, 1), repeat=len(r_support)):
                        r_matrix = np.zeros((4, 4), dtype=np.int64)
                        for (a, b), sign in zip(r_support, r_signs):
                            r_matrix[a, b] = r_matrix[b, a] = sign
                        if int(r_matrix.sum()) != 8:
                            continue
                        if not np.array_equal(
                            involution @ r_matrix, r_matrix @ involution
                        ):
                            continue

                        rk_values = {
                            (
                                int(v @ r_matrix @ v),
                                int(v @ k_matrix @ v),
                            )
                            for v in fibre_spins
                        }
                        lower_bound = max(
                            abs(r * q_base + k * q_companion)
                            for r, k in rk_values
                            for q_base, q_companion in joint_profile
                        )

                        for e_signs in product((-1, 1), repeat=k_size):
                            e_matrix = np.zeros((4, 4), dtype=np.int64)
                            for index, sign in zip(k_indices, e_signs):
                                a, b = pairs[index]
                                e_matrix[a, b] = sign
                                e_matrix[b, a] = -sign
                            if not np.array_equal(
                                involution @ e_matrix,
                                -e_matrix @ involution,
                            ):
                                continue
                            lower_histogram[lower_bound] += 1
                            templates.append(
                                (
                                    involution,
                                    r_matrix,
                                    k_matrix,
                                    e_matrix,
                                    frozenset(k_support),
                                    lower_bound,
                                )
                            )
    return templates, lower_histogram


def permutation_canonical_key(template):
    matrices = template[:4]
    keys = []
    identity = np.eye(4, dtype=np.int64)
    for permutation in permutations(range(4)):
        reorder = identity[list(permutation)]
        keys.append(
            b"".join(
                (reorder @ matrix @ reorder.T)
                .astype(np.int8)
                .tobytes()
                for matrix in matrices
            )
        )
    return min(keys)


def macro_pairs(permutation):
    answer = []
    used = set()
    for u in range(12):
        if u in used:
            continue
        v = permutation[u]
        answer.append((u, v))
        used.update((u, v))
    assert len(answer) == 6
    return answer


def diagonal_variables(template, macro_pairing):
    involution, _, _, _, k_support, _ = template
    micro_pairs = list(combinations(range(4), 2))
    variables = []
    for u, v in macro_pairing:
        for index, (a, b) in enumerate(micro_pairs):
            if index in k_support:
                continue
            local = np.zeros((4, 4), dtype=np.int64)
            local[a, b] = local[b, a] = 1
            variable = np.zeros((48, 48), dtype=np.int64)
            variable[4 * u : 4 * u + 4, 4 * u : 4 * u + 4] = local
            variable[4 * v : 4 * v + 4, 4 * v : 4 * v + 4] = (
                -involution @ local @ involution
            )
            variables.append(variable)
    assert len(variables) == 24
    return variables


def base_lift_matrix(template, base, chiral, companion):
    _, r_matrix, k_matrix, e_matrix, _, _ = template
    return (
        np.kron(base, r_matrix)
        + np.kron(companion, k_matrix)
        + np.kron(chiral, e_matrix)
    )


def validate_all_plus_matrix(
    matrix, template, base, chiral, macro_permutation
):
    involution = template[0]
    assert np.array_equal(matrix, matrix.T)
    assert np.all(np.diag(matrix) == 0)
    assert set(matrix[~np.eye(48, dtype=bool)]) == {-1, 1}
    for i, j in combinations(range(12), 2):
        block_sum = int(
            matrix[4 * i : 4 * i + 4, 4 * j : 4 * j + 4].sum()
        )
        assert block_sum == 8 * int(base[i, j])
    inherited = np.kron(chiral, involution)
    assert np.array_equal(inherited @ inherited, -np.eye(48, dtype=np.int64))
    assert np.array_equal(inherited @ matrix, -matrix @ inherited)
    assert all(macro_permutation[macro_permutation[i]] == i for i in range(12))


def two_pattern_separator(matrix):
    """Exact maximum when every macro fibre uses one of two common states."""
    fibre_spins = [
        np.asarray(spin, dtype=np.int64)
        for spin in product((-1, 1), repeat=4)
    ]
    macro_spins = np.asarray(
        list(product((-1, 1), repeat=12)), dtype=np.int64
    )
    best = (-1, None, None, None)

    for first_index in range(16):
        for second_index in range(first_index, 16):
            first = fibre_spins[first_index]
            second = fibre_spins[second_index]
            centre = (first + second) // 2
            direction = (second - first) // 2
            offset = np.tile(centre, 12)
            embedding = np.kron(
                np.eye(12, dtype=np.int64), direction[:, None]
            )
            linear = embedding.T @ matrix @ offset
            quadratic = embedding.T @ matrix @ embedding
            constant = int(offset @ matrix @ offset)
            energies = (
                constant
                + 2 * (macro_spins @ linear)
                + np.einsum(
                    "bi,ij,bj->b",
                    macro_spins,
                    quadratic,
                    macro_spins,
                )
            )
            location = int(np.argmax(np.abs(energies)))
            absolute = int(abs(energies[location]))
            if absolute > best[0]:
                spin = offset + embedding @ macro_spins[location]
                best = (
                    absolute,
                    spin.astype(np.int64),
                    int(energies[location]),
                    (first_index, second_index),
                )
    return best


def minimum_absolute_over_diagonal_signs(base_energy, coefficients):
    """Exact min over theta_i in {+-1} by integer subset-sum DP."""
    reachable = {0}
    for coefficient in coefficients:
        reachable = {
            old + sign * int(coefficient)
            for old in reachable
            for sign in (-1, 1)
        }
    return min(abs(base_energy + correction) for correction in reachable)


def main():
    base, chiral, companion, permutation = base_data()
    certify_full_fibre_algebra(base)
    certify_squarefree_characteristic_polynomial(base)
    templates, lower_histogram = enumerate_templates(base, companion)
    assert lower_histogram == Counter(
        {320: 96, 400: 288, 416: 336, 448: 288}
    )

    survivors = [template for template in templates if template[-1] == TARGET]
    orbit_representatives = {}
    for template in survivors:
        orbit_representatives.setdefault(
            permutation_canonical_key(template), template
        )
    assert len(orbit_representatives) == 4

    pairing = macro_pairs(permutation)
    certificates = []
    for orbit, template in enumerate(orbit_representatives.values()):
        base_lift = base_lift_matrix(template, base, chiral, companion)
        variables = diagonal_variables(template, pairing)
        all_plus = base_lift + sum(variables, np.zeros((48, 48), dtype=np.int64))
        validate_all_plus_matrix(
            all_plus, template, base, chiral, permutation
        )

        separated_value, witness, witness_energy, fibre_pair = (
            two_pattern_separator(all_plus)
        )
        base_energy = int(witness @ base_lift @ witness)
        coefficients = [
            int(witness @ variable @ witness) for variable in variables
        ]
        assert witness_energy == base_energy + sum(coefficients)
        minimum = minimum_absolute_over_diagonal_signs(
            base_energy, coefficients
        )
        assert minimum > TARGET
        certificates.append(
            {
                "orbit": orbit,
                "all_plus_witness_value": separated_value,
                "fibre_pair": fibre_pair,
                "base_energy": base_energy,
                "nonzero_coefficients": sorted(
                    coefficient for coefficient in coefficients if coefficient
                ),
                "minimum_over_all_2^24_diagonal_fillings": minimum,
                "witness": tuple(map(int, witness)),
            }
        )

    print("two cross-fibre blocks generate M_4(R): certified over Q")
    print("characteristic polynomial is square-free: certified over Q")
    print("template lower-bound histogram:", dict(sorted(lower_histogram.items())))
    print("320-survivor permutation orbits:", len(orbit_representatives))
    for certificate in certificates:
        print(certificate)
    print("all 1,008 inherited-chiral Clifford templates are ruled out at 320")


if __name__ == "__main__":
    main()
