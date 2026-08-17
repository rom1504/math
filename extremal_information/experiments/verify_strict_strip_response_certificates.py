#!/usr/bin/env python3
"""Exact checks for strict-strip response and quantization certificates.

The script verifies three claims using only integer/Fraction arithmetic:

1. Strict width three has a full-dimensional boundary-response cube.  The
   certificate uses nearest-neighbour vertical edges and same-row horizontal
   edges only.
2. A two-letter width-two column alphabet has seven reachable normalized
   messages but only two weighted residual states.
3. Repeated normalized rounding of an antiferromagnetic Ising transfer can
   keep the projective control state periodic while its scalar reward drifts
   linearly.

No numerical tolerance is used.
"""

from __future__ import annotations

from collections import deque
from fractions import Fraction
from itertools import product
import json


Q = Fraction


def spin_states(width: int) -> list[tuple[int, ...]]:
    return list(product((-1, 1), repeat=width))


def determinant(matrix: list[list[Fraction]]) -> Fraction:
    """Return the determinant by exact Gaussian elimination."""

    work = [row[:] for row in matrix]
    size = len(work)
    sign = 1
    for column in range(size):
        pivot_row = next(
            (row for row in range(column, size) if work[row][column]), None
        )
        if pivot_row is None:
            return Q(0)
        if pivot_row != column:
            work[column], work[pivot_row] = work[pivot_row], work[column]
            sign *= -1
        pivot = work[column][column]
        for row in range(column + 1, size):
            multiplier = work[row][column] / pivot
            for index in range(column, size):
                work[row][index] -= multiplier * work[column][index]
    answer = Q(sign)
    for index in range(size):
        answer *= work[index][index]
    return answer


def inverse(matrix: list[list[Fraction]]) -> list[list[Fraction]]:
    """Return the inverse by exact Gauss-Jordan elimination."""

    size = len(matrix)
    work = [
        row[:]
        + [Q(1 if source == target else 0) for target in range(size)]
        for source, row in enumerate(matrix)
    ]
    for column in range(size):
        pivot_row = next(
            row for row in range(column, size) if work[row][column]
        )
        work[column], work[pivot_row] = work[pivot_row], work[column]
        pivot = work[column][column]
        work[column] = [entry / pivot for entry in work[column]]
        for row in range(size):
            if row == column:
                continue
            multiplier = work[row][column]
            work[row] = [
                left - multiplier * right
                for left, right in zip(work[row], work[column])
            ]
    return [row[size:] for row in work]


def matrix_vector(
    matrix: list[list[Fraction]], vector: list[Fraction]
) -> list[Fraction]:
    return [
        sum((coefficient * value for coefficient, value in zip(row, vector)), Q(0))
        for row in matrix
    ]


WIDTH3_STATES = spin_states(3)
BASE_FIRST_FIELDS = (Q(-3), Q(-7), Q(0))
BASE_FIRST_VERTICALS = (Q(-8), Q(-8))
BASE_HORIZONTALS = (Q(-2), Q(0), Q(5))


def first_column_energy(
    state: tuple[int, int, int],
    fields: tuple[Fraction, Fraction, Fraction] = BASE_FIRST_FIELDS,
) -> Fraction:
    return (
        sum((fields[index] * state[index] for index in range(3)), Q(0))
        + BASE_FIRST_VERTICALS[0] * state[0] * state[1]
        + BASE_FIRST_VERTICALS[1] * state[1] * state[2]
    )


def active_transition_score(
    old: tuple[int, int, int],
    new: tuple[int, int, int],
    fields: tuple[Fraction, Fraction, Fraction] = BASE_FIRST_FIELDS,
    horizontals: tuple[Fraction, Fraction, Fraction] = BASE_HORIZONTALS,
) -> Fraction:
    return first_column_energy(old, fields) + sum(
        (
            horizontals[index] * old[index] * new[index]
            for index in range(3)
        ),
        Q(0),
    )


def base_width3_data() -> tuple[
    list[Fraction], list[tuple[int, int, int]], list[Fraction]
]:
    responses: list[Fraction] = []
    predecessors: list[tuple[int, int, int]] = []
    margins: list[Fraction] = []
    for new in WIDTH3_STATES:
        scored = sorted(
            ((active_transition_score(old, new), old) for old in WIDTH3_STATES),
            reverse=True,
        )
        # A decoupled field of value one in a prior column adds one uniformly.
        responses.append(Q(1) + scored[0][0])
        predecessors.append(scored[0][1])
        margins.append(scored[0][0] - scored[1][0])
    return responses, predecessors, margins


def width3_feature_matrix(
    predecessors: list[tuple[int, int, int]],
) -> list[list[Fraction]]:
    """Columns are 1,y1,y2,y3,y1y2,y2y3,x1,x2*y2."""

    rows: list[list[Fraction]] = []
    for new, old in zip(WIDTH3_STATES, predecessors):
        rows.append(
            [
                Q(1),
                Q(new[0]),
                Q(new[1]),
                Q(new[2]),
                Q(new[0] * new[1]),
                Q(new[1] * new[2]),
                Q(old[0]),
                Q(old[1] * new[1]),
            ]
        )
    return rows


def perturbed_width3_response(
    parameters: list[Fraction],
) -> tuple[list[Fraction], list[set[tuple[int, int, int]]], list[Fraction]]:
    """Evaluate the legal three-column strip for eight parameter changes.

    The parameters correspond to the feature columns in
    ``width3_feature_matrix``.  The constant coordinate is supplied by
    changing a positive decoupled field in the first (dummy) column.
    """

    dummy_field = Q(1) + parameters[0]
    assert dummy_field > 0
    old_fields = (
        BASE_FIRST_FIELDS[0] + parameters[6],
        BASE_FIRST_FIELDS[1],
        BASE_FIRST_FIELDS[2],
    )
    horizontals = (
        BASE_HORIZONTALS[0],
        BASE_HORIZONTALS[1] + parameters[7],
        BASE_HORIZONTALS[2],
    )

    def new_column_energy(new: tuple[int, int, int]) -> Fraction:
        return (
            parameters[1] * new[0]
            + parameters[2] * new[1]
            + parameters[3] * new[2]
            + parameters[4] * new[0] * new[1]
            + parameters[5] * new[1] * new[2]
        )

    responses: list[Fraction] = []
    maximizers: list[set[tuple[int, int, int]]] = []
    margins: list[Fraction] = []
    for new in WIDTH3_STATES:
        scores = {
            old: dummy_field
            + active_transition_score(old, new, old_fields, horizontals)
            + new_column_energy(new)
            for old in WIDTH3_STATES
        }
        optimum = max(scores.values())
        winners = {old for old, value in scores.items() if value == optimum}
        distinct_values = sorted(set(scores.values()), reverse=True)
        margin = (
            distinct_values[0] - distinct_values[1]
            if len(distinct_values) > 1
            else Q(0)
        )
        responses.append(optimum)
        maximizers.append(winners)
        margins.append(margin)
    return responses, maximizers, margins


def verify_width3_cube() -> dict[str, object]:
    base_response, predecessors, base_margins = base_width3_data()
    assert base_response == [Q(v) for v in (18, 28, 18, 28, 20, 24, 20, 24)]
    assert predecessors == [
        (1, -1, 1),
        (1, -1, 1),
        (1, -1, 1),
        (1, -1, 1),
        (-1, 1, -1),
        (1, -1, 1),
        (-1, 1, -1),
        (1, -1, 1),
    ]
    assert base_margins == [Q(v) for v in (2, 14, 2, 14, 6, 6, 6, 6)]

    feature_matrix = width3_feature_matrix(predecessors)
    feature_determinant = determinant(feature_matrix)
    feature_inverse = inverse(feature_matrix)
    inverse_norm = max(sum(map(abs, row), Q(0)) for row in feature_inverse)
    assert feature_determinant == -1024
    assert inverse_norm == 2

    # Find the exact largest centered sup-cube on which each displayed base
    # predecessor remains a maximizer.  For a desired response perturbation d,
    # only parameters 6 and 7 affect predecessor comparisons, and they equal
    # rows 6 and 7 of M^{-1} times d.
    stability_radius: Fraction | None = None
    for new, winner in zip(WIDTH3_STATES, predecessors):
        winner_score = active_transition_score(winner, new)
        for competitor in WIDTH3_STATES:
            if competitor == winner:
                continue
            gap = winner_score - active_transition_score(competitor, new)
            coefficients = [
                (winner[0] - competitor[0]) * feature_inverse[6][index]
                + (winner[1] - competitor[1])
                * new[1]
                * feature_inverse[7][index]
                for index in range(8)
            ]
            sensitivity = sum(map(abs, coefficients), Q(0))
            if sensitivity:
                candidate = gap / sensitivity
                stability_radius = (
                    candidate
                    if stability_radius is None
                    else min(stability_radius, candidate)
                )
    # Positivity of the dummy field permits radius one, so predecessor
    # stability is the active constraint.
    assert stability_radius == Q(1, 2)

    # Linear predecessor-gap inequalities attain their minima at cube
    # corners.  Checking all corners therefore certifies the entire cubes,
    # not merely these 256 samples.
    unique_radius = Q(1, 8)
    minimum_unique_margin: Fraction | None = None
    for signs in product((-1, 1), repeat=8):
        desired = [unique_radius * sign for sign in signs]
        parameters = matrix_vector(feature_inverse, desired)
        response, maximizers, _ = perturbed_width3_response(parameters)
        assert response == [
            base_response[index] + desired[index] for index in range(8)
        ]
        for index, winner in enumerate(predecessors):
            assert maximizers[index] == {winner}

        # Compute the margin specifically against the best competitor; the
        # generic helper's distinct-value margin agrees while maxima are unique.
        for new, winner in zip(WIDTH3_STATES, predecessors):
            old_fields = (
                BASE_FIRST_FIELDS[0] + parameters[6],
                BASE_FIRST_FIELDS[1],
                BASE_FIRST_FIELDS[2],
            )
            horizontals = (
                BASE_HORIZONTALS[0],
                BASE_HORIZONTALS[1] + parameters[7],
                BASE_HORIZONTALS[2],
            )
            winner_score = active_transition_score(
                winner, new, old_fields, horizontals
            )
            competitor_score = max(
                active_transition_score(old, new, old_fields, horizontals)
                for old in WIDTH3_STATES
                if old != winner
            )
            margin = winner_score - competitor_score
            minimum_unique_margin = (
                margin
                if minimum_unique_margin is None
                else min(minimum_unique_margin, margin)
            )
    assert minimum_unique_margin == Q(3, 2)

    # At the maximal radius 1/2 certified by this fixed affine inverse, ties
    # are allowed, but the displayed predecessor is always still a maximizer.
    # Therefore the whole closed response cube is realized exactly.
    closed_radius = stability_radius
    tie_count = 0
    for signs in product((-1, 1), repeat=8):
        desired = [closed_radius * sign for sign in signs]
        parameters = matrix_vector(feature_inverse, desired)
        response, maximizers, _ = perturbed_width3_response(parameters)
        assert response == [
            base_response[index] + desired[index] for index in range(8)
        ]
        for index, winner in enumerate(predecessors):
            assert winner in maximizers[index]
            tie_count += len(maximizers[index]) > 1
    assert tie_count > 0

    return {
        "base_response": [int(value) for value in base_response],
        "base_minimum_margin": int(min(base_margins)),
        "feature_determinant": int(feature_determinant),
        "inverse_infinity_norm": int(inverse_norm),
        "unique_cube_radius": str(unique_radius),
        "unique_cube_minimum_margin": str(minimum_unique_margin),
        "certificate_closed_cube_radius": str(closed_radius),
        "certificate_boundary_tie_incidents": tie_count,
    }


WIDTH2_STATES = spin_states(2)
Column = tuple[Fraction, Fraction, Fraction, Fraction, Fraction]


def width2_column_update(
    message: tuple[Fraction, ...], column: Column
) -> tuple[Fraction, ...]:
    h1, h2, vertical, j1, j2 = map(Q, column)
    updated = []
    for new in WIDTH2_STATES:
        local = h1 * new[0] + h2 * new[1] + vertical * new[0] * new[1]
        updated.append(
            local
            + max(
                message[index]
                + j1 * old[0] * new[0]
                + j2 * old[1] * new[1]
                for index, old in enumerate(WIDTH2_STATES)
            )
        )
    return tuple(updated)


def normalized_width2_update(
    message: tuple[Fraction, ...], column: Column
) -> tuple[tuple[Fraction, ...], Fraction]:
    updated = width2_column_update(message, column)
    reward = max(updated)
    return tuple(value - reward for value in updated), reward


def column_span_bound(column: Column) -> Fraction:
    h1, h2, vertical, j1, j2 = map(Q, column)

    def local(state: tuple[int, int]) -> Fraction:
        return h1 * state[0] + h2 * state[1] + vertical * state[0] * state[1]

    return max(
        local(first)
        - local(second)
        + 2
        * sum(
            abs(coupling)
            for changed, coupling in zip(
                (first[0] != second[0], first[1] != second[1]), (j1, j2)
            )
            if changed
        )
        for first in WIDTH2_STATES
        for second in WIDTH2_STATES
    )


def refine_weighted_partition(
    transitions: list[list[tuple[int, Fraction]]],
) -> list[int]:
    labels = [0] * len(transitions)
    while True:
        signatures = [
            tuple((reward, labels[target]) for target, reward in row)
            for row in transitions
        ]
        signature_labels: dict[tuple[tuple[Fraction, int], ...], int] = {}
        refined = []
        for signature in signatures:
            if signature not in signature_labels:
                signature_labels[signature] = len(signature_labels)
            refined.append(signature_labels[signature])
        if refined == labels:
            return labels
        labels = refined


def verify_seven_to_two_quotient() -> dict[str, object]:
    columns: tuple[Column, ...] = (
        (-1, -1, -1, -1, -1),
        (-1, -1, 0, 1, 0),
    )
    seed = (Q(0), Q(0), Q(0), Q(0))
    states = [seed]
    state_index = {seed: 0}
    queue = deque([seed])
    transition_data: dict[tuple[int, int], tuple[int, Fraction]] = {}

    while queue:
        state = queue.popleft()
        source = state_index[state]
        for letter, column in enumerate(columns):
            target_state, reward = normalized_width2_update(state, column)
            if target_state not in state_index:
                state_index[target_state] = len(states)
                states.append(target_state)
                queue.append(target_state)
            transition_data[(source, letter)] = (
                state_index[target_state],
                reward,
            )

    expected_states = {
        (Q(0), Q(0), Q(0), Q(0)),
        (Q(0), Q(0), Q(0), Q(-4)),
        (Q(0), Q(-2), Q(-2), Q(-4)),
        (Q(-2), Q(0), Q(0), Q(-4)),
        (Q(-2), Q(0), Q(0), Q(-2)),
        (Q(0), Q(-2), Q(-4), Q(-6)),
        (Q(-2), Q(0), Q(0), Q(-6)),
    }
    assert set(states) == expected_states

    transitions = [
        [transition_data[(source, letter)] for letter in range(len(columns))]
        for source in range(len(states))
    ]
    labels = refine_weighted_partition(transitions)
    blocks: dict[int, set[tuple[Fraction, ...]]] = {}
    for state, label in zip(states, labels):
        blocks.setdefault(label, set()).add(state)

    expected_a = {
        (Q(0), Q(0), Q(0), Q(0)),
        (Q(0), Q(0), Q(0), Q(-4)),
        (Q(-2), Q(0), Q(0), Q(-4)),
        (Q(-2), Q(0), Q(0), Q(-2)),
        (Q(-2), Q(0), Q(0), Q(-6)),
    }
    expected_b = expected_states - expected_a
    assert set(map(frozenset, blocks.values())) == {
        frozenset(expected_a),
        frozenset(expected_b),
    }

    label_a = next(label for label, block in blocks.items() if block == expected_a)
    label_b = next(label for label, block in blocks.items() if block == expected_b)
    quotient_rows: dict[int, tuple[tuple[int, Fraction], ...]] = {}
    for label, block in blocks.items():
        signatures = {
            tuple((labels[target], reward) for target, reward in transitions[state_index[state]])
            for state in block
        }
        assert len(signatures) == 1
        quotient_rows[label] = signatures.pop()
    assert quotient_rows[label_a] == ((label_a, Q(3)), (label_b, Q(3)))
    assert quotient_rows[label_b] == ((label_a, Q(1)), (label_b, Q(3)))

    bounds = [column_span_bound(column) for column in columns]
    assert bounds == [Q(8), Q(6)]
    universal_bound = (8 + 1) ** 4 - 8**4
    assert universal_bound == 2465

    return {
        "reachable_normalized_states": len(states),
        "coarsest_weighted_residual_states": len(blocks),
        "column_span_bounds": [int(value) for value in bounds],
        "universal_lattice_state_bound": universal_bound,
        "quotient": {
            "A": {"c0": ["A", 3], "c1": ["B", 3]},
            "B": {"c0": ["A", 1], "c1": ["B", 3]},
        },
    }


def nearest_grid(value: Fraction, mesh: Fraction) -> Fraction:
    """Round to a nearest mesh point; tested examples avoid ties."""

    scaled = value / mesh
    lower_integer = scaled.numerator // scaled.denominator
    lower = Q(lower_integer) * mesh
    upper = Q(lower_integer + 1) * mesh
    lower_distance = abs(value - lower)
    upper_distance = abs(value - upper)
    assert lower_distance != upper_distance
    return lower if lower_distance < upper_distance else upper


def verify_repeated_rounding_drift() -> dict[str, object]:
    checks = 0
    representative: dict[str, object] | None = None
    for coupling in (Q(1), Q(2), Q(5)):
        for field in (Q(1, 20), Q(1, 10), Q(1, 5)):
            mesh = Q(1)
            assert 0 < field < coupling
            assert field < mesh / 4
            column: Column = (field, 0, 0, -coupling, 0)

            exact_message = (Q(0),) * 4
            approximate_shape = (Q(0),) * 4
            approximate_baseline = Q(0)
            for depth in range(1, 17):
                exact_message = width2_column_update(exact_message, column)
                exact_optimum = max(exact_message)
                expected_exact = depth * coupling + (depth % 2) * field
                assert exact_optimum == expected_exact

                raw_shape, reward = normalized_width2_update(
                    approximate_shape, column
                )
                expected_raw_shape = (
                    -2 * field,
                    -2 * field,
                    Q(0),
                    Q(0),
                )
                assert raw_shape == expected_raw_shape
                rounded_shape = tuple(
                    nearest_grid(value, mesh) for value in raw_shape
                )
                assert rounded_shape == (Q(0),) * 4
                approximate_shape = rounded_shape
                approximate_baseline += reward
                assert reward == coupling + field
                assert approximate_baseline == depth * (coupling + field)

                expected_error = (
                    depth * field if depth % 2 == 0 else (depth - 1) * field
                )
                assert approximate_baseline - exact_optimum == expected_error
                checks += 1

            # The exact projective state returns to flat after two steps, as
            # does the rounded state after every step.  Nevertheless the
            # two-step scalar-cycle discrepancy is 2s.
            first_shape, first_reward = normalized_width2_update(
                (Q(0),) * 4, column
            )
            second_shape, second_reward = normalized_width2_update(
                first_shape, column
            )
            assert first_shape == (
                -2 * field,
                -2 * field,
                Q(0),
                Q(0),
            )
            assert second_shape == (Q(0),) * 4
            assert first_reward == coupling + field
            assert second_reward == coupling - field
            exact_cycle_reward = first_reward + second_reward
            rounded_cycle_reward = 2 * (coupling + field)
            assert exact_cycle_reward == 2 * coupling
            assert rounded_cycle_reward - exact_cycle_reward == 2 * field

            assert max(map(abs, first_shape)) == 2 * field < mesh / 2
            projective_rounding_error = (
                max(first_shape) - min(first_shape)
            ) / 2
            assert projective_rounding_error == field

            if coupling == 2 and field == Q(1, 10):
                representative = {
                    "K": int(coupling),
                    "s": str(field),
                    "mesh": int(mesh),
                    "local_sup_rounding_error": str(2 * field),
                    "local_projective_rounding_error": str(field),
                    "two_step_reward_holonomy": str(2 * field),
                    "depth_16_absolute_error": str(16 * field),
                }
    assert representative is not None
    return {"exact_depth_checks": checks, "representative": representative}


def main() -> None:
    print(
        json.dumps(
            {
                "repeated_rounding": verify_repeated_rounding_drift(),
                "seven_to_two_quotient": verify_seven_to_two_quotient(),
                "strict_width3_cube": verify_width3_cube(),
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
