#!/usr/bin/env python3
"""Exact finite checks for predictor-induced behavioral recoupling.

The first example makes the factor two and the strict behavioral-margin
condition sharp.  The second checks the pair-graph fixed point on a genuine
two-letter predictor.
"""

from collections import deque
from fractions import Fraction as Q


def behavioral_distance(states, letters, transition, output):
    """Return max future output separation for every ordered state pair."""

    answer = {}
    for source in ((s, t) for s in states for t in states):
        seen = {source}
        queue = deque([source])
        value = Q(0)
        while queue:
            left, right = queue.popleft()
            value = max(value, abs(output[left] - output[right]))
            for letter in letters:
                target = (
                    transition[letter][left],
                    transition[letter][right],
                )
                if target not in seen:
                    seen.add(target)
                    queue.append(target)
        answer[source] = value
    return answer


def sharp_factor_two_check():
    # One physical fixed point has response zero.  Its predictor starts with
    # output -1 and moves to output +1 forever, so epsilon=1 is valid while
    # the encoder/transition mismatch has behavioral size exactly 2 epsilon.
    states = (0, 1)
    letters = (0,)
    transition = {0: {0: 1, 1: 1}}
    output = {0: Q(-1), 1: Q(1)}
    distance = behavioral_distance(states, letters, transition, output)
    epsilon = Q(1)
    assert distance[0, 1] == 2 * epsilon

    physical_response = Q(0)
    encoded = 0
    for depth in range(20):
        assert abs(output[encoded] - physical_response) <= epsilon
        encoded = transition[0][encoded]
    return {
        "epsilon": str(epsilon),
        "behavioral_defect": str(distance[0, 1]),
        "strict_margin_is_necessary": True,
    }


def two_letter_pair_graph_check():
    states = (0, 1, 2, 3)
    letters = ("a", "b")
    transition = {
        "a": {0: 1, 1: 1, 2: 3, 3: 3},
        "b": {0: 2, 1: 0, 2: 2, 3: 1},
    }
    output = {0: Q(0), 1: Q(1), 2: Q(-2), 3: Q(3)}
    distance = behavioral_distance(states, letters, transition, output)

    # Bellman fixed-point identity on the finite synchronized pair graph.
    for left in states:
        for right in states:
            expected = max(
                abs(output[left] - output[right]),
                *(distance[
                    transition[letter][left], transition[letter][right]
                ] for letter in letters),
            )
            assert distance[left, right] == expected

    # Every transition is nonexpansive in the future-behavior metric.
    checks = 0
    for letter in letters:
        for left in states:
            for right in states:
                assert distance[
                    transition[letter][left], transition[letter][right]
                ] <= distance[left, right]
                checks += 1
    return {
        "nonexpansive_pair_checks": checks,
        "distinct_behavioral_distances": sorted(
            {str(value) for value in distance.values()}
        ),
    }


def main():
    import json

    print(
        json.dumps(
            {
                "sharp_factor_two": sharp_factor_two_check(),
                "two_letter_pair_graph": two_letter_pair_graph_check(),
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
