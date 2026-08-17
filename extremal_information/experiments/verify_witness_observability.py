#!/usr/bin/env python3
"""Finite checks for the witness-observability identities.

The verifier exhausts small deterministic de Bruijn systems.  It checks the
uniform-gap formula directly over every endpoint, the rotation witness for
every pair of equal-length histories, and the finite-gap leakage ceiling.
"""

from __future__ import annotations

from itertools import combinations, product
import json


def words(q: int, length: int):
    return list(product(range(q), repeat=length))


def rotate(state, amount: int):
    amount %= len(state)
    return state[amount:] + state[:amount]


def tight_endpoints(states, history):
    length = len(history)
    if length == 0:
        return set(states)
    return {state for state in states if state[-length:] == history}


def rooted_response(states, tight, gap, repeat, rotation, q):
    def readout(state):
        return rotate(state, rotation)[0] / (q - 1)

    return max(
        (0 if state in tight else -gap) + repeat * readout(state)
        for state in states
    )


def verify():
    filtering_checks = 0
    pair_checks = 0
    leakage_checks = 0
    samples = []

    for q in (2, 3):
        for m in range(1, 6):
            states = words(q, m)
            gap = 2 * m + 3
            repeat = max(1, gap // 2)
            assert repeat < gap

            for length in range(1, m + 1):
                histories = words(q, length)
                endpoint_sets = {
                    history: tight_endpoints(states, history)
                    for history in histories
                }

                for history, tight in endpoint_sets.items():
                    for rotation in range(m):
                        observed_coordinate = rotation
                        if observed_coordinate >= m - length:
                            expected_symbol = history[
                                observed_coordinate - (m - length)
                            ]
                            expected = repeat * expected_symbol / (q - 1)
                            actual = rooted_response(
                                states, tight, gap, repeat, rotation, q
                            )
                            assert abs(actual - expected) < 1e-12
                            filtering_checks += 1

                for left, right in combinations(histories, 2):
                    differing = next(
                        index for index in range(length)
                        if left[index] != right[index]
                    )
                    rotation = m - length + differing
                    left_value = rooted_response(
                        states, endpoint_sets[left], gap, repeat, rotation, q
                    )
                    right_value = rooted_response(
                        states, endpoint_sets[right], gap, repeat, rotation, q
                    )
                    expected_gap = (
                        repeat * abs(left[differing] - right[differing])
                        / (q - 1)
                    )
                    assert abs(abs(left_value - right_value) - expected_gap) < 1e-12
                    assert expected_gap >= repeat / (q - 1)
                    pair_checks += 1

            # At long probe horizons every history has the same leading rate.
            long_repeat = 10 * gap
            rates = []
            for history in words(q, m):
                tight = {history}
                value = rooted_response(
                    states, tight, gap, long_repeat, 0, q
                )
                rates.append(value / long_repeat)
            assert max(rates) - min(rates) <= gap / long_repeat + 1e-12
            leakage_checks += len(rates)

            samples.append(
                {
                    "q": q,
                    "m": m,
                    "phases": q**m,
                    "gap": gap,
                    "repeat": repeat,
                }
            )

    return {
        "filtering_checks": filtering_checks,
        "pair_separation_checks": pair_checks,
        "leakage_rate_checks": leakage_checks,
        "samples": samples,
    }


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, sort_keys=True))
