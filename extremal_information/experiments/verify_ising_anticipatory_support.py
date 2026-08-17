#!/usr/bin/env python3
"""Exact checks for the strict width-two Ising anticipatory quotient."""

from fractions import Fraction as Q
from itertools import permutations, product


LETTERS = ("a", "b", "c")
STATES = tuple(product((1, -1), repeat=2))


def baseline(letter, source, target, penalty):
    s1, s2 = source
    t1, _ = target
    if letter == "a":
        return penalty * Q(s2 - 1, 2) + penalty * Q(t1 - 1, 2)
    if letter == "b":
        return penalty * Q(s1 * s2 - 1, 2) + penalty * Q(t1 - 1, 2)
    return -penalty * Q(1 + s1 * s2, 2) - penalty * Q(1 + t1, 2)


def matrix(letter, penalty, interacting=False):
    return tuple(
        tuple(
            baseline(letter, source, target, penalty)
            - (source[0] * target[0] if interacting and letter == "a" else 0)
            for target in STATES
        )
        for source in STATES
    )


def maxplus_product(left, right):
    size = len(left)
    return tuple(
        tuple(
            max(left[i][k] + right[k][j] for k in range(size))
            for j in range(size)
        )
        for i in range(size)
    )


def word_product(matrices, word):
    answer = matrices[word[0]]
    for letter in word[1:]:
        answer = maxplus_product(answer, matrices[letter])
    return answer


def max_cycle_mean(matrix_value):
    size = len(matrix_value)
    best = None
    for length in range(1, size + 1):
        for cycle in permutations(range(size), length):
            total = sum(
                matrix_value[cycle[index]][cycle[(index + 1) % length]]
                for index in range(length)
            )
            mean = total / length
            best = mean if best is None else max(best, mean)
    return best


def tau(letter):
    return 1 if letter in ("a", "b") else -1


def carrier_response(word):
    previous = tau(word[-1])
    answer = 0
    for letter in word:
        if letter == "a":
            answer -= previous
        previous = tau(letter)
    return Q(answer)


def main():
    penalty = Q(5)
    plain = {letter: matrix(letter, penalty) for letter in LETTERS}
    weighted = {
        letter: matrix(letter, penalty, interacting=True) for letter in LETTERS
    }

    supports = {
        sign: {index for index, state in enumerate(STATES) if state[0] == sign}
        for sign in (1, -1)
    }
    support_checks = 0
    for sign in (1, -1):
        for letter in LETTERS:
            reached = {
                target
                for source in supports[sign]
                for target in range(len(STATES))
                if plain[letter][source][target] == 0
            }
            assert reached == supports[tau(letter)]
            support_checks += 1

    # No single nonempty raw support is stable under every zero relation.
    stable_single_supports = 0
    for mask in range(1, 1 << len(STATES)):
        support = {index for index in range(len(STATES)) if (mask >> index) & 1}
        if all(
            support
            <= {
                target
                for source in support
                for target in range(len(STATES))
                if plain[letter][source][target] == 0
            }
            for letter in LETTERS
        ):
            stable_single_supports += 1
    assert stable_single_supports == 0

    signatures = tuple(
        tuple(int(any(plain[letter][source][target] == 0
                      for target in range(len(STATES)))) for letter in LETTERS)
        for source in range(len(STATES))
    )
    assert len(set(signatures)) == 4

    word_checks = 0
    for depth in range(1, 7):
        for word in product(LETTERS, repeat=depth):
            assert max_cycle_mean(word_product(plain, word)) == 0
            expected = carrier_response(word)
            assert max_cycle_mean(word_product(weighted, word)) == expected
            cyclic_ca = sum(
                word[(index - 1) % depth] == "c" and letter == "a"
                for index, letter in enumerate(word)
            )
            assert expected == 2 * cyclic_ca - word.count("a")
            word_checks += 1

    left, right = tuple("aabccb"), tuple("abbcac")
    assert sorted(left) == sorted(right)
    assert carrier_response(left) == -2
    assert carrier_response(right) == 2

    print(f"anticipatory support checks: {support_checks}")
    print(f"plain/weighted word-response checks: {word_checks}")
    print("width-two Ising: scalar output, 2 support states, 4 path states")


if __name__ == "__main__":
    main()
