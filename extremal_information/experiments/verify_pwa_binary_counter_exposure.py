#!/usr/bin/env python3
"""Exact checks for the one-PWA-generator exponential counter orbit."""

from fractions import Fraction as Q
from itertools import combinations, product


def counter_map(state):
    """Continuous min/max lattice circuit on dual rails (u_0..,v_0..)."""

    half = len(state) // 2
    u, v = state[:half], state[half:]
    out_u = [v[0]]
    out_v = [u[0]]
    for index in range(1, half):
        carry = min(u[:index])
        no_carry = max(v[:index])
        out_u.append(max(min(u[index], no_carry),
                         min(v[index], carry)))
        out_v.append(max(min(v[index], no_carry),
                         min(u[index], carry)))
    return tuple(out_u + out_v)


def dual_rail(integer, bits):
    u = tuple(Q((integer >> index) & 1) for index in range(bits))
    return u + tuple(1 - value for value in u)


def probe(state):
    half = len(state) // 2
    return (state[half - 1] - state[-1]) / 2


def constant_weight_words(length):
    weight = length // 2
    supports = tuple(combinations(range(length), weight))
    # Coordinate zero absent first, present second, as required by the probe.
    supports = tuple(sorted(supports, key=lambda support: (0 in support, support)))
    return tuple(
        tuple(Q(int(index in support)) for index in range(length))
        for support in supports
    )


def local_successor(state, alphabet):
    activations = tuple(
        min(state[index] for index, bit in enumerate(word) if bit)
        for word in alphabet
    )
    return tuple(
        max(activations[index] for index, word in enumerate(alphabet)
            if alphabet[(index + 1) % len(alphabet)][coordinate])
        for coordinate in range(len(state))
    )


def block_counter(state, alphabet):
    length = len(alphabet[0])
    blocks = tuple(
        state[start:start + length]
        for start in range(0, len(state), length)
    )
    maximal = alphabet[-1]
    updated = [local_successor(blocks[0], alphabet)]
    for index in range(1, len(blocks)):
        carry = min(
            min(block[coordinate] for coordinate, bit in enumerate(maximal) if bit)
            for block in blocks[:index]
        )
        no_carry = max(
            max(block[coordinate] for coordinate, bit in enumerate(maximal) if not bit)
            for block in blocks[:index]
        )
        successor = local_successor(blocks[index], alphabet)
        updated.append(tuple(
            max(min(value, no_carry), min(next_value, carry))
            for value, next_value in zip(blocks[index], successor)
        ))
    return tuple(value for block in updated for value in block)


def main():
    transition_checks = 0
    separation_checks = 0
    for bits in range(1, 11):
        period = 1 << bits
        orbit = tuple(dual_rail(value, bits) for value in range(period))
        assert len(set(orbit)) == period
        for value, state in enumerate(orbit):
            assert counter_map(state) == orbit[(value + 1) % period]
            transition_checks += 1

        word = tuple(probe(state) for state in orbit)
        assert set(word) == {Q(-1, 2), Q(1, 2)}
        # The half-negative, half-positive word has primitive period 2^m.
        for shift in range(1, period):
            assert any(word[index] != word[(index + shift) % period]
                       for index in range(period))
            separation_checks += 1

    # Finite rational-grid checks of projective homogeneity, monotonicity's
    # range consequence, and span-one invariance beyond the Boolean orbit.
    structural_checks = 0
    values = (Q(0), Q(1, 2), Q(1))
    for bits in range(1, 4):
        for state in product(values, repeat=2 * bits):
            if max(state) - min(state) > 1:
                continue
            image = counter_map(state)
            assert min(state) <= min(image) <= max(image) <= max(state)
            shifted = tuple(value + Q(3, 7) for value in state)
            assert counter_map(shifted) == tuple(
                value + Q(3, 7) for value in image
            )
            structural_checks += 1

    print(f"exact binary-counter transitions: {transition_checks}")
    print(f"distinct cyclic probe-shift checks: {separation_checks}")
    print(f"rational structural checks: {structural_checks}")
    print("one PWA generator plus one identity probe exposes 2^(r/2) phases")

    block_checks = 0
    for length in (2, 4, 6):
        alphabet = constant_weight_words(length)
        base = len(alphabet)
        for blocks in (1, 2):
            period = base ** blocks
            orbit = []
            for value in range(period):
                digits = []
                remainder = value
                for _ in range(blocks):
                    digits.append(alphabet[remainder % base])
                    remainder //= base
                orbit.append(tuple(entry for digit in digits for entry in digit))
            for value, state in enumerate(orbit):
                assert block_counter(state, alphabet) == orbit[(value + 1) % period]
                block_checks += 1
            probe_word = tuple(
                state[(blocks - 1) * length] - Q(1, 2)
                for state in orbit
            )
            assert probe_word.count(Q(-1, 2)) == period // 2
            assert probe_word.count(Q(1, 2)) == period // 2
            for shift in range(1, period):
                assert any(probe_word[index] != probe_word[(index + shift) % period]
                           for index in range(period))
                block_checks += 1
    print(f"constant-weight block-counter checks: {block_checks}")
    print("growing local blocks give polynomial-presentation 2^(r-o(r)) exposure")


if __name__ == "__main__":
    main()
