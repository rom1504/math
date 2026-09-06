"""Exact integer checks for cubic affine-MM phase and index modules."""

import itertools
import numpy as np


def parity(value):
    return bin(int(value)).count("1") & 1


def fwht(array):
    out = np.array(array, dtype=np.int64, copy=True)
    step = 1
    while step < len(out):
        block = out.reshape(-1, 2 * step, *out.shape[1:])
        left, right = block[:, :step].copy(), block[:, step:].copy()
        block[:, :step], block[:, step:] = left + right, left - right
        step *= 2
    return out


def characters(q):
    return np.array([[1 - 2 * parity(a & b) for a in range(q)]
                     for b in range(q)], dtype=np.int64)


def polynomial_and_field(r, chosen):
    q = 1 << r
    polynomial = np.zeros(q, dtype=np.int64)
    field = np.zeros(q, dtype=np.int64)
    for term in chosen:
        support = [i for i in range(r) if (term >> i) & 1]
        for v in range(q):
            polynomial[v] ^= int((v & term) == term)
        if len(support) == 1:
            contributions = [(support[0], term)]
        elif len(support) == 2:
            contributions = [(support[0], 1 << support[1])]
        else:
            contributions = [(i, term ^ (1 << i)) for i in support]
        for component, mask in contributions:
            for v in range(q):
                if (v & mask) == mask:
                    field[v] ^= 1 << component
    for a in range(q):
        for v in range(q):
            assert parity(a & (field[v ^ a] ^ field[v])) == polynomial[a]
    return polynomial, field


def verify_phase(r, polynomial, field):
    q = 1 << r
    chars = characters(q)
    labels, carrier, swapped = [], [], []
    for z in range(q * q):
        x, y = z & (q - 1), z >> r
        labels.append(x ^ int(field[y]))
        carrier.append(1 - 2 * parity(x & y))
        swapped.append((x << r) | y)
    labels = np.array(labels)
    carrier = np.array(carrier)
    transformed = fwht(carrier[:, None] * chars[labels])
    expected = q * carrier[:, None] * chars[labels[swapped]] * (1 - 2 * polynomial)[None, :]
    assert np.array_equal(transformed, expected)
    assert np.array_equal(np.bincount(labels, minlength=q), np.full(q, q))


def verify_controlled(r, polynomial, field):
    state_bits = r + 2
    state_count = 1 << state_bits
    label_count = 1 << (r + 1)
    chars = characters(label_count)
    sigma = np.array([a ^ (int(polynomial[a & ((1 << r) - 1)]) << r)
                      for a in range(label_count)])
    assert np.array_equal(sigma[sigma], np.arange(label_count))
    permutations = np.empty((label_count, state_count), dtype=np.int64)
    for index in range(label_count):
        a, c = index & ((1 << r) - 1), index >> r
        for y in range(state_count):
            u, v, w = y & 1, (y >> 1) & ((1 << r) - 1), y >> (r + 1)
            vp = v ^ (u * a)
            wp = w ^ (c * u) ^ parity(a & int(field[v]))
            permutations[index, y] = u | (vp << 1) | (wp << (r + 1))
        assert sorted(permutations[index].tolist()) == list(range(state_count))
    for index in range(label_count):
        assert np.array_equal(permutations[sigma[index], permutations[index]], np.arange(state_count))

    labels, carrier, swapped = [], [], []
    for z in range(state_count * state_count):
        x, y = z & (state_count - 1), z >> state_bits
        u, v = y & 1, (y >> 1) & ((1 << r) - 1)
        xv, xw = (x >> 1) & ((1 << r) - 1), x >> (r + 1)
        labels.append((u * xv) ^ (xw * int(field[v])) ^ ((xw * u) << r))
        carrier.append(1 - 2 * parity(x & y))
        swapped.append((x << state_bits) | y)
    labels, carrier = np.array(labels), np.array(carrier)
    transformed = fwht(carrier[:, None] * chars[labels])
    expected = state_count * carrier[:, None] * chars[labels[swapped]][:, sigma]
    assert np.array_equal(transformed, expected)
    label_characters = fwht(np.bincount(labels, minlength=label_count))
    assert max(abs(label_characters[1:])) <= state_count * state_count // 2


if __name__ == "__main__":
    checked = 0
    r = 3
    terms = list(range(1, 1 << r))
    for chosen_bits in range(1 << len(terms)):
        chosen = [term for i, term in enumerate(terms) if (chosen_bits >> i) & 1]
        polynomial, field = polynomial_and_field(r, chosen)
        verify_phase(r, polynomial, field)
        verify_controlled(r, polynomial, field)
        checked += 1
    print(f"PASS all {checked} zero-constant polynomials on 3 bits: "
          "64-point exact phase modules and 1024-point controlled-index modules.")

    r = 4
    terms = [sum(1 << i for i in support) for degree in (1, 2, 3)
             for support in itertools.combinations(range(r), degree)]
    selections = [terms, terms[::2], terms[1::2]] + [[term] for term in terms]
    for chosen in selections:
        polynomial, field = polynomial_and_field(r, chosen)
        verify_phase(r, polynomial, field)
        verify_controlled(r, polynomial, field)
    print(f"PASS {len(selections)} 4-bit polynomials, including overlapping cubics: "
          "256-point phase modules and 4096-point controlled-index modules.")
