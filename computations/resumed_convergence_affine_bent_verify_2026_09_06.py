"""Exact integer replay of unitriangular affine bent profile modules."""

import numpy as np


def fwht(a):
    out = np.array(a, dtype=np.int64, copy=True)
    width = 1
    while width < len(out):
        block = out.reshape(-1, 2 * width, *out.shape[1:])
        left = block[:, :width].copy()
        right = block[:, width:].copy()
        block[:, :width] = left + right
        block[:, width:] = left - right
        width *= 2
    return out


def bits(value, width):
    return np.array([(value >> i) & 1 for i in range(width)], dtype=np.int64)


def labels_and_carrier(d, edges):
    count = 1 << (2 * d)
    label = np.empty(count, dtype=np.int64)
    carrier = np.empty(count, dtype=np.int64)
    swapped = np.empty(count, dtype=np.int64)
    for z in range(count):
        x, y = bits(z & ((1 << d) - 1), d), bits(z >> d, d)
        label[z] = sum(int(x[i] * y[j]) << t for t, (i, j) in enumerate(edges))
        carrier[z] = 1 - 2 * (int(x @ y) & 1)
        swapped[z] = ((z & ((1 << d) - 1)) << d) | (z >> d)
    return label, carrier, swapped


def inversion_permutation(d, edges):
    out = []
    ident = np.eye(d, dtype=np.int64)
    for a in range(1 << len(edges)):
        nil = np.zeros((d, d), dtype=np.int64)
        for t, (i, j) in enumerate(edges):
            nil[i, j] = (a >> t) & 1
        inv, power = ident.copy(), ident.copy()
        for _ in range(1, d):
            power = (power @ nil) & 1
            inv ^= power
        assert np.array_equal(((ident ^ nil) @ inv) & 1, ident)
        out.append(sum(int(inv[i, j]) << t for t, (i, j) in enumerate(edges)))
    out = np.array(out, dtype=np.int64)
    assert sorted(out.tolist()) == list(range(len(out)))
    assert np.array_equal(out[out], np.arange(len(out)))
    return out


def verify(d):
    edges = [(i, j) for i in range(d) for j in range(i + 1, d)]
    r, physical_count = len(edges), 1 << (2 * d)
    q = 1 << r
    sigma = inversion_permutation(d, edges)
    label, carrier, swapped = labels_and_carrier(d, edges)
    chars = np.array([[1 - 2 * (bin(a & b).count("1") & 1) for a in range(q)]
                      for b in range(q)], dtype=np.int64)
    components = carrier[:, None] * chars[label]
    transformed = fwht(components)
    expected = (1 << d) * carrier[:, None] * chars[label[swapped]][:, sigma]
    assert np.array_equal(transformed, expected)

    kernel = chars[:, sigma] @ chars.T
    assert np.array_equal(kernel @ kernel.T, q * q * np.eye(q, dtype=np.int64))
    if q <= 8:
        profiles = np.array([[1 - 2 * ((f >> t) & 1) for f in range(1 << q)]
                             for t in range(q)], dtype=np.int64)
    else:
        rng = np.random.default_rng(20260906)
        profiles = np.column_stack((np.ones(q, dtype=np.int64),
                                    -np.ones(q, dtype=np.int64),
                                    1 - 2 * rng.integers(0, 2, size=(q, 64))))
    physical_inputs = carrier[:, None] * profiles[label]
    actual = q * fwht(physical_inputs)
    target_profiles = kernel @ profiles
    target = (1 << d) * carrier[:, None] * target_profiles[label[swapped]]
    assert np.array_equal(actual, target)

    counts = np.bincount(label, minlength=q)
    character_sums = fwht(counts)
    assert character_sums[0] == physical_count
    assert np.max(np.abs(character_sums[1:])) <= physical_count // 2
    print(f"PASS d={d}, labels={q}, physical={physical_count}, "
          f"components={q}, profiles={profiles.shape[1]}, "
          f"max_nontrivial_label_character={np.max(np.abs(character_sums[1:]))}/{physical_count}")


if __name__ == "__main__":
    for dimension in (2, 3, 4):
        verify(dimension)
