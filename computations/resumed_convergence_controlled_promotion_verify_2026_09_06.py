"""Exact integer verification of arbitrary vector phase-to-index promotion."""

import numpy as np


def parity(x):
    return bin(int(x)).count("1") & 1


def characters(q):
    return np.array([[1 - 2 * parity(a & b) for a in range(q)]
                     for b in range(q)], dtype=np.int64)


def verify(r, k, control):
    na, ns = 1 << r, 1 << k
    q, selector_count = na * ns, ns * ns
    chars = characters(q)
    phase = np.array([1 - 2 * parity((index >> r) & control[index & (na - 1)])
                      for index in range(q)], dtype=np.int64)
    kernel = (chars * phase[None, :]) @ chars.T
    selector_chars = characters(selector_count)
    inputs = np.empty((selector_count, q, q), dtype=np.int64)
    expected = np.empty_like(inputs)
    sigma = np.array([index ^ (int(control[index & (na - 1)]) << r)
                      for index in range(q)], dtype=np.int64)
    assert np.array_equal(sigma[sigma], np.arange(q))

    for selector in range(selector_count):
        s, t = selector & (ns - 1), selector >> k
        u, b = s, t
        for x in range(q):
            phi1, phi2 = x & (na - 1), x >> r
            label_in = phi1 | (s << r)
            h_in = 1 - 2 * parity(s & (t ^ phi2))
            inputs[selector, x] = h_in * chars[label_in]
            label_out = phi1 | (b << r)
            h_out = 1 - 2 * parity(b & (u ^ phi2))
            expected[selector, x] = (q * ns) * h_out * chars[label_out, sigma]

    after_base = np.einsum("yx,sxa->sya", kernel, inputs)
    actual = np.einsum("os,sya->oya", selector_chars, after_base)
    assert np.array_equal(actual, expected)


if __name__ == "__main__":
    for r, k in ((2, 1), (2, 2)):
        na, ns = 1 << r, 1 << k
        for table in range(ns ** na):
            control = np.array([(table >> (k * a)) & (ns - 1)
                                for a in range(na)], dtype=np.int64)
            verify(r, k, control)
        print(f"PASS r={r}, k={k}, every one of {ns ** na} control truth tables, "
              "all components, exact phase/index/carrier/selector orientation.")
    rng = np.random.default_rng(20260906)
    for _ in range(17):
        verify(3, 2, rng.integers(0, 4, size=8))
    print("PASS 17 additional nonlinear vector controls r=3,k=2.")
