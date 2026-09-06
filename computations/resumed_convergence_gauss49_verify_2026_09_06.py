"""Exact finite-field certificate for the two Gauss phase orientations.

No external finite-field package or unproved primitive-polynomial table is used.
Running the complete orbit proves that F2[x]/(x^21+x^2+1) is a field and x is
primitive: all 2^21-1 nonzero elements occur as units in this orbit.
"""
import numpy as np

N, V = 21, 49
MODULUS = (1 << N) | (1 << 2) | 1
ORDER = (1 << N) - 1


def multiply(a, b):
    out = 0
    while b:
        if b & 1:
            out ^= a
        b >>= 1
        a <<= 1
        if a >> N:
            a ^= MODULUS
    return out


def main():
    trace_mask = 0
    for i in range(N):
        a, trace = 1 << i, 0
        for _ in range(N):
            trace ^= a
            a = multiply(a, a)
        assert trace in (0, 1)
        trace_mask |= trace << i
    bins = np.zeros(V, dtype=np.int64)
    a = 1
    for power in range(ORDER):
        if power:
            assert a != 1
        bins[power % V] += 1 - 2*(bin(a & trace_mask).count('1') % 2)
        a <<= 1
        if a >> N:
            a ^= MODULUS
    assert a == 1
    qr = {1, 2, 4}
    for j, count in enumerate(bins):
        if j % 7:
            expected = -81 if j % 7 in qr else 47
        else:
            expected = -337 if j == 0 or j//7 in qr else 687
        assert count == expected, (j, count, expected)
    assert bins.sum() == -1
    # Exact polynomial identity: (-1-i sqrt(7))^7 = 832-448 i sqrt(7).
    real, imag = 1, 0
    for _ in range(7):
        real, imag = -real + 7*imag, -real-imag
    assert (real, imag) == (832, -448)
    # Numerical display only; all certificate assertions above are integral.
    gauss = np.fft.ifft(bins) * V / 2**(N/2)
    print('PASS: complete primitive orbit, trace, and all 49 integer bins')
    print('G(chi)/sqrt(q):', gauss[1])
    print('G(chi^7)/sqrt(q):', gauss[7])
    print('bins:', bins.tolist())


if __name__ == '__main__':
    main()
