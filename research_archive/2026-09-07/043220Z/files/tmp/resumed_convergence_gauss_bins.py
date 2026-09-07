"""Exact integer Gauss-period bins for an explicit binary field."""
import numpy as np

n, v = 21, 49
modulus = (1 << n) | (1 << 2) | 1
order = (1 << n) - 1


def multiply(a, b):
    out = 0
    while b:
        if b & 1:
            out ^= a
        b >>= 1
        a <<= 1
        if a >> n:
            a ^= modulus
    return out


trace_mask = 0
for i in range(n):
    a = 1 << i
    trace = 0
    for _ in range(n):
        trace ^= a
        a = multiply(a, a)
    assert trace in (0, 1), (i, trace)
    trace_mask |= trace << i

bins = np.zeros(v, dtype=np.int64)
a = 1
for power in range(order):
    if power:
        assert a != 1, ('not primitive', power)
    bins[power % v] += 1 - 2*(bin(a & trace_mask).count('1') % 2)
    a <<= 1
    if a >> n:
        a ^= modulus
assert a == 1
assert bins.sum() == -1
gauss = np.fft.ifft(bins) * v / 2**(n/2)
zeta = (-1 + 1j*np.sqrt(7)) / np.sqrt(8)
for j in [1, 2, 3, 7, 14, 21]:
    power = 7 if j % 7 == 0 else 1
    print(j, gauss[j], 'sign + error', abs(gauss[j]-zeta**power),
          'sign - error', abs(gauss[j]-zeta.conjugate()**power))
print('BINS', bins.tolist())
