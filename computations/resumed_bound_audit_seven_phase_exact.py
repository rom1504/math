"""Exact integer checks of the seven-coset involution and even lift."""

import numpy as np

chi = np.array([0, 1, 1, -1, 1, -1, -1], dtype=np.int64)
P = np.array([[int((j + k) % 7 == 0) for k in range(7)] for j in range(7)], dtype=np.int64)
Cnum = np.array([[chi[(j + k) % 7] for k in range(7)] for j in range(7)], dtype=np.int64)
assert np.array_equal(Cnum @ Cnum, 7 * np.eye(7, dtype=np.int64) - np.ones((7, 7), dtype=np.int64))
assert np.array_equal(P @ Cnum, -(Cnum @ P))
assert np.array_equal(Cnum @ np.ones(7, dtype=np.int64), np.zeros(7, dtype=np.int64))


def field_bins(degree, primitive_polynomial):
    size = 1 << degree

    def mul(a, b):
        ans = 0
        while b:
            if b & 1:
                ans ^= a
            b >>= 1
            a <<= 1
            if a & size:
                a ^= primitive_polynomial
        return ans

    powers, val = [], 1
    for _ in range(size - 1):
        powers.append(val)
        val = mul(val, 2)
    assert val == 1 and len(set(powers)) == size - 1
    bins = np.zeros(7, dtype=np.int64)
    for exponent, val in enumerate(powers):
        tr, cur = 0, val
        for _ in range(degree):
            tr ^= cur
            cur = mul(cur, cur)
        assert tr in (0, 1)
        bins[exponent % 7] += 1 - 2 * tr
    return bins


b3 = field_bins(3, 0b1011)
b6 = field_bins(6, 0b1000011)
assert b3[0] == -1
assert any(np.array_equal(b3[1:], sign * chi[1:]) for sign in (-1, 1))
assert b6[0] == 5
assert any(np.array_equal(b6[1:], -np.ones(6, dtype=np.int64) + sign * 2 * chi[1:]) for sign in (-1, 1))
print("PASS: C^2=I-E, PC=-CP, and exact order-seven F8/F64 Gauss-bin shapes.")
print("F8 bins:", b3.tolist())
print("F64 bins:", b6.tolist())

# Floating-point sanity test only: exact carrier identity with deliberately
# incoherent intermediate-layer signs. The proof is the frequency split.
e, r = 3, 2
v, M, m = 7 ** e, 7 ** (e - 1), 7 ** (e - r)
theta = 2 * np.pi / (7 ** r)
lam = np.zeros(v, dtype=complex)
for j in range(1, v):
    unit, valuation = j, 0
    while unit % 7 == 0:
        unit //= 7
        valuation += 1
    orientation = -1 if valuation == 1 else 1
    lam[j] = np.exp(1j * orientation * chi[unit % 7] * (7 ** valuation) * theta)
lam7 = np.array([0] + [np.exp(1j * chi[j] * theta) for j in range(1, 7)])
f = np.array([1, -1, 1, 1, -1, -1, 1.], dtype=float)
g = np.array([-1, -1, 1, -1, 1, 1, 1.], dtype=float)
h = np.array([1, 1, 1, 1, -1, -1, -1.], dtype=float)
xx, yy = np.zeros(v), np.zeros(v)
for a in range(v):
    b, t = a % M, a // M
    xx[a] = h[b % m] * f[t]
    yy[a] = h[((-b) % M) % m] * g[(t + int(b > 0)) % 7]
Kx = v * np.fft.ifft(lam * np.fft.ifft(xx))
Uf = f.mean() + 7 * np.fft.ifft(lam7 * np.fft.ifft(f))
predicted = np.zeros(v, dtype=complex)
for a in range(v):
    b, t = a % M, a // M
    predicted[a] = h[((-b) % M) % m] * Uf[(t + int(b > 0)) % 7] - h.mean() * f.mean()
assert np.max(np.abs(Kx - predicted)) < 1e-12
assert abs(np.mean(yy * Kx) - (np.mean(g * Uf) - h.mean() ** 2 * g.mean() * f.mean())) < 1e-12
print("Numerical sanity only: cyclic carry and resonant carrier identity agree at v=343.")
