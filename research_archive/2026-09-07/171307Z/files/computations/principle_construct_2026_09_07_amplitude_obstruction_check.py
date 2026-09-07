"""Exact finite-field/Fourier checks for the amplitude-damping obstruction."""
from fractions import Fraction
import json
import math


def mul(a, b, q, polynomial):
    z = 0
    while b:
        if b & 1:
            z ^= a
        b >>= 1
        a <<= 1
        if a & q:
            a ^= polynomial
    return z


def trace(a, q, polynomial):
    z, b = 0, a
    for _ in range(q.bit_length()-1):
        z ^= b
        b = mul(b, b, q, polynomial)
    assert z in (0, 1)
    return z


def transform(v):
    v = v[:]
    h = 1
    while h < len(v):
        for start in range(0, len(v), 2*h):
            for j in range(start, start+h):
                a, b = v[j], v[j+h]
                v[j], v[j+h] = a+b, a-b
        h *= 2
    return v


records = []
for q, polynomial in [(4, 7), (8, 11), (16, 19), (32, 37)]:
    n = q*q
    connection = []
    for u in range(q):
        for v in range(q):
            z = mul(u, v, q, polynomial)
            connection.append(int(trace(z, q, polynomial) == 1 and
                                  trace(mul(2, z, q, polynomial), q, polynomial) == 0))
    assert connection[0] == 0
    degree = sum(connection)
    assert degree == (q*q-q)//4
    spectrum = transform(connection)
    assert spectrum[0] == degree
    assert max(map(abs, spectrum[1:])) <= Fraction(3*q, 4)
    theta = Fraction(degree, n-1)
    a2, b2 = theta/(1-theta), (1-theta)/theta
    assert (1-theta)*a2+theta*b2 == 1
    assert b2 == Fraction(3*q+4, q)
    a, b = math.sqrt(a2), math.sqrt(b2)
    op = max(abs(-a-(a+b)*x) for x in spectrum[1:])
    records.append({"q": q, "n": n, "negative_degree": degree,
                    "max_nonconstant_D_eigenvalue_abs": max(map(abs, spectrum[1:])),
                    "a_squared": str(a2), "b_squared": str(b2),
                    "weighted_spectral_cap_upper_normalized_DISPLAY": op/(2*math.sqrt(n)),
                    "sign_cap_all_ones_exact": n*(n-1-2*degree)//2})
print(json.dumps({"status": "EXACT_DEGREES_AND_FOURIER_BOUNDS_PASS", "records": records}, indent=2))
