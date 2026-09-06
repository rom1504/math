"""Exact finite identities and labeled numerical checks for the permanent floor.

Run with .venv/bin/python. No input/output files or random seeds are needed.
Rational permanent/projection identities are exact. Transcendental tests
below are 80-digit sanity checks, not replacements for the analytic proof.
"""

from fractions import Fraction as F
from itertools import combinations_with_replacement, product
from math import factorial
import json
import mpmath as mp

mp.mp.dps = 80


def kernel(a, b):
    # Temperature log(2), so every entry is an EXACT dyadic rational.
    return (F(1, 2 ** ((a - b) ** 2)) + F(1, 2 ** ((a + b) ** 2))) / 2


def permanent(matrix):
    n = len(matrix)
    vals = {0: F(1)}
    for row in matrix:
        nxt = {}
        for mask, value in vals.items():
            for j in range(n):
                if not mask >> j & 1:
                    target = mask | (1 << j)
                    nxt[target] = nxt.get(target, F(0)) + value * row[j]
        vals = nxt
    return vals[(1 << n) - 1]


def orbital(v):
    return permanent([[kernel(a, b) for b in v] for a in v]) / factorial(len(v))


def real(q):
    return mp.mpf(q.numerator) / q.denominator


def gaussian(t, variance=1):
    if variance == 0:
        return mp.mpf(0)
    z = t * variance
    rho = 4 * z / (1 + mp.sqrt(1 + 16 * z * z))
    return -z * (1 - rho) + mp.log1p(-rho * rho) / 4


def main():
    count = projections = 0
    smallest_gap = mp.inf
    t = mp.log(2)
    for m in range(1, 8):
        for v in combinations_with_replacement(range(3), m):
            p = orbital(v)
            assert 0 < p <= 1
            count += 1
            for i in range(m):
                # Squared split projection bound, exact rational arithmetic.
                assert p <= orbital(v[:i] + v[i + 1:]) * orbital((v[i],))
                projections += 1
            energy = F(sum(a * a for a in v), m)
            gap = mp.log(real(p)) / (2 * m) - gaussian(t, real(energy))
            assert gap >= -mp.mpf("1e-70")
            if energy:
                smallest_gap = min(smallest_gap, gap)

    # All 2-by-3 sign bases, without any orthogonality assumption.
    bases = 0
    tilted_gap = mp.inf
    for entries in product((-1, 1), repeat=6):
        rows = (entries[:3], entries[3:])
        total_energy = 0
        total = mp.mpf(0)
        for spin in product((-1, 1), repeat=2):
            h = tuple(sum(spin[i] * rows[i][j] for i in range(2)) for j in range(3))
            total_energy += sum(a * a for a in h)
            # v=h/sqrt(2); use temperature log(2)/2 on integer h.
            vv = [mp.mpf(a) / mp.sqrt(2) for a in h]
            mat = [[(mp.exp(-t*(a-b)**2)+mp.exp(-t*(a+b)**2))/2
                    for b in vv] for a in vv]
            pp = sum(mp.fprod(mat[i][perm[i]] for i in range(3))
                     for perm in __import__("itertools").permutations(range(3))) / 6
            total += mp.exp(t * sum(a*a for a in vv)) * mp.sqrt(pp)
        assert total_energy == 4 * 2 * 3
        gap = mp.log(total / 4) / 3 - t - gaussian(t)
        assert gap >= -mp.mpf("1e-70")
        tilted_gap = min(tilted_gap, gap)
        bases += 1

    floor = mp.sqrt(15) / 8
    tstar = 2 * mp.sqrt(15)
    attained = (mp.log(2) + tstar + gaussian(tstar)) / (2 * tstar)
    assert abs(attained - floor) < mp.mpf("1e-70")
    for i in range(1, 1001):
        p = mp.mpf(i) / 1000
        rho = mp.sqrt(1 - 2 ** (-4*p))
        temp = rho / (2 * (1 - rho*rho))
        value = (p*mp.log(2)+temp+gaussian(temp))/(2*temp*mp.sqrt(p))
        assert abs(value-rho/(2*mp.sqrt(p))) < mp.mpf("1e-70")
        assert value >= floor-mp.mpf("1e-70")

    print(json.dumps({
        "exact": {"dyadic_permanents": count, "split_projection_checks": projections,
                  "nonorthogonal_sign_base_energy_checks": bases},
        "numerical_80_digit_sanity_checks_only": {
            "minimum_nonzero_orbital_log_gap": str(smallest_gap),
            "minimum_tilted_sign_base_log_gap": str(tilted_gap),
            "parameter_values_tested": 1000,
            "universal_certificate_floor": str(floor),
            "attaining_temperature_at_p1": str(tstar)},
        "scope": "Finite arithmetic guards. Uniform/asymptotic assertions are analytic."
    }, indent=2))


if __name__ == "__main__":
    main()
