"""Exact checks for the independent 2609.12427v1 audit.

No inputs, downloads, solver calls, or filesystem writes.  Rational checks
certify the finite part of an all-degree bound; the audit supplies the
analytic tail argument.  This is not a numerical search certificate.
"""

from collections import Counter
from fractions import Fraction
from math import ceil, comb, exp, factorial, floor, log, sqrt


def balance_probability(r):
    return Fraction(
        sum(comb(r, d) for d in range(ceil(r / 3), floor(2 * r / 3) + 1)),
        2**r,
    )


def check(name, condition):
    assert condition, name
    print(f"PASS: {name}")


print("Finite balance probabilities for the analytic r >= 32 tail:")
for r in range(8, 32):
    probability = balance_probability(r)
    check(f"p_{r} = {probability} >= 21/32", probability >= Fraction(21, 32))

# Hoeffding gives p_r >= 1 - 2 exp(-r/18).  At r >= 32 this
# exceeds 21/32 because exp(16/9) > 64/11, certified below by
# a positive-term truncation of the exponential series.
taylor = sum((Fraction(16, 9) ** j) / factorial(j) for j in range(6))
check("sum_{j=0}^5 (16/9)^j/j! > 64/11", taylor > Fraction(64, 11))
print(f"Taylor lower bound = {taylor}")

# With the paper's exact damping bound 3, exact probabilities,
# actual block counts, and the actual endpoint theta, these are
# c_5^5, c_6^12, and c_7^14 for weight exponent s=5/2.
small_powers = {
    5: (5, Fraction(2**17 * 3**15, 5**18)),
    6: (12, Fraction(2**67, 3**23 * 5**14)),
    7: (14, Fraction(2**110 * 3**29, 5**8 * 7**50)),
}
for r, (power, value) in small_powers.items():
    check(f"c_{r} < 9/10", value < Fraction(9, 10) ** power)
    print(f"c_{r} upper bound = {float(value) ** (1 / power):.15f}")

# c_4^8 = 27/32.  The coarse r >= 8 factor has
# c_tail^48 = 2^335 / (3^165 7^27).  Compare exact integers.
check("9/10 < c_4 < 1", Fraction(9, 10) ** 8 < Fraction(27, 32) < 1)
check("c_tail < c_4", 2**365 < 3**183 * 7**27)
c4 = (27 / 32) ** (1 / 8)
tail = (32 / 21) ** (9 / 16) * 9 ** (1 / 16) * 2 ** (3 / 2)
tail *= (4 / 27) ** (5 / 6) * sqrt(4 / 3)
print(f"c_4 = {c4:.15f}")
print(f"c_tail = {tail:.15f}")
print(f"1/(1-c_4) = {1/(1-c4):.15f}")

original = (1 - 2 * exp(-8 / 3)) ** (-13 / 24)
original *= 13 ** (1 / 24) * 3 * sqrt(2) * (2 / 3) ** 5
print(f"Paper's displayed contraction c = {original:.15f}")

# A concrete six-variable witness obstructs competitive unrestricted
# degree-two BH constants.  No claim of global optimality is needed.
matrix = [
    [0, 1, 1, 1, 1, 1],
    [1, 0, 1, 1, -1, -1],
    [1, 1, 0, -1, -1, 1],
    [1, 1, -1, 0, 1, -1],
    [1, -1, -1, 1, 0, 1],
    [1, -1, 1, -1, 1, 0],
]
energies = Counter()
for mask in range(64):
    spins = [1 if (mask >> i) & 1 else -1 for i in range(6)]
    energy = sum(matrix[i][j] * spins[i] * spins[j] for i in range(6) for j in range(i + 1, 6))
    energies[energy] += 1
check("H6 energy distribution", energies == Counter({-5: 12, -3: 20, 3: 20, 5: 12}))
print(f"H6 histogram = {dict(sorted(energies.items()))}")
print(f"B_2 >= {15**0.75/5:.15f}")
print(f"Direct asymptotic BH-constant ceiling = {5/30**0.75:.15f}")
