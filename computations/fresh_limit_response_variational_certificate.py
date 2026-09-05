"""Exact interval certificate for the optimized smooth-mask response.

The only infinite expression is a positive Hermite series, bounded on both
sides by Parseval and an explicit geometric tail.  All numerical arithmetic
is rational and rounds outwards; no floating point or quadrature is used.
"""

from fractions import Fraction as F
from math import factorial
import json

from fresh_limit_rooted_lower_certificate import I, phi, Phi, DIGITS


def main():
    rho = F(47, 50)
    alpha = F(81, 100)
    cutoff = 100
    hermite = [F(1), alpha]
    for degree in range(1, 2 * cutoff):
        hermite.append(alpha * hermite[-1] - degree * hermite[-2])
    rational_series = sum((rho ** (4 * j) * hermite[2 * j - 1] ** 2
                           / factorial(2 * j)
                           for j in range(1, cutoff + 1)), F(0))
    density = phi(I(alpha))
    p = 2 * Phi(I(alpha)) - 1
    partial_variance = 4 * density * density * I(rational_series)
    tail = I(rho ** (4 * (cutoff + 1))) * p * (1 - p)
    variance = I(partial_variance.lo, (partial_variance + tail).hi)
    s = variance.sqrt()
    result = 2 * density * (rho * p + I(1 - rho * rho).sqrt() * s)
    target = F(19289, 50000)  # 0.38578 exactly.
    assert result.lo > target
    print(json.dumps({
        "method": "exact_fraction_outward_intervals_positive_Hermite_series",
        "grid_decimal_digits": DIGITS,
        "rho": "47/50", "alpha": "81/100", "Hermite_cutoff": cutoff,
        "p": p.json(), "partial_variance": partial_variance.json(),
        "variance_tail_upper": tail.json(), "s": s.json(),
        "lower_bound": result.json(), "target": "19289/50000",
        "margin_above_target": (result - target).json(),
        "verified": True,
    }, indent=2))


if __name__ == "__main__":
    main()
