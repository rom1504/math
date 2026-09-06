"""Exact rational checks for both-sign Gibbs-cluster cap extraction."""

from fractions import Fraction as F
from itertools import product
import json

from continued_feedback_ternary_latent_interval_certificate_2026_09_06 import (
    entropy_interval, log_interval, sqrt_interval,
)


def certificate_check():
    p = F(31, 32)
    a = F(91470529542342299, 20460000000000000000)
    delta = F(1, 2**24)
    q2 = (1-2*delta)**2
    h_low, h_high = entropy_interval(delta)
    s_low, s_high = sqrt_interval(p)
    old_low = F(1, 2)-a/(8*s_low)
    old_high = F(1, 2)-a/(8*s_high)
    # This form keeps the correlated occurrences of sqrt(p) together.
    new_low = (F(1, 2)-(a+p*h_high)/(8*s_low))/q2
    new_high = (F(1, 2)-(a+p*h_low)/(8*s_high))/q2
    improvement_left_low = p*h_low
    improvement_right_high = 4*delta*(1-delta)*(4*s_high-a)
    assert improvement_left_low > improvement_right_high
    assert new_high < F(499432211, 10**9) < old_low
    return {
        "p": str(p), "t": "4", "a": str(a), "delta": str(delta),
        "old_interval": [str(old_low), str(old_high)],
        "new_interval": [str(new_low), str(new_high)],
        "new_upper_display_only": float(new_high),
        "strict_improvement_lower": str(old_low-new_high),
        "outward_decimal_cap": "0.499432211",
    }


def finite_checks():
    log2_low, log2_high = log_interval(2)
    checked_matrices = 0
    checked_flip_laws = 0
    smallest_margin = None
    for n in range(2, 5):
        edges = [(i, j) for i in range(n) for j in range(i+1, n)]
        spins = list(product((-1, 1), repeat=n))
        for coefficients in product((-1, 1), repeat=len(edges)):
            energies = {
                x: sum(c*x[i]*x[j] for c, (i, j) in zip(coefficients, edges))
                for x in spins
            }
            x_star = max(spins, key=lambda x: abs(energies[x]))
            extreme_energy = energies[x_star]
            cap = abs(extreme_energy)
            z_both = sum((F(2)**v + F(2)**(-v)
                          for v in energies.values()), F(0))
            log_z_low, _ = log_interval(z_both)
            for delta in (F(1, 8), F(1, 4), F(1, 3)):
                q2 = (1-2*delta)**2
                weights = {}
                for xi in spins:
                    flips = xi.count(-1)
                    weight = delta**flips*(1-delta)**(n-flips)
                    y = tuple(x_star[i]*xi[i] for i in range(n))
                    weights[y] = weight
                assert sum(weights.values(), F(0)) == 1
                expected = sum((weights[y]*energies[y] for y in spins), F(0))
                assert expected == q2*extreme_energy
                # E[-log mu] = n*h(delta), checked by coefficients of
                # the two logarithms, without transcendental arithmetic.
                expected_flips = sum((weights[y]*sum(y[i] != x_star[i]
                                     for i in range(n)) for y in spins), F(0))
                assert expected_flips == n*delta
                _, entropy_high = entropy_interval(delta)
                margin = log_z_low-n*entropy_high-log2_high*q2*cap
                assert margin > 0
                smallest_margin = margin if smallest_margin is None else min(
                    smallest_margin, margin)
                checked_flip_laws += 1
            checked_matrices += 1
    assert checked_matrices == 74
    return {
        "all_hollow_sign_matrices_orders_2_to_4": checked_matrices,
        "rational_flip_laws": checked_flip_laws,
        "inverse_temperature": "log(2)",
        "smallest_certified_both_sign_log_margin": str(smallest_margin),
        "arithmetic": "integer energies, rational partitions, outward rational logs",
    }


def main():
    print(json.dumps({"certificate": certificate_check(),
                      "finite_regression": finite_checks()}, indent=2))


if __name__ == "__main__":
    main()
