"""Exact arithmetic constants for the general-mask Hermite witness.

The Gaussian inequalities are proved in the linked theorem; this program
checks only the rational substitutions. No floating point or solver.
"""
from fractions import Fraction as F
import json


def verify():
    edge_floor = F(4, 3) * (F(43, 100)**2-F(4, 27))
    assert edge_floor == F(9923, 202500)
    assert edge_floor > F(49, 1000)
    degree = 417
    inverse_sqrt_upper = F(48971, 10**6)
    assert degree*inverse_sqrt_upper**2 > 1
    margin = F(49, 1000)-inverse_sqrt_upper
    assert margin == F(29, 10**6)
    mass_norm_lower = margin/F(2**208)
    assert mass_norm_lower > F(7, 10**68)
    # There are 208 odd degrees from 3 through 417, and sqrt(208)<15.
    assert 208 < 15**2
    individual_coefficient_lower = mass_norm_lower/15
    assert individual_coefficient_lower > F(4, 10**69)
    result = dict(status="exact_rational_constant_certificate",
                  certificate_threshold="43/100",
                  edge_coefficient_floor=str(edge_floor),
                  maximum_odd_degree=degree,
                  nonlinear_projection_norm_lower="7/10^68",
                  one_coefficient_absolute_lower="4/10^69",
                  local_slack_proved=False,
                  original_convergence_proved=False)
    print(json.dumps(result, indent=2))
    return result


if __name__ == "__main__":
    verify()
