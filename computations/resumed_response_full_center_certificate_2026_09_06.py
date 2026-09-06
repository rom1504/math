"""Exact rational full-noise center gain for the canonical 21-anchor mask.

Reconstructs the finite canonical response, then uses only two-dimensional
Gaussian formulas and weighted Jensen. No inverse-polynomial distribution,
optimizer, floating point, or numerical quadrature enters the certificate.
"""

import argparse
import json
from collections import Counter
from fractions import Fraction as F
from math import comb, factorial
from pathlib import Path

from fresh_finite_anchor_fixed_point_certificate import CHILDREN
from fresh_limit_hierarchical_fixed_point_certificate import Phi_extended
from fresh_limit_mask_ascent_certificate import phi_large
from fresh_limit_rooted_lower_certificate import DIGITS, I, INV_SQRT_2PI, Phi, phi


def canonical_pair():
    alpha, resolvent, degree = F(361, 500), F(17, 5), 200
    rho = [F(t, 10000) for t in
           [8108, -3110, 1441, 1662, -691, -1088, -758, -874, -444,
            331, 638, 496, 572, 502, 357, 563, 392, 452, 230, 286, 330]]
    assert len(rho) == len(CHILDREN) == len(set(CHILDREN))
    assert all(len(ch) % 2 == 0 and all(j < t for j in ch)
               for t, ch in enumerate(CHILDREN))
    total = sum((x*x for x in rho), F(0))
    innovation = 1-total
    assert innovation > 0
    density = phi(alpha)
    mass = 2*Phi(alpha)-1
    he = [F(1), alpha]
    for d in range(1, degree):
        he.append(alpha*he[-1]-d*he[-2])
    grouped = [F(0) for _ in range(degree+1)]
    for d in range(2, degree+1, 2):
        raw = he[d-1]**2/factorial(d)
        for ell in range(d+1):
            grouped[ell] += raw*comb(d, ell)*total**(d-ell)*innovation**ell
    covariance_base = rho[0]*mass
    removed = F(0)
    for t, children in enumerate(CHILDREN[1:], 1):
        denominator, monomial = 1, F(1)
        for j, count in Counter(children).items():
            denominator *= factorial(count)
            monomial *= rho[j]**count
        raw = he[len(children)-1]*monomial
        removed += raw*raw/denominator
        coefficient = -2*density*raw/I(denominator).sqrt()
        covariance_base += rho[t]*coefficient
    grouped[0] -= removed
    assert all(x >= 0 for x in grouped)
    norm = sum((x/(resolvent+ell)**2 for ell, x in enumerate(grouped)), F(0))
    derivative = sum((ell*x/(resolvent+ell)**2 for ell, x in enumerate(grouped)), F(0))
    weighted = sum((x/(resolvent+ell) for ell, x in enumerate(grouped)), F(0))
    derivative_energy = I(derivative/norm)
    assert derivative_energy.hi < 1
    covariance = covariance_base+I(innovation).sqrt()*2*density*weighted/I(norm).sqrt()
    residual_variance = mass-covariance*covariance
    assert residual_variance.lo > 0 and covariance.lo > 0
    return alpha, density, mass, covariance, residual_variance.sqrt(), derivative_energy


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    alpha, density, mass, covariance, residual_sd, derivative = canonical_pair()
    root_mass = mass.sqrt()
    sign_mean = 2*Phi_extended(covariance*alpha/residual_sd)-1
    outer_tail = 1-Phi_extended(root_mass*alpha/residual_sd)
    coefficient_v = 2*density*sign_mean+4*covariance/root_mass*INV_SQRT_2PI*outer_tail
    coefficient_n = 4*residual_sd/root_mass*INV_SQRT_2PI*outer_tail
    projection_variance = coefficient_v*coefficient_v+coefficient_n*coefficient_n
    nonlinear_variance = 1-mass-projection_variance
    assert nonlinear_variance.lo > 0
    nonlinear_sd = nonlinear_variance.sqrt()
    old_value = covariance*coefficient_v+residual_sd*coefficient_n
    signed_center_mean = old_value/mass
    gaussian_argument = signed_center_mean/nonlinear_sd
    gain = 2*mass*(nonlinear_sd*phi_large(gaussian_argument)
                   -signed_center_mean*(1-Phi_extended(gaussian_argument)))
    assert gain.lo > 0
    full_bound = old_value+gain
    target = F(431, 1000)
    assert full_bound.lo > target
    result = {
        "method": "exact_fraction_outward_intervals",
        "digits": DIGITS,
        "alpha": str(alpha),
        "canonical_derivative_energy": derivative.json(),
        "mask_mass": mass.json(),
        "V_W_covariance": covariance.json(),
        "W_conditional_sd": residual_sd.json(),
        "first_projection_V_coefficient": coefficient_v.json(),
        "first_projection_N_coefficient": coefficient_n.json(),
        "first_projection_variance": projection_variance.json(),
        "full_nonlinear_variance": nonlinear_variance.json(),
        "full_nonlinear_sd": nonlinear_sd.json(),
        "old_value": old_value.json(),
        "signed_center_mean_on_mask": signed_center_mean.json(),
        "gaussian_tail_argument": gaussian_argument.json(),
        "weighted_Jensen_gain": gain.json(),
        "full_center_lower_bound": full_bound.json(),
        "comparison_target": str(target),
        "margin_above_target": (full_bound-target).json(),
        "verified": True,
    }
    rendered = json.dumps(result, indent=2)+"\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")


if __name__ == "__main__":
    main()
