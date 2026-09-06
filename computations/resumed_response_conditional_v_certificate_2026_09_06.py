"""Exact conditional-V bin Jensen certificate for the full center response.

The 22-dimensional inverse polynomial is projected exactly onto the single
Gaussian V. Every bin mass and first moment is an exact Hermite endpoint
integral enclosed by outward rational arithmetic. No quadrature is used.
"""

import argparse
import json
from collections import Counter
from fractions import Fraction as F
from math import comb, factorial
from pathlib import Path

from fresh_finite_anchor_fixed_point_certificate import CHILDREN
from fresh_limit_mask_ascent_certificate import cdf_point, phi_large
from fresh_limit_rooted_lower_certificate import DIGITS, I, INV_SQRT_2PI, Phi, phi
from resumed_response_full_center_certificate_2026_09_06 import canonical_pair


def conditional_coefficients(alpha, resolvent, density, mass, degree=200):
    rho = [F(t, 10000) for t in
           [8108, -3110, 1441, 1662, -691, -1088, -758, -874, -444,
            331, 638, 496, 572, 502, 357, 563, 392, 452, 230, 286, 330]]
    total = sum((r*r for r in rho), F(0))
    innovation = 1-total
    he = [F(1), alpha]
    for d in range(1, degree):
        he.append(alpha*he[-1]-d*he[-2])
    first_resolvent = [F(0) for _ in range(degree+1)]
    second_resolvent = [F(0) for _ in range(degree+1)]
    for d in range(0, degree+1, 2):
        for ell in range(d+1):
            probability = F(comb(d, ell))*total**(d-ell)*innovation**ell
            first_resolvent[d] += probability/(resolvent+ell)
            second_resolvent[d] += probability/(resolvent+ell)**2
    old_coefficients = [I(0) for _ in range(degree+1)]
    removed_v_squared = [F(0) for _ in range(degree+1)]
    removed_raw = F(0)
    for t, children in enumerate(CHILDREN):
        d = len(children)
        denominator, monomial = 1, F(1)
        for j, count in Counter(children).items():
            denominator *= factorial(count)
            monomial *= rho[j]**count
        v_t = I(F(factorial(d), denominator)).sqrt()*monomial
        old_coefficients[d] += rho[t]*v_t
        removed_v_squared[d] += F(factorial(d), denominator)*monomial*monomial
        if d:
            removed_raw += he[d-1]**2*monomial*monomial/denominator
    norm_stripped = sum((he[d-1]**2/F(factorial(d))*second_resolvent[d]
                         for d in range(2, degree+1, 2)), F(0))
    norm_stripped -= removed_raw/resolvent**2
    assert norm_stripped > 0
    beta = [I(0) for _ in range(degree+1)]
    beta[0] = mass
    projection = old_coefficients.copy()
    for d in range(2, degree+1, 2):
        # Normalize BEFORE converting to fixed absolute-grid intervals:
        # 1/sqrt(d!) would underflow that grid at high degree.
        normalized_hermite = I(he[d-1]**2/F(factorial(d))).sqrt()
        if he[d-1] < 0:
            normalized_hermite = -normalized_hermite
        beta[d] = -2*density*normalized_hermite
        correction = (-I(innovation).sqrt()*normalized_hermite
                      /I(norm_stripped).sqrt()
                      *(first_resolvent[d]-removed_v_squared[d]/resolvent))
        projection[d] += correction
    assert projection[0].lo <= rho[0] <= projection[0].hi
    covariance_reconstructed = sum((projection[d]*beta[d]
                                    for d in range(0, degree+1, 2)), I(0))
    return projection, covariance_reconstructed, I(norm_stripped)


def hermite_endpoint(x, degree):
    values = [I(1), I(x)]
    for d in range(1, degree):
        values.append(x/I(d+1).sqrt()*values[-1]
                      -I(F(d, d+1)).sqrt()*values[-2])
    return Phi(x), phi(x), values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--alpha", type=F, default=F(3623, 5000))
    parser.add_argument("--resolvent", type=F, default=F(3479, 1000))
    parser.add_argument("--bins", type=int, default=64)
    parser.add_argument("--target", type=F, default=F(8629, 20000))
    args = parser.parse_args()
    assert args.bins >= 1
    alpha, density, mass, covariance, conditional_sd, derivative = canonical_pair(args.alpha, args.resolvent)
    projection, reconstructed, norm = conditional_coefficients(alpha, args.resolvent, density, mass)
    check = reconstructed-covariance
    assert check.lo <= 0 <= check.hi
    assert max(abs(check.lo), abs(check.hi)) < F(1, 10**35)
    root_mass = mass.sqrt()
    sign_mean = 2*cdf_point(covariance*alpha/conditional_sd)-1
    tail = 1-cdf_point(root_mass*alpha/conditional_sd)
    coefficient_v = 2*density*sign_mean+4*covariance/root_mass*INV_SQRT_2PI*tail
    coefficient_n = 4*conditional_sd/root_mass*INV_SQRT_2PI*tail
    nonlinear_variance = 1-mass-coefficient_v*coefficient_v-coefficient_n*coefficient_n
    assert nonlinear_variance.lo > 0
    tau = nonlinear_variance.sqrt()
    lam = coefficient_v-coefficient_n*covariance/conditional_sd
    gamma = coefficient_n/conditional_sd
    old_value = covariance*coefficient_v+conditional_sd*coefficient_n
    assert lam.lo > 0 and gamma.lo > 0
    degree = len(projection)-1
    width = alpha/args.bins
    left_data = hermite_endpoint(F(0), degree)
    total_bound, total_mass, total_mean = I(0), I(0), I(0)
    records = []
    for index in range(args.bins):
        left, right = index*width, (index+1)*width
        right_data = hermite_endpoint(right, degree)
        bin_mass = 2*(right_data[0]-left_data[0])
        assert bin_mass.lo > 0
        bin_g = projection[0]*bin_mass
        for d in range(2, degree+1, 2):
            bin_g += 2*projection[d]/I(d).sqrt()*(left_data[1]*left_data[2][d-1]
                                                 -right_data[1]*right_data[2][d-1])
        bin_k = lam*bin_g+gamma*bin_mass
        assert bin_k.lo > 0
        argument = bin_k/(bin_mass*tau)
        assert 0 < argument.lo <= argument.hi < 8
        bin_bound = bin_k*(2*cdf_point(argument)-1)+2*bin_mass*tau*phi_large(argument)
        total_bound += bin_bound
        total_mass += bin_mass
        total_mean += bin_k
        records.append({"left": str(left), "right": str(right),
                        "symmetric_bin_mass": bin_mass.json(),
                        "integral_K": bin_k.json(),
                        "Jensen_lower_bound": bin_bound.json()})
        left_data = right_data
    for difference in [total_mass-mass, total_mean-old_value]:
        assert difference.lo <= 0 <= difference.hi
        assert max(abs(difference.lo), abs(difference.hi)) < F(1, 10**18)
    assert total_bound.lo > args.target
    result = {
        "method": "exact_fraction_conditional_Hermite_endpoint_bin_Jensen",
        "digits": DIGITS, "alpha": str(alpha), "resolvent": str(args.resolvent),
        "degree": degree, "bins_on_positive_half": args.bins,
        "canonical_derivative_energy": derivative.json(),
        "projection_norm_stripped": norm.json(),
        "conditional_projection_covariance_check": check.json(),
        "conditional_projection_coefficients": [x.json() for x in projection],
        "mask_mass": mass.json(), "V_W_covariance": covariance.json(),
        "K_polynomial_scale": lam.json(), "K_mask_scale": gamma.json(),
        "nonlinear_variance": nonlinear_variance.json(), "old_value": old_value.json(),
        "sum_bin_masses": total_mass.json(), "sum_bin_K_means": total_mean.json(),
        "conditional_V_lower_bound": total_bound.json(),
        "target": str(args.target), "margin_above_target": (total_bound-args.target).json(),
        "bin_certificates": records, "verified": True,
    }
    rendered = json.dumps(result, indent=2)+"\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
        print(json.dumps({key: value for key, value in result.items()
                          if key not in ["bin_certificates", "conditional_projection_coefficients"]}, indent=2))
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
