"""Exact rational norm/moment/bin certificate for one causal inverse shape.

The existing rich core is unchanged. A new unit inverse u(V) uses the
center, centered conditional inverse on the center, and its tail part.
Every transcendental or algebraic quantity is enclosed outwards; no
floating quadrature, optimizer output, or numerical tail is used.
"""

import argparse
import json
from fractions import Fraction as F
from pathlib import Path

from fresh_limit_rooted_lower_certificate import I, DIGITS, Phi, phi, INV_SQRT_2PI
from fresh_limit_mask_ascent_certificate import cdf_point, phi_large
from resumed_response_full_center_certificate_2026_09_06 import canonical_pair
from resumed_response_conditional_v_certificate_2026_09_06 import conditional_coefficients, hermite_endpoint


def gaussian_power_moments(alpha, degree, order=80):
    """Central even moments by a decreasing alternating integrated series."""
    assert 0 < alpha < 1
    answer = []
    for exponent in range(0, degree + 1, 2):
        term = alpha**(exponent + 1) / (exponent + 1)
        total = term
        for ell in range(order):
            term *= (-alpha * alpha / 2 / (ell + 1)
                     * F(exponent + 2 * ell + 1, exponent + 2 * ell + 3))
            total += term
        next_term = term * (-alpha * alpha / 2 / (order + 1)
                            * F(exponent + 2 * order + 1, exponent + 2 * order + 3))
        assert term >= 0 and next_term < 0
        answer.append(2 * INV_SQRT_2PI * I(total + next_term, total))
    return answer


def to_power_coefficients(hermite_coefficients):
    degree = len(hermite_coefficients) - 1
    result = [I(0) for _ in range(degree + 1)]
    previous, current = [I(1)], [I(0), I(1)]
    result[0] = hermite_coefficients[0]
    for k in range(1, degree):
        following = [I(0) for _ in range(k + 2)]
        first_factor = 1 / I(k + 1).sqrt()
        second_factor = I(F(k, k + 1)).sqrt()
        for j, value in enumerate(current):
            following[j + 1] += first_factor * value
        for j, value in enumerate(previous):
            following[j] -= second_factor * value
        if (k + 1) % 2 == 0:
            for j, value in enumerate(following):
                result[j] += hermite_coefficients[k + 1] * value
        previous, current = current, following
    return result


def central_second_moment(projection, alpha):
    powers = to_power_coefficients(projection)
    degree = len(projection) - 1
    squared = [I(0) for _ in range(2 * degree + 1)]
    for i in range(0, degree + 1, 2):
        for j in range(i, degree + 1, 2):
            squared[i + j] += (1 if i == j else 2) * powers[i] * powers[j]
    moments = gaussian_power_moments(alpha, 2 * degree)
    norm = sum((squared[2 * j] * moments[j] for j in range(degree + 1)), I(0))
    mean = sum((powers[2 * j] * moments[j] for j in range(degree // 2 + 1)), I(0))
    return norm, mean, moments[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--bins", type=int, default=64)
    parser.add_argument("--target", type=F, default=F(431471, 1000000))
    args = parser.parse_args()
    core_alpha, core_resolvent = F(3623, 5000), F(3479, 1000)
    _, core_density, core_mass, core_covariance, _, derivative = canonical_pair(core_alpha, core_resolvent)
    projection, reconstructed, _ = conditional_coefficients(core_alpha, core_resolvent, core_density, core_mass)
    check = reconstructed - core_covariance
    assert check.lo <= 0 <= check.hi
    assert max(abs(check.lo), abs(check.hi)) < F(1, 10**35)
    assert derivative.hi < 1

    alpha = F(91, 125)
    density, mass = phi(alpha), 2 * Phi(alpha) - 1
    degree = len(projection) - 1
    endpoint = hermite_endpoint(alpha, degree)
    beta = [I(0) for _ in range(degree + 1)]
    beta[0] = mass
    for k in range(2, degree + 1, 2):
        beta[k] = -2 * density * endpoint[2][k - 1] / I(k).sqrt()
    covariance = sum((projection[k] * beta[k] for k in range(0, degree + 1, 2)), I(0))
    full_second = sum((z * z for z in projection), I(0))
    center_second, power_mean, power_mass = central_second_moment(projection, alpha)
    for difference in [power_mean - covariance, power_mass - mass]:
        assert difference.lo <= 0 <= difference.hi
        assert max(abs(difference.lo), abs(difference.hi)) < F(1, 10**25)
    assert center_second.hi - center_second.lo < F(1, 10**25)
    center_variance = center_second - covariance * covariance / mass
    tail_second = full_second - center_second
    assert center_variance.lo > 0 and tail_second.lo > 0
    root_mass = mass.sqrt()
    root_center = center_variance.sqrt()
    root_tail = tail_second.sqrt()
    x1, x2 = F(229, 20000), -F(647, 50000)
    x0 = I(1 - x1 * x1 - x2 * x2).sqrt()
    # The orthonormal basis makes ||u||=1 EXACTLY, independently of rounding.
    c = x0 * covariance / root_mass + x1 * root_center + x2 * root_tail
    d = (1 - c * c).sqrt()
    assert 0 < c.lo < c.hi < 1 and d.lo > 0
    lam = 2 * density * (2 * cdf_point(c * alpha / d) - 1)
    gam = 4 * INV_SQRT_2PI * (1 - cdf_point(alpha / d))
    variance = 1 - mass - lam * lam - 2 * c * lam * gam - gam * gam
    assert variance.lo > 0
    tau = variance.sqrt()
    scale = lam + gam * x1 / root_center
    shift = gam * (x0 / root_mass - x1 * covariance / (mass * root_center))

    left_data = hermite_endpoint(F(0), degree)
    total, total_mass = I(0), I(0)
    records = []
    for index in range(args.bins):
        left, right = alpha * index / args.bins, alpha * (index + 1) / args.bins
        right_data = hermite_endpoint(right, degree)
        bin_mass = 2 * (right_data[0] - left_data[0])
        assert bin_mass.lo > 0
        bin_mean = projection[0] * bin_mass
        for k in range(2, degree + 1, 2):
            bin_mean += 2 * projection[k] / I(k).sqrt() * (
                left_data[1] * left_data[2][k - 1]
                - right_data[1] * right_data[2][k - 1])
        bin_K = scale * bin_mean + shift * bin_mass
        argument = bin_K / (bin_mass * tau)
        assert 0 < argument.lo <= argument.hi < 8
        bound = bin_K * (2 * cdf_point(argument) - 1) + 2 * bin_mass * tau * phi_large(argument)
        total += bound
        total_mass += bin_mass
        records.append({"left": str(left), "right": str(right),
                        "mass": bin_mass.json(), "integral_K": bin_K.json(),
                        "Jensen_lower": bound.json()})
        left_data = right_data
    mass_check = total_mass - mass
    assert mass_check.lo <= 0 <= mass_check.hi
    assert total.lo > args.target
    result = {
        "method": "exact_fraction_causal_unit_inverse_shape_conditional_bin_Jensen",
        "verified": True, "digits": DIGITS,
        "core_alpha": str(core_alpha), "core_resolvent": str(core_resolvent),
        "outer_mask_alpha": str(alpha), "degree": degree, "bins": args.bins,
        "core_derivative_energy": derivative.json(),
        "basis_coordinate_x0": x0.json(), "basis_coordinate_x1": str(x1),
        "basis_coordinate_x2": str(x2),
        "mask_mass": mass.json(), "conditional_inverse_norm_squared": full_second.json(),
        "center_inverse_mean_integral": covariance.json(),
        "center_inverse_second_integral": center_second.json(),
        "center_inverse_variance": center_variance.json(),
        "tail_inverse_second_integral": tail_second.json(),
        "output_V_covariance": c.json(), "output_conditional_sd": d.json(),
        "lambda": lam.json(), "gamma": gam.json(), "residual_variance": variance.json(),
        "conditional_K_scale": scale.json(), "conditional_K_shift": shift.json(),
        "lower_bound": total.json(), "target": str(args.target),
        "margin_above_target": (total - args.target).json(),
        "bin_certificates": records,
    }
    rendered = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(rendered)
        print(json.dumps({key: value for key, value in result.items() if key != "bin_certificates"}, indent=2))
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
