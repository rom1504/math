"""Exact interval certificate for the hierarchical Gaussian fixed-point bound.

No floating-point arithmetic is used.  The intermediate response is the
finite Hermite resolvent of the central indicator, through degree 200.
The proof permits this unbounded polynomial: only its final mask is rounded.

Run: .venv/bin/python -B computations/fresh_limit_hierarchical_fixed_point_certificate.py
"""

from fractions import Fraction as F
from math import factorial
import json

from fresh_limit_rooted_lower_certificate import (
    I, DIGITS, PI, INV_SQRT_2PI, phi, Phi,
)


def gaussian_integral_fraction(x, degree=256):
    """Enclose integral_0^x exp(-z*z/2) dz for a nonnegative rational x.

    The integrated alternating Taylor terms need not decrease initially.
    Beyond degree 256 they do decrease for x<=8, so the first omitted
    term bounds the alternating tail.  The finite sum is exact rational.
    """
    assert F(0) <= x <= F(8)
    u = x * x / 2
    term = F(1)
    total = term
    for k in range(1, degree + 1):
        term = term * (-u) / k
        total += term / (2 * k + 1)
    assert u < degree + 2
    rem = x * u ** (degree + 1) / (
        factorial(degree + 1) * (2 * degree + 3))
    value = x * total
    return I(value - rem, value + rem)


def Phi_extended(x):
    x = I.get(x)
    assert 0 <= x.lo <= x.hi <= 8
    lo = I(F(1, 2)) + INV_SQRT_2PI * gaussian_integral_fraction(x.lo)
    hi = I(F(1, 2)) + INV_SQRT_2PI * gaussian_integral_fraction(x.hi)
    # Monotonicity of the Gaussian distribution function gives this hull.
    return I(lo.lo, hi.hi)


def main():
    alpha = F(37, 50)
    resolvent = F(97, 10)
    max_degree = 200
    density = phi(I(alpha))
    mass = 2 * Phi(I(alpha)) - 1

    # He_r(alpha), with the probabilists' normalization.
    hermite = [F(1), alpha]
    for r in range(1, max_degree):
        hermite.append(alpha * hermite[-1] - r * hermite[-2])

    norm_tail = F(0)
    derivative_tail = F(0)
    covariance_tail = F(0)
    for r in range(2, max_degree + 1, 2):
        square = hermite[r - 1] ** 2 / factorial(r)
        norm_tail += square / (resolvent + r) ** 2
        derivative_tail += r * square / (resolvent + r) ** 2
        covariance_tail += square / (resolvent + r)

    factor = 4 * density * density
    norm_squared = mass * mass / (resolvent * resolvent) + factor * norm_tail
    derivative_squared = factor * derivative_tail
    unnormalized_covariance = mass * mass / resolvent + factor * covariance_tail
    derivative_energy = derivative_squared / norm_squared
    mean = mass / resolvent / norm_squared.sqrt()
    covariance = unnormalized_covariance / norm_squared.sqrt()
    residual_variance = mass - covariance * covariance

    assert mean.lo > 0
    assert derivative_energy.hi < 1
    assert residual_variance.lo > 0
    residual_sd = residual_variance.sqrt()
    first_arg = covariance * alpha / residual_sd
    second_arg = alpha * mass.sqrt() / residual_sd
    value = (2 * covariance * density * (2 * Phi_extended(first_arg) - 1)
             + 4 * mass.sqrt() * INV_SQRT_2PI * (1 - Phi_extended(second_arg)))

    target = F(213, 500)  # 0.426 exactly.
    assert value.lo > target
    print(json.dumps({
        "method": "exact_fraction_outward_intervals",
        "grid_decimal_digits": DIGITS,
        "gaussian_tail_series_degree": 256,
        "alpha": "37/50",
        "resolvent_a": "97/10",
        "maximum_even_Hermite_degree": max_degree,
        "intermediate_response": "g_r=C*c_r/(a+r), r even 0..200",
        "normalization": "C=(sum c_r^2/(a+r)^2)^(-1/2)",
        "pi": PI.json(),
        "mask_mass": mass.json(),
        "unnormalized_norm_squared": norm_squared.json(),
        "normalized_mean": mean.json(),
        "normalized_derivative_energy": derivative_energy.json(),
        "strict_stability_margin": (1 - derivative_energy).json(),
        "fixed_point_mask_covariance": covariance.json(),
        "residual_variance": residual_variance.json(),
        "first_cdf_argument": first_arg.json(),
        "second_cdf_argument": second_arg.json(),
        "lower_bound": value.json(),
        "target": "213/500",
        "margin_above_target": (value - target).json(),
        "verified": True,
    }, indent=2))


if __name__ == "__main__":
    main()
