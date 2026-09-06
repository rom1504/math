"""Exact fixed-(alpha,c) upper dual for ALL even unit causal inverse shapes.

The primal candidate certificate is read as rational interval data. The
dual uses rational multipliers, exact strongly-concave pointwise bounds,
and convex chords over certified ranges of the conditional polynomial.
No floating quadrature or root solver enters the proof.
"""

import argparse
import json
from fractions import Fraction as F
from pathlib import Path

from fresh_limit_rooted_lower_certificate import I, INV_SQRT_2PI
from fresh_limit_mask_ascent_certificate import cdf_point, phi_large
from resumed_response_full_center_certificate_2026_09_06 import canonical_pair
from resumed_response_conditional_v_certificate_2026_09_06 import conditional_coefficients, hermite_endpoint


def from_json(value):
    return I(F(value["lower"]), F(value["upper"]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primal", type=Path, default=Path(
        "computations/results/resumed_response_feedforward_inverse_certificate_2026_09_06.json"))
    parser.add_argument("--bins", type=int, default=128)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    primal = json.loads(args.primal.read_text())
    assert primal["verified"]
    alpha = F(primal["outer_mask_alpha"])
    core_alpha, resolvent = F(primal["core_alpha"]), F(primal["core_resolvent"])
    _, density0, mass0, _, _, _ = canonical_pair(core_alpha, resolvent)
    projection, _, _ = conditional_coefficients(core_alpha, resolvent, density0, mass0)
    degree = len(projection) - 1
    lam, gam = from_json(primal["lambda"]), from_json(primal["gamma"])
    tau = from_json(primal["residual_variance"]).sqrt()
    covariance = from_json(primal["output_V_covariance"])
    tail_second = from_json(primal["tail_inverse_second_integral"])
    eta, beta = F(852051, 400000000), F(1259877, 5000000000)
    # eta=.0021301275; beta=.0002519754.
    concavity = 2 * eta - 2 * gam * gam * INV_SQRT_2PI / tau
    assert concavity.lo > 0
    # Mehler: sum_(j<=degree-2) h_j(v)^2 <192 for |v|<=1.
    # Use r=99/100, r^(-(degree-2))<8, (1-r²)^(-1/2)<8, exp(...)<3.
    assert degree <= 200 and alpha < 1
    assert F(100, 99) ** (degree - 2) < 8
    assert F(10000, 199) < 64
    derivative2_energy = sum((k * (k - 1) * projection[k] * projection[k]
                              for k in range(2, degree + 1, 2)), I(0))
    second_derivative_bound = (192 * derivative2_energy).sqrt()

    def pointwise_upper(m):
        y = I(0)
        for _ in range(4):
            argument = (lam * m + gam * y) / tau
            assert 0 < argument.lo <= argument.hi < 8
            y = (gam * (2 * cdf_point(argument) - 1) - beta * m) / (2 * eta)
        y0 = (y.lo + y.hi) / 2
        k = lam * m + gam * y0
        argument = k / tau
        s = 2 * cdf_point(argument) - 1
        objective = 2 * tau * phi_large(argument) + k * s - eta * y0 * y0 - beta * m * y0
        residual = gam * s - 2 * eta * y0 - beta * m
        residual_abs = max(abs(residual.lo), abs(residual.hi))
        return objective.hi + residual_abs**2 / (2 * concavity.lo)

    integral_upper, records = I(0), []
    left_data = hermite_endpoint(F(0), degree)
    for index in range(args.bins):
        left, right = alpha * index / args.bins, alpha * (index + 1) / args.bins
        midpoint, radius = (left + right) / 2, (right - left) / 2
        right_data = hermite_endpoint(right, degree)
        midpoint_data = hermite_endpoint(midpoint, degree)
        m_mid = sum((projection[k] * midpoint_data[2][k]
                     for k in range(0, degree + 1, 2)), I(0))
        derivative_mid = sum((projection[k] * I(k).sqrt() * midpoint_data[2][k - 1]
                              for k in range(2, degree + 1, 2)), I(0))
        deviation = radius * max(abs(derivative_mid.lo), abs(derivative_mid.hi))
        deviation += radius * radius * second_derivative_bound.hi / 2
        lower, upper = m_mid.lo - deviation, m_mid.hi + deviation
        assert 0 < lower < upper
        bin_mass = 2 * (right_data[0] - left_data[0])
        bin_mean = projection[0] * bin_mass
        for k in range(2, degree + 1, 2):
            bin_mean += 2 * projection[k] / I(k).sqrt() * (
                left_data[1] * left_data[2][k - 1]
                - right_data[1] * right_data[2][k - 1])
        weight_left = (upper * bin_mass - bin_mean) / (upper - lower)
        weight_right = (bin_mean - lower * bin_mass) / (upper - lower)
        assert weight_left.lo >= 0 and weight_right.lo >= 0
        f_left, f_right = pointwise_upper(lower), pointwise_upper(upper)
        chord = f_left * weight_left + f_right * weight_right
        integral_upper += chord
        records.append({"left": str(left), "right": str(right),
                        "polynomial_range_lower": str(lower), "polynomial_range_upper": str(upper),
                        "dual_integral_chord": chord.json()})
        left_data = right_data
    dual = eta + beta * covariance + beta * beta * tail_second / (4 * eta) + integral_upper
    lower = from_json(primal["lower_bound"])
    assert dual.hi > lower.lo
    result = {
        "method": "exact_fraction_global_fixed_covariance_causal_shape_dual",
        "verified": True, "scope": "fixed outer alpha and fixed output covariance; all even unit u(V)",
        "outer_alpha": str(alpha), "covariance": covariance.json(),
        "eta": str(eta), "beta": str(beta), "bins": args.bins,
        "pointwise_strong_concavity": concavity.json(),
        "conditional_polynomial_second_derivative_bound": second_derivative_bound.json(),
        "dual_upper_enclosure": dual.json(), "primal_lower_enclosure": lower.json(),
        "certified_dual_gap_upper": I(dual.hi - lower.lo).json(),
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
