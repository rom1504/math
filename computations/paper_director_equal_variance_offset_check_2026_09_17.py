"""Floating-point diagnostic of the proved equal-variance offset theorem.

X is zero with probability 1-1/K and +/-sqrt(K) with probability 1/(2K).
Its variance is exactly one and Hoeffding gives subGaussian proxy K.
Y is standard Gaussian. This calculation is NOT a theorem certificate.
"""

import argparse
import json
import math
from pathlib import Path

from scipy.integrate import quad
from scipy.optimize import minimize_scalar


def run():
    kappa = math.sqrt(2 / math.pi)
    results = []
    for proxy in (2.0, 4.0, 9.0, 13.0):
        amplitude = math.sqrt(proxy)

        def response_x(offset):
            return ((1 - 1 / proxy) * abs(offset)
                    + (abs(offset + amplitude) + abs(offset - amplitude))
                    / (2 * proxy))

        def response_y(offset):
            return (kappa * math.exp(-offset * offset / 2)
                    + offset * math.erf(offset / math.sqrt(2)))

        def difference(offset):
            return response_x(offset) - response_y(offset)

        discount = kappa - 1 / amplitude
        radius = math.sqrt(2 * proxy * math.log(32 * proxy / discount**2))
        bound = discount**2 / (8 * radius)
        optimum = minimize_scalar(lambda s: -difference(s), bounds=(0, amplitude),
                                  method="bounded", options={"xatol": 1e-13})
        positive = difference(optimum.x)
        # Integrate the separately nonnegative stop-loss functions, avoiding
        # cancellation between them. The tail beyond 16 is negligible here.
        integral_x = 2 * quad(lambda s: response_x(s) - s, 0, amplitude,
                              epsabs=1e-11)[0]
        integral_y = 2 * quad(lambda s: response_y(s) - s, 0, 16,
                              epsabs=1e-11)[0]
        assert positive >= bound
        assert abs(integral_x - 1) < 1e-9
        assert abs(integral_y - 1) < 1e-9
        results.append(dict(proxy=proxy, discount=discount, radius=radius,
                            theorem_loss_lower=bound, diagnostic_max_loss=positive,
                            diagnostic_offset=optimum.x,
                            stop_loss_integral_x=integral_x,
                            stop_loss_integral_gaussian=integral_y))
    positive_examples = []
    cold, hot = 0.25, 4.0
    cold_weight = (hot - 1) / (hot - cold)
    variance_of_variance = (cold_weight * (cold - 1)**2
                           + (1 - cold_weight) * (hot - 1)**2)
    for offset_radius in (0.0, 0.5, 2.0, 5.0):
        def averaged_absolute(variance):
            if offset_radius == 0:
                return kappa * math.sqrt(variance)

            def integrand(s):
                return (kappa * math.sqrt(variance) * math.exp(-s*s / (2*variance))
                        + s * math.erf(s / math.sqrt(2*variance)))

            return quad(integrand, 0, offset_radius, epsabs=1e-12)[0] / offset_radius

        gain = (averaged_absolute(1)
                - cold_weight * averaged_absolute(cold)
                - (1 - cold_weight) * averaged_absolute(hot))
        lower = (math.exp(-offset_radius**2 / (2*cold)) * variance_of_variance
                 / (4 * math.sqrt(2*math.pi) * hot**1.5))
        assert gain + 1e-12 >= lower
        positive_examples.append(dict(offset_radius=offset_radius,
                                      actual_averaged_gain=gain,
                                      curvature_lower=lower))
    return dict(status="floating-point diagnostic PASS; analytic proof is separate",
                examples=results, uniform_offset_examples=positive_examples)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    encoded = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")
