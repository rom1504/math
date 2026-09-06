"""Exact lower witness against the annealed Gaussian two-replica inequality."""

from fractions import Fraction as F
import json

from continued_feedback_ternary_latent_interval_certificate_2026_09_06 import (
    entropy_interval, log_interval,
)


def main():
    # q=1-2x. For every p in [31/32,1], rho^2=1-2^(-4p)<=15/16.
    # Choose the FEASIBLE cross-covariance parameter r=rho, not a numerical
    # optimizer. Its Gaussian information is an upper bound on the minimum.
    x = F(1, 2000)
    p_low = F(31, 32)
    alpha_low = F(1, 16)
    ratio = (alpha_low - 2*x + x*x) / (alpha_low * (1-x)**2)
    assert ratio == F(3936016, 3996001)
    assert (1-x)**2 > 1-alpha_low  # Joint covariance is positive definite.
    entropy_low, _ = entropy_interval(x)
    logarithm_low, _ = log_interval(ratio)
    p_one_lower = entropy_low + logarithm_low / 4
    uniform_lower = p_low * entropy_low + logarithm_low / 4
    assert p_one_lower == F(1038119009001, 2000000000000000)
    assert uniform_lower == F(24619155871531, 64000000000000000)
    assert uniform_lower > 0
    print(json.dumps({
        "overlap_q": str(1-2*x),
        "feasible_gaussian_trial": "r=rho",
        "retention_interval": [str(p_low), "1"],
        "annealed_relation": "rho^2=1-2^(-4p)",
        "logarithm_argument": str(ratio),
        "p_one_rate_lower": str(p_one_lower),
        "uniform_rate_lower": str(uniform_lower),
        "strictly_positive": True,
        "arithmetic": "outward rational intervals; no cubic optimizer",
        "scope": "Gaussian two-replica rate only, not an actual second-moment theorem",
    }, indent=2))


if __name__ == "__main__":
    main()
