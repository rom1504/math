"""Exact scoped counterexample/attenuation bounds, plus labeled Gaussian regressions."""

from fractions import Fraction as F
from itertools import product
import json
import math

from scipy.integrate import quad
from scipy.special import ndtr, ndtri

from continued_feedback_ternary_latent_interval_certificate_2026_09_06 import log_interval


def exact_checks():
    log2_low, log2_high = log_interval(2)
    assert log2_high < 1
    assert 2*log2_low > 1
    ratio = F(65536, 64827)
    violation_lower = log_interval(ratio)[0]/2
    assert violation_lower > 0
    delta = F(1, 2**24)
    assert F(5, 2)**18 > F(2**24, 6)
    assert delta < F(1, 27)
    assert F(1369, 36)/(1-delta) < 51
    v = 4*delta*(1-delta)
    beta_c = F(8)/F(31, 32)
    assert beta_c*v < F(1, 2)
    ratio_upper = 64*(51*delta)**2
    assert ratio_upper < F(1, 10**8)
    return {"negative_triangle_actual_partition": "49/4",
            "unattenuated_covariance_witness_determinant": "16/27",
            "strict_violation_lower": str(violation_lower),
            "current_copula_to_clipped_credit_ratio_upper": str(ratio_upper),
            "all_claims_in_this_section": "exact outward rational arithmetic"}


def gaussian_pair_checks():
    def density(z):
        return math.exp(-z*z/2)/math.sqrt(2*math.pi)
    cases = 0
    worst_slack = math.inf
    for first_mean, second_mean in product((-.75, 0., .5, .9), repeat=2):
        first_t = float(ndtri((1-first_mean)/2))
        second_t = float(ndtri((1-second_mean)/2))
        first_a, second_a = 2*density(first_t), 2*density(second_t)
        first_w = 1-first_mean**2-first_a**2
        second_w = 1-second_mean**2-second_a**2
        for correlation in (-.7, -.2, .2, .7):
            joint, error = quad(lambda z: density(z)*ndtr(
                (second_t-correlation*z)/math.sqrt(1-correlation**2)),
                -math.inf, first_t, epsabs=2e-12, epsrel=2e-12)
            covariance = 4*(joint-ndtr(first_t)*ndtr(second_t))
            remainder = abs(covariance-first_a*second_a*correlation)
            bound = math.sqrt(first_w*second_w)*correlation**2
            assert remainder <= bound+1e-10
            worst_slack = min(worst_slack, bound-remainder)
            cases += 1
    delta = 2**-24
    beta = 8/math.sqrt(31/32)
    c_bound = 1/math.sqrt(31/32)
    threshold = float(ndtri(delta))
    a_squared = 4*density(threshold)**2
    variance = 4*delta*(1-delta)
    kappa = beta*c_bound*a_squared
    copula = beta**2*a_squared**2/(4*(1+kappa))
    clipped = beta**2*variance**2/128
    return {"status": "numerical regression only, not the proof",
            "bivariate_threshold_cases": cases,
            "minimum_hermite_tail_slack": worst_slack,
            "current_threshold": threshold, "current_a_squared": a_squared,
            "current_attenuation": a_squared/variance,
            "current_copula_credit_per_spin": copula,
            "current_clipped_credit_per_spin": clipped,
            "current_credit_ratio": copula/clipped}


def main():
    print(json.dumps({"exact": exact_checks(), "diagnostics": gaussian_pair_checks()}, indent=2))


if __name__ == '__main__':
    main()
