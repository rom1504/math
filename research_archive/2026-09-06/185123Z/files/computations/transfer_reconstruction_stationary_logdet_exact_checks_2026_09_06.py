"""Exact bounds and labeled numerical identity checks for raw logdet failure."""

from fractions import Fraction as F
from itertools import product
import json
import math

from continued_feedback_ternary_latent_interval_certificate_2026_09_06 import (
    log_interval,
)


def exact_counterexample_bounds():
    beta, c_bound = F(101, 100), F(101, 100)
    x_low, x_high = F(1, 40), F(3, 100)
    assert x_low/(3*(1-x_low)) < beta-1
    assert x_high/3+x_high*x_high/5 > beta-1
    v_low, v_high = 1-x_high, 1-x_low
    assert beta*c_bound > 1
    weighted_norm_upper = beta*c_bound*v_high
    assert weighted_norm_upper == F(397839, 400000) < 1
    log_argument = 1-(beta*v_low)**2
    assert log_argument == F(4018791, 100000000)
    _, log_upper = log_interval(log_argument)
    gaussian_lower = -log_upper/4
    actual_gap_upper = beta/2
    assert gaussian_lower > F(4, 5)
    assert gaussian_lower-actual_gap_upper > F(1, 4)
    return {
        "beta": str(beta), "C": str(c_bound),
        "limiting_squared_mean_bracket": [str(x_low), str(x_high)],
        "variance_bracket": [str(v_low), str(v_high)],
        "weighted_norm_upper": str(weighted_norm_upper),
        "gaussian_logdet_per_spin_lower": str(gaussian_lower),
        "actual_stationary_chaos_per_spin_upper": str(actual_gap_upper),
        "discrepancy_per_spin_lower": str(gaussian_lower-actual_gap_upper),
    }


def numerical_identity_check():
    # The positive triangle is H_4+I restricted to three coordinates.
    # Set a mean exactly in floating point, then choose b from stationarity.
    n, mean = 3, 0.5
    field = math.atanh(mean)
    b = field/((n-1)*mean)
    variance = 1-mean*mean
    partition = 0.0
    centered_mgf = 0.0
    for x in product((-1, 1), repeat=n):
        energy = sum(x[i]*x[j] for i in range(n) for j in range(i+1, n))
        partition += math.exp(b*energy)
        probability = math.prod((1+mean*x_i)/2 for x_i in x)
        chaos = sum((x[i]-mean)*(x[j]-mean)
                    for i in range(n) for j in range(i+1, n))
        centered_mgf += probability*math.exp(b*chaos)
    entropy = -(1+mean)/2*math.log((1+mean)/2)-(1-mean)/2*math.log((1-mean)/2)
    product_value = n*entropy+b*n*(n-1)*mean*mean/2
    residual = abs(math.log(partition)-product_value-math.log(centered_mgf))
    assert residual < 1e-12
    return {"status": "numerical regression only; proof is algebraic",
            "order": n, "mean": mean, "variance": variance,
            "identity_residual": residual}


def main():
    print(json.dumps({"exact_asymptotic_falsifier": exact_counterexample_bounds(),
                      "stationary_identity_regression": numerical_identity_check()}, indent=2))


if __name__ == "__main__":
    main()
