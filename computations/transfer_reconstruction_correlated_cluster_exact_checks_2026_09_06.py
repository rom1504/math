"""Exact finite checks for the fixed-marginal correlated-cluster law."""

from fractions import Fraction as F
from itertools import product
import json

from continued_feedback_ternary_latent_interval_certificate_2026_09_06 import (
    entropy_interval,
)


def clip(value, radius):
    return max(-radius, min(value, radius))


def check_law(matrix, means, beta, kappa, c_variance_bound):
    n, root_n = 4, F(2)
    u_side, v_side = (0, 1), (2, 3)
    variances = [1-m*m for m in means]
    radius = 1/(8*kappa)
    u_laws = []
    for xu in product((-1, 1), repeat=2):
        probability = F(1)
        for a, j in enumerate(u_side):
            probability *= (1+means[j]*xu[a])/2
        fields = {
            i: sum(matrix[i][j]*(xu[a]-means[j])
                   for a, j in enumerate(u_side))/root_n
            for i in v_side
        }
        u_laws.append((xu, probability, fields))
    clip_means = {i: sum(prob*clip(fields[i], radius)
                        for _, prob, fields in u_laws) for i in v_side}
    joint = {}
    entropy_low = sum(entropy_interval((1-means[j])/2)[0] for j in u_side)
    entropy_high = sum(entropy_interval((1-means[j])/2)[1] for j in u_side)
    field_shift_covariance = F(0)
    vv_shift_quadratic = F(0)
    for xu, probability, fields in u_laws:
        shifts = {i: kappa*variances[i]*(clip(fields[i], radius)-clip_means[i])
                  for i in v_side}
        for i in v_side:
            assert abs(shifts[i]) <= variances[i]/4
            h_low, h_high = entropy_interval((1-means[i]-shifts[i])/2)
            entropy_low += probability*h_low
            entropy_high += probability*h_high
            field_shift_covariance += probability*fields[i]*shifts[i]
        vv_shift_quadratic += probability*sum(
            shifts[i]*matrix[i][j]*shifts[j] for i in v_side for j in v_side)
        for xv in product((-1, 1), repeat=2):
            weight = probability
            for a, i in enumerate(v_side):
                weight *= (1+(means[i]+shifts[i])*xv[a])/2
            joint[xu+xv] = weight
    assert sum(joint.values(), F(0)) == 1
    for i in range(n):
        assert sum(weight*x[i] for x, weight in joint.items()) == means[i]
    expected_energy = sum(weight*sum(matrix[i][j]*x[i]*x[j]
                                    for i in range(n) for j in range(i+1, n))
                          for x, weight in joint.items())
    baseline_energy = sum(matrix[i][j]*means[i]*means[j]
                          for i in range(n) for j in range(i+1, n))
    b = beta/root_n
    assert b*(expected_energy-baseline_energy) == (
        beta*field_shift_covariance+beta/(2*root_n)*vv_shift_quadratic)
    baseline_entropy_high = sum(entropy_interval((1-m)/2)[1] for m in means)
    actual_gain_low = entropy_low-baseline_entropy_high+b*(expected_energy-baseline_energy)
    s_max = sum(variances[j] for j in u_side)/n
    assert 64*kappa*kappa*(3*s_max+F(4, n)) <= F(1, 2)
    variance_weight = sum(variances[i]*variances[j]
                          for i in v_side for j in u_side)/n
    predicted_gain = (beta*kappa/2-kappa*kappa*(1+beta*c_variance_bound/2))*variance_weight
    assert actual_gain_low >= predicted_gain > 0
    return actual_gain_low, predicted_gain


def main():
    n = 4
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    beta, c_bound = F(1, 8), F(3, 2)
    delta = F(1, 8)
    q, v = 1-2*delta, 4*delta*(1-delta)
    assert beta*beta*v <= F(1, 48)
    assert beta*c_bound*v <= 1
    assert n >= 64*beta*beta
    smallest_uniform = None
    smallest_nonuniform = None
    for signs in product((-1, 1), repeat=len(edges)):
        matrix = [[F(0) for _ in range(n)] for _ in range(n)]
        for a, (i, j) in zip(signs, edges):
            matrix[i][j] = matrix[j][i] = F(a)
        gain, _ = check_law(matrix, [q]*n, beta, beta/4, c_bound*v)
        uniform_lower = beta*beta*v*v*F(n*n-1, 128*n)
        assert gain >= uniform_lower
        smallest_uniform = gain if smallest_uniform is None else min(smallest_uniform, gain)
        gain, _ = check_law(matrix, [F(1, 2), F(3, 4), F(-1, 2), F(-3, 4)],
                            beta, F(1, 64), c_bound)
        smallest_nonuniform = gain if smallest_nonuniform is None else min(smallest_nonuniform, gain)
    p = F(31, 32)
    current_delta = F(1, 2**24)
    current_variance = 4*current_delta*(1-current_delta)
    beta_squared = 64/p
    beta_c = 8/p
    assert beta_squared*current_variance < F(1, 48)
    assert beta_c*current_variance < 1
    print(json.dumps({
        "all_order_four_signings": 64,
        "fully_enumerated_joint_laws": 128,
        "marginals_preserved_exactly": True,
        "smallest_homogeneous_gain_lower": str(smallest_uniform),
        "smallest_nonuniform_gain_lower": str(smallest_nonuniform),
        "homogeneous_theorem_credit": str(uniform_lower),
        "current_delta": str(current_delta),
        "current_beta_squared_times_variance": str(beta_squared*current_variance),
        "current_beta_C_times_variance": str(beta_c*current_variance),
        "current_asymptotic_credit_per_spin": str(beta_squared*current_variance**2/128),
        "arithmetic": "rational joint laws and outward rational entropy intervals",
    }, indent=2))


if __name__ == "__main__":
    main()
