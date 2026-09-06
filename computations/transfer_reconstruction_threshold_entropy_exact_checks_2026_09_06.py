"""Exact score normalization and numerical rare-pair entropy diagnostics."""

from fractions import Fraction as F
from itertools import product
import json
import math

from scipy.integrate import quad
from scipy.special import ndtr


def exact_score_checks():
    n = 4
    mean, a = F(3, 4), F(1, 4)
    variance = 1-mean*mean
    attenuation = a*a/variance
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    checked = 0
    for signs in product((-1, 1), repeat=len(edges)):
        expected_score, expected_square = F(0), F(0)
        for x in product((-1, 1), repeat=n):
            probability = F(1)
            for value in x:
                probability *= (1+mean*value)/2
            score = a*a/(variance*variance)*sum(
                sign*(x[i]-mean)*(x[j]-mean)/2
                for sign, (i, j) in zip(signs, edges))
            expected_score += probability*score
            expected_square += probability*score*score
        assert expected_score == 0
        assert expected_square == attenuation**2*F(n-1, 2) == F(3, 98)
        checked += 1
    return {"all_order_four_signings": checked,
            "score_mean": "0", "score_second_moment": "3/98",
            "qualification": "exact Bernoulli score algebra with symbolic a=1/4; not a Gaussian-quantile certificate"}


def pair_information(z, correlation):
    density = lambda x: math.exp(-x*x/2)/math.sqrt(2*math.pi)
    delta = float(ndtr(-z))
    joint, error = quad(lambda x: density(x)*ndtr(
        (correlation*x-z)/math.sqrt(1-correlation*correlation)),
        z, math.inf, epsabs=1e-32, epsrel=1e-10, limit=150)
    single = delta-joint
    common = 1-2*delta+joint
    independent_rare = delta*delta
    information = (joint*math.log(joint/independent_rare)
                   +2*single*math.log1p((independent_rare-joint)/(delta*(1-delta)))
                   +common*math.log1p((joint-independent_rare)/(1-delta)**2))
    eta = density(z)**2/(delta*(1-delta))
    tangent_information = eta**2*correlation**2/2
    return {"z": z, "rho": correlation, "delta": delta,
            "pair_information": information,
            "quadratic_tangent": tangent_information,
            "ratio_to_quadratic_tangent": information/tangent_information,
            "quadrature_reported_error": error}


def main():
    diagnostics = [pair_information(z, rho) for z in (3., 4., 5., 6.)
                   for rho in (.1, .5)]
    print(json.dumps({"exact": exact_score_checks(),
                      "pair_resonance_diagnostics": diagnostics,
                      "diagnostic_scope": "numerical only; replicated-pair covariance is not a flat signing"},
                     indent=2))


if __name__ == '__main__':
    main()
