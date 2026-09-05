"""Exact anchored Gaussian fixed-point lower certificate.

The intermediate centered bivariate polynomial has total Hermite degree
200. All arithmetic is exact rational arithmetic with outward intervals.
No infinite Hermite tail or numerical quadrature is used.
"""
from fractions import Fraction as F
from math import comb, factorial
import json

from fresh_limit_rooted_lower_certificate import I, DIGITS, INV_SQRT_2PI, phi, Phi
from fresh_limit_hierarchical_fixed_point_certificate import Phi_extended


def main():
    alpha = F(73, 100)
    rho = F(83, 100)
    resolvent = F(9, 2)
    degree = 200
    rho2 = rho*rho
    innovation_variance = 1-rho2
    density = phi(I(alpha))
    mass = 2*Phi(I(alpha))-1

    hermite = [F(1), alpha]
    for r in range(1, degree):
        hermite.append(alpha*hermite[-1]-r*hermite[-2])
    rho_powers = [rho2**r for r in range(degree+1)]
    noise_powers = [innovation_variance**r for r in range(degree+1)]
    grouped = [F(0) for _ in range(degree+1)]
    for r in range(2, degree+1, 2):
        base = hermite[r-1]**2/factorial(r)
        for ell in range(r+1):
            grouped[ell] += (base*comb(r,ell)*rho_powers[r-ell]
                             *noise_powers[ell])

    norm_rational = sum((q/(resolvent+ell)**2
                         for ell,q in enumerate(grouped)),F(0))
    derivative_rational = sum((ell*q/(resolvent+ell)**2
                               for ell,q in enumerate(grouped)),F(0))
    covariance_rational = sum((q/(resolvent+ell)
                               for ell,q in enumerate(grouped)),F(0))

    # The common positive factor 4*phi(alpha)^2 cancels exactly from D.
    derivative_energy = I(derivative_rational/norm_rational)
    norm_squared = 4*density*density*norm_rational
    centered_covariance = 2*density*covariance_rational/I(norm_rational).sqrt()
    source_mass = I((grouped[0]/resolvent**2)/norm_rational)
    covariance = rho*mass+I(innovation_variance).sqrt()*centered_covariance
    residual_variance = mass-covariance*covariance

    assert derivative_energy.hi < 1
    assert source_mass.lo > 0
    assert residual_variance.lo > 0
    residual_sd = residual_variance.sqrt()
    first_argument = covariance*alpha/residual_sd
    second_argument = alpha*mass.sqrt()/residual_sd
    value = (2*covariance*density*(2*Phi_extended(first_argument)-1)
             +4*mass.sqrt()*INV_SQRT_2PI*(1-Phi_extended(second_argument)))

    target = F(4283,10000)
    assert value.lo > target
    record = {
        'method':'exact_fraction_outward_intervals',
        'grid_decimal_digits':DIGITS,
        'alpha':'73/100',
        'anchor_correlation':'83/100',
        'resolvent_a':'9/2',
        'maximum_total_Hermite_degree':degree,
        'intermediate_response':'h_(k,l)=C*c_(k,l)/(a+l), 2<=k+l<=200 even',
        'mean_of_h':'0 exactly',
        'squared_norm_of_h':'1 exactly after the displayed normalization',
        'mask_mass':mass.json(),
        'unnormalized_squared_norm':norm_squared.json(),
        'conditional_derivative_energy':derivative_energy.json(),
        'strict_stability_margin':(1-derivative_energy).json(),
        'zero_innovation_degree_mass':source_mass.json(),
        'centered_threshold_covariance':centered_covariance.json(),
        'fixed_point_mask_covariance':covariance.json(),
        'residual_variance':residual_variance.json(),
        'first_cdf_argument':first_argument.json(),
        'second_cdf_argument':second_argument.json(),
        'lower_bound':value.json(),
        'target':'4283/10000',
        'margin_above_target':(value-target).json(),
        'verified':True,
    }
    rendered = json.dumps(record,indent=2)
    print(rendered)
    with open('computations/results/fresh_anchored_fixed_point_certificate.json',
              'w',encoding='utf-8') as handle:
        handle.write(rendered+'\n')


if __name__ == '__main__':
    main()
