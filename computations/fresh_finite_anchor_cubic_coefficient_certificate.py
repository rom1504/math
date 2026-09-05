"""Rational interval certificate for the banked mask's cubic edge response.

Run the finite-anchor certificate first.  This file uses its certified
Gaussian covariance intervals, not an optimization or numerical quadrature.
"""
from fractions import Fraction as F
import json

from fresh_limit_rooted_lower_certificate import I, DIGITS, INV_SQRT_2PI, phi
from fresh_limit_hierarchical_fixed_point_certificate import Phi_extended
from fresh_limit_mask_ascent_certificate import phi_large


def main():
    source = 'computations/results/fresh_finite_anchor_fixed_point_certificate.json'
    with open(source, encoding='utf-8') as handle:
        data = json.load(handle)
    assert data['verified'] and data['alpha'] == '361/500'
    assert data['rho'][0] == '2027/2500'

    def interval(key):
        item = data[key]
        return I(F(item['lower']), F(item['upper']))

    alpha, rho = F(361,500), F(2027,2500)
    p, w = interval('mask_mass'), interval('covariance')
    sigma = (p-w*w).sqrt()
    c = w/sigma
    k = (p-rho*w)/sigma
    # G0 = rho V + k Z + independent centered Gaussian;
    # W = w V + sigma Z.  Verify the residual variance is positive.
    residual = 1-rho*rho-k*k
    assert residual.lo > 0
    root = (1+c*c).sqrt()
    a = alpha*root
    tail = 1-Phi_extended(a)
    j0 = INV_SQRT_2PI*tail/root
    j2 = INV_SQRT_2PI*(a*phi_large(a)+tail)/(root*root*root)
    aa = 6*rho*rho*k-6*rho*k*k*c+2*k*k*k*c*c
    dd = -6*rho*rho*k-2*k*k*k
    value = (2/I(6).sqrt())*(
        rho**3*phi(I(alpha))*(alpha*alpha-1)*(2*Phi_extended(c*alpha)-1)
        +(2*c*rho**3+aa)*j2+(dd-2*c*rho**3)*j0)
    assert value.hi < F(-6,100)
    report = {
        'method': 'exact_fraction_outward_intervals_no_quadrature',
        'digits': DIGITS,
        'source_certificate': source,
        'coefficient': 'E[sign(W) 1{|V|>alpha} He_3(G0)/sqrt(6)]',
        'cubic_coefficient': value.json(),
        'conditional_G0_residual_variance': residual.json(),
        'conditional_V_variance_given_W': (1-w*w/p).json(),
        'verified_below_minus_0_06': True,
    }
    rendered = json.dumps(report, indent=2)
    print(rendered)
    with open('computations/results/fresh_finite_anchor_cubic_coefficient_certificate.json',
              'w', encoding='utf-8') as handle:
        handle.write(rendered+'\n')


if __name__ == '__main__':
    main()
