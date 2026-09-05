"""Exact rational interval certificate for 21 finite tree anchors.

Uses no optimizer, floating point, infinite polynomial tail, or quadrature.
"""
from collections import Counter
from fractions import Fraction as F
from math import comb, factorial
import json

from fresh_limit_rooted_lower_certificate import I, DIGITS, INV_SQRT_2PI, phi, Phi
from fresh_limit_hierarchical_fixed_point_certificate import Phi_extended


CHILDREN = [(), (0,0), (0,0,0,0), (0,1), (0,0,0,0,0,0),
            (0,0,0,1), (0,2), (0,3), (1,1), (0,0,0,0,0,0,0,0),
            (0,0,0,0,0,1), (0,0,0,2), (0,0,0,3), (0,0,1,1),
            (0,4), (0,5), (0,6), (0,7), (0,8), (1,2), (1,3)]


def main():
    alpha, a, degree = F(722,1000), F(17,5), 200
    rho = [F(t,10000) for t in
           [8108,-3110,1441,1662,-691,-1088,-758,-874,-444,
            331,638,496,572,502,357,563,392,452,230,286,330]]
    assert len(rho)==len(CHILDREN)==len(set(CHILDREN))
    assert all(len(ch)%2==0 and all(j<t for j in ch)
               for t,ch in enumerate(CHILDREN))
    total = sum((r*r for r in rho),F(0))
    innovation = 1-total
    assert innovation>0
    density=phi(I(alpha))
    mass=2*Phi(I(alpha))-1
    hermite=[F(1),alpha]
    for d in range(1,degree):
        hermite.append(alpha*hermite[-1]-d*hermite[-2])
    grouped=[F(0) for _ in range(degree+1)]
    for d in range(2,degree+1,2):
        base=hermite[d-1]**2/factorial(d)
        for ell in range(d+1):
            grouped[ell] += base*comb(d,ell)*total**(d-ell)*innovation**ell

    # The edge constant is excluded already by starting at total degree 2.
    covariance_base=rho[0]*mass
    removed=F(0)
    for t,ch in enumerate(CHILDREN[1:],1):
        mult=Counter(ch)
        denom=1
        monomial=F(1)
        for j,k in mult.items():
            denom *= factorial(k)
            monomial *= rho[j]**k
        raw=hermite[len(ch)-1]*monomial
        removed += raw*raw/denom
        coefficient=-2*density*raw/I(denom).sqrt()
        covariance_base += rho[t]*coefficient
    grouped[0] -= removed
    assert all(v>=0 for v in grouped)

    norm=sum((v/(a+ell)**2 for ell,v in enumerate(grouped)),F(0))
    derivative=sum((ell*v/(a+ell)**2 for ell,v in enumerate(grouped)),F(0))
    cov=sum((v/(a+ell) for ell,v in enumerate(grouped)),F(0))
    derivative_energy=I(derivative/norm)
    covariance=covariance_base+I(innovation).sqrt()*2*density*cov/I(norm).sqrt()
    residual_variance=mass-covariance*covariance
    assert derivative_energy.hi<1
    assert residual_variance.lo>0 and covariance.lo>0
    sd=residual_variance.sqrt()
    arg1=covariance*alpha/sd
    arg2=mass.sqrt()*alpha/sd
    value=(2*covariance*density*(2*Phi_extended(arg1)-1)
           +4*mass.sqrt()*INV_SQRT_2PI*(1-Phi_extended(arg2)))
    target=F(4306,10000)
    assert value.lo>target
    report={'method':'exact_fraction_outward_intervals','digits':DIGITS,
            'alpha':str(alpha),'resolvent_a':str(a),'degree':degree,
            'rho':[str(r) for r in rho],'anchor_children':CHILDREN,
            'innovation_variance':str(innovation),
            'conditional_derivative_energy':derivative_energy.json(),
            'strict_contraction_margin':(1-derivative_energy).json(),
            'mask_mass':mass.json(),'covariance':covariance.json(),
            'residual_variance':residual_variance.json(),
            'lower_bound':value.json(),'verified':True}
    rendered=json.dumps(report,indent=2)
    print(rendered)
    with open('computations/results/fresh_finite_anchor_fixed_point_certificate.json',
              'w',encoding='utf-8') as handle:
        handle.write(rendered+'\n')


if __name__=='__main__':
    main()
