"""Exact scalar certificate for a strict ceiling on the tree-mask method.

This is NOT an upper bound on the original Boolean minimax constant.
It certifies only the analytic consequence of the tree-height residual and
Gaussian rearrangement-stability arguments.
"""
from fractions import Fraction as F
from math import factorial
import json

from fresh_limit_rooted_lower_certificate import I, DIGITS, INV_SQRT_2PI, phi, Phi


def main():
    alpha = F(27,40)  # .675 exceeds the median absolute-Gaussian threshold.
    q = F(19,20)
    degree = 400
    density = phi(alpha)
    mass = 2*Phi(alpha)-1
    assert mass.lo > F(1,2)
    hermite = [F(1),alpha]
    for r in range(1,degree):
        hermite.append(alpha*hermite[-1]-r*hermite[-2])
    series = sum((q**r*hermite[r-1]**2/factorial(r)
                  for r in range(2,degree+1,2)),F(0))
    partial = mass*mass+4*density*density*series
    tail = q**(degree+2)*mass*(1-mass)
    kernel_upper = ((partial+tail)/mass).hi
    assert 0 < kernel_upper < q
    kernel = I(kernel_upper)
    defect_squared = 2-2*(q*kernel).sqrt()-2*((1-q)*(1-kernel)).sqrt()
    defect_floor = F(37,200)
    assert defect_squared.lo > defect_floor**2

    # Unique maximum of R(p)=2 sqrt(p) phi(Phi^-1((1+p)/2)).
    root_lo,root_hi = F(65730655,10**8),F(65730656,10**8)
    for t,expected in ((root_lo,-1),(root_hi,1)):
        derivative_equation = t*(2*Phi(t)-1)-phi(t)
        assert derivative_equation.hi < 0 if expected < 0 else derivative_equation.lo > 0
    root = I(root_lo,root_hi)
    root_mass = 2*Phi(root)-1
    assert F(12,25) < root_mass.lo < root_mass.hi < F(1,2)
    envelope = 2*root_mass.sqrt()*phi(root)

    left_threshold = F(6433,10000)
    right_threshold = F(6744,10000)
    assert (2*Phi(left_threshold)-1).hi < F(12,25)
    assert (2*Phi(right_threshold)-1).hi < F(1,2)
    left_cap = 2*I(F(12,25)).sqrt()*phi(left_threshold)
    right_cap = 2*I(F(1,2)).sqrt()*phi(right_threshold)

    # Uniform rearrangement deficit on p in [.48,.5].
    p_min = F(12,25)
    gap = I(p_min).sqrt()*p_min*p_min*defect_floor**4/(8*INV_SQRT_2PI)
    assert (2*phi(alpha)).lo > (F(1,4)*defect_floor**4/(8*INV_SQRT_2PI)).hi
    inside_cap = envelope-gap
    upper = max(left_cap.hi,right_cap.hi,inside_cap.hi)
    target = F(899,2000)
    assert upper < target
    record = {
        'method':'exact_fraction_outward_intervals',
        'grid_decimal_digits':DIGITS,
        'scope':'ceiling of the hierarchical one-mask certificate, not original minimax upper bound',
        'noise_kernel_threshold':'27/40',
        'height_test_correlation':'19/20',
        'Hermite_degree':degree,
        'noise_kernel_upper':I(kernel_upper).json(),
        'residual_squared_lower':defect_squared.json(),
        'uniform_residual_floor':'37/200',
        'scalar_envelope_maximum':envelope.json(),
        'left_mass_region_cap':left_cap.json(),
        'right_mass_region_cap':right_cap.json(),
        'middle_region_deficit':gap.json(),
        'middle_region_cap':inside_cap.json(),
        'global_method_cap':I(upper).json(),
        'strict_target':'899/2000',
        'verified':True,
    }
    rendered = json.dumps(record,indent=2)
    print(rendered)
    with open('computations/results/fresh_tree_mask_ceiling_certificate.json','w',encoding='utf-8') as handle:
        handle.write(rendered+'\n')


if __name__ == '__main__':
    main()
