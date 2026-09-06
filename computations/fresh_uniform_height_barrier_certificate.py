"""Exact elementary constants for uniform central-mask height separation."""
from fractions import Fraction as F
import json
from fresh_limit_rooted_lower_certificate import I, phi, Phi
from fresh_limit_hierarchical_fixed_point_certificate import Phi_extended


def verify():
    low_mass=2*Phi(I(F(1,5)))-1
    high_mass=2*Phi(I(F(9,10)))-1
    density=phi(I(F(9,10)))
    tail=1-Phi_extended(I(2))
    assert low_mass.hi<F(9,50)
    assert high_mass.lo>F(5,8)
    assert density.lo>F(1,4)
    assert tail.lo>F(1,50)
    width=F(1,1000)
    epsilon=width**2
    assert 1+F(9,10)*width<2
    noise_defect=width/F(200)
    kernel_gap=noise_defect-epsilon
    residual=kernel_gap/2
    assert residual==F(1,500000)
    print(json.dumps(dict(mask_mass_interval=['9/50','5/8'],
                         q=str(1-epsilon),kernel_gap_lower=str(kernel_gap),
                         height_residual_lower=str(residual),
                         gaussian_inequalities_verified=True),sort_keys=True))


if __name__=='__main__':
    verify()
