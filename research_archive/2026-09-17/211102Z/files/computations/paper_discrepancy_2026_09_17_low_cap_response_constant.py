"""Exact rational constants for the no-operator-assumption response theorem.

The analytic inputs (uniform Gaussian-sign comparison, Grothendieck, and
the eventual half-range lower bound) are explicit in the audit artifact.
This script certifies arithmetic, not those imported analytic theorems.
"""
from fractions import Fraction as F
import json

C = F(1,2)
c0 = F(3,10)
ell = F(43,100)
groth = F(2)
epsilon = F(1,2**14)
alpha = F(1,2**23)
b0 = F(1,2**16)
t = F(1,2**13)
assert alpha/epsilon == F(1,512)
sqrt_epsilon = F(1,128)
sqrt_b0 = F(1,256)
assert sqrt_epsilon**2 == epsilon and sqrt_b0**2 == b0
# (5/sqrt(2))*sqrt(alpha)=5/4096, an exact rational simplification.
marked_cross = F(5,4096)
assert marked_cross**2 == F(25,2)*alpha
range_loss = F(3,2)*ell*alpha/epsilon
b_loss = sqrt_b0/(4*groth)
alpha_loss = marked_cross/(4*groth)
groth_loss = C*(sqrt_epsilon/2+epsilon)
total_loss = range_loss+b_loss+alpha_loss+groth_loss
lower_cap = ell+c0/(2*groth)-total_loss
assert lower_cap > C
# t*||H||<=1/2 and the Schur-arcsine remainder is at most half tb0.
assert 6*t <= F(1,2)
assert 216*t**3 <= t*b0/2
# pi<22/7 implies kappa^2>1/2 and 1/2-pi/12>5/21.
d_high_lower = F(5,42)*alpha
d_low_lower = t*b0/4
assert d_low_lower == F(1,2**31)
assert d_high_lower >= d_low_lower
# kappa>1/2, f(v)<=1-v^2/8.
discount_before_error = d_low_lower**2/16
assert discount_before_error == F(1,2**66)
discount_eventual = discount_before_error/2
# Elementary safe Grothendieck constant: asinh(1)>=5/6 and pi<22/7.
assert F(22,7)/(2*F(5,6)) < groth
print(json.dumps({
    "status":"PASS exact rational arithmetic",
    "cap_upper":str(C), "protected_energy_lower":str(c0),
    "eventual_half_range_lower":str(ell), "grothendieck_upper":str(groth),
    "epsilon":str(epsilon), "covariance_mass_threshold":str(alpha),
    "low_signed_covariance_frobenius_threshold":str(b0), "low_branch_t":str(t),
    "losses":{"range_deletion":str(range_loss),"low_frobenius":str(b_loss),
              "marked_cross":str(alpha_loss),"grothendieck_cross_and_tail":str(groth_loss)},
    "total_loss_exact":str(total_loss),
    "contradictory_cap_lower_exact":str(lower_cap),
    "positive_margin_exact":str(lower_cap-C),
    "discount_before_uniform_comparison_error":str(discount_before_error),
    "eventual_uniform_response_discount":str(discount_eventual),
    "scope":"Arithmetic certificate only; finite order cutoff supplied by uniform analytic limits."
},indent=2))
