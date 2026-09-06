"""Exact T lower witness separated from the archived exact E upper."""
import json
from fractions import Fraction as F
from pathlib import Path
from continued_feedback_ternary_latent_interval_certificate_2026_09_06 import log_interval, entropy_interval, sqrt_interval


def gaussian_lower(v, t):
    if v == 0:
        return F(0)
    x = t * v
    rho = (sqrt_interval(1 + 16 * x * x)[0] - 1) / (4 * x)
    assert 0 <= rho < 1
    return -x * (1-rho) + log_interval(1-rho*rho)[0] / 4


def main():
    p,t = F(31,32),F(4)
    z1,z2 = F(1363,1500),F(374,375)
    w = (z2-p)/(z2-z1)
    v1,v2 = z1/p,(z2-z2*z2)/p
    q1 = entropy_interval(z1)[0] + z1*log_interval(2)[0] + gaussian_lower(v1,t)
    q2 = entropy_interval(z2)[0] + gaussian_lower(v2,t)
    lower = w*q1 + (1-w)*q2 - entropy_interval(p)[1] + t*(1-sqrt_interval(p)[1])
    e_upper = F(-19678127864847,800000000000000)
    gap = lower-e_upper
    assert gap>F(9,1000)
    result = {"status":"exact rational lower-channel certificate against archived exact E upper",
              "p":str(p),"t":str(t),"posterior_z":[str(z1),str(z2)],
              "posterior_bias":["0","1"],"first_orbit_weight":str(w),
              "variances":[str(v1),str(v2)],"T_offset_lower":str(lower),
              "T_offset_lower_display":float(lower),"E_offset_upper":str(e_upper),
              "gap_lower":str(gap),"gap_lower_display":float(gap),
              "gap_exceeds":"9/1000",
              "E_dependency":"continued_feedback_ternary_latent_interval_certificate_2026_09_06.py default certificate"}
    Path("computations/decisive_bridge_supersolution_ternary_strict_certificate_2026_09_06.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__ == "__main__":
    main()
