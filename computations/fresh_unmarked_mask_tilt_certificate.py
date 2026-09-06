"""Exact rational certificate of a strictly positive conference-only gain.

The base scalar integrals are read as previously certified outward intervals
from fresh_limit_mask_ascent_certificate.json. No numerical quadrature or
floating point is used here. The new proof uses Gaussian tilt tail formulas,
a certified rectangle contained in the mask, and elementary tail bounds.
"""

from fractions import Fraction as F
from math import factorial
from pathlib import Path
import json

from fresh_limit_rooted_lower_certificate import I, INV_SQRT_2PI
from fresh_limit_mask_ascent_certificate import Certificate, cdf_point, phi_large


def read_interval(data):
    return I(F(data["lower"]), F(data["upper"]))


def output_interval(x):
    return {"lower": str(x.lo), "upper": str(x.hi)}


def certify():
    base_path = Path(__file__).resolve().parent / "results/fresh_limit_mask_ascent_certificate.json"
    base = json.loads(base_path.read_text())
    assert base["verified"] is True
    assert base["alpha"] == "37/50" and base["Hermite_degree"] == 200
    cert = Certificate()
    p, w, s = cert.mass, cert.w, cert.s
    p1 = read_interval(base["new_mask_mass"])
    r1 = read_interval(base["covariance_V_new"])
    r0 = read_interval(base["covariance_W_new"])
    mv = 2*INV_SQRT_2PI*w/p.sqrt()-read_interval(base["new_V_signed_moment"])
    mw = 2*INV_SQRT_2PI*p.sqrt()-read_interval(base["new_absolute_moment"])
    c, d = (p*r1-w*r0)/(s*s), (r0-w*r1)/(s*s)
    af, bf = (p*mv-w*mw)/(s*s), (mw-w*mv)/(s*s)
    residual_squared = 1-p1-af*mv-bf*mw
    assert 0 <= c.lo <= c.hi < F(9,50)
    assert 0 <= d.lo <= d.hi < F(3,4)
    assert -F(1,2) < af.lo <= af.hi < F(1,2)
    assert -F(1,5) < bf.lo <= bf.hi < F(1,5)
    assert residual_squared.lo > F(9,100)

    # The mask contains {|V|<=17/20, |W|<=1/6}. Ignore its nonnegative
    # b*H0 term; certify g>=13/20 on the whole V interval by Taylor bins.
    assert cert.a.lo > F(149,250)
    assert cert.b.lo > 0
    assert F(5,2)*F(149,250)*F(13,20)-F(4,5) > F(1,6)
    min_g = F(10)
    for k in range(34):
        lo, hi = F(k,40), F(k+1,40)
        center, radius = (lo+hi)/2, (hi-lo)/2
        enclosure = cert.g_range(cert.g_at(center), radius).v
        min_g = min(min_g, enclosure.lo)
        assert enclosure.lo > F(13,20)

    tau, v0, t0 = F(1,125), F(17,20), F(1,6)
    sv = ((s*s+tau*tau)/(p+tau*tau)).sqrt()
    sw = (p*tau*tau/(p+tau*tau)).sqrt()
    normalizer = INV_SQRT_2PI*tau/(p+tau*tau).sqrt()
    alpha = v0/sv
    prob_v = 2*(1-cdf_point(alpha))
    # sw<=tau. The two-sided Gaussian tail beyond t0/sw is <10^-50.
    # Indeed 2Phi(-x)<=exp(-x²/2); lower-bound exp((t0/tau)²/2)
    # by a finite positive Taylor sum, using exact integers only.
    exponent = (t0/tau)**2/2
    exp_lower = sum((exponent**k/F(factorial(k)) for k in range(201)), F(0))
    assert exp_lower > 10**50
    tail = I(F(1,10**50))
    tail_sqrt = I(F(1,10**25))

    density_upper = normalizer*(prob_v+tail)
    v_abs_upper = normalizer*(2*sv*phi_large(alpha)+sv*tail_sqrt)
    v_squared_upper = normalizer*(
        2*sv*sv*(alpha*phi_large(alpha)+1-cdf_point(alpha))
        +I(3).sqrt()*sv*sv*tail_sqrt)
    w_abs_upper = normalizer*sw*(prob_v.sqrt()+tail_sqrt)
    density_lower = normalizer*2*(1-cdf_point(1/sv))

    feature_correlation_lower = I(F(3,10))-(
        density_upper+F(1,2)*v_abs_upper+F(1,5)*w_abs_upper
    )/(F(3,10)*INV_SQRT_2PI)
    assert feature_correlation_lower.lo > 0
    paired_loss_upper = 2*tau*(
        F(9,50)*w/(s*s)*v_squared_upper+F(3,4)*density_upper)
    gain = 2*density_lower*feature_correlation_lower-paired_loss_upper
    assert gain.lo > F(1,10**6)
    baseline = read_interval(base["lower_bound"])
    return {
        "method": "exact_fraction_gaussian_tilt_and_mask_rectangle",
        "scope": "conference_or_explicit_all_power_delocalization_only",
        "tau": str(tau),
        "rectangle": {"V_bound": str(v0), "W_bound": str(t0)},
        "certified_min_g": str(min_g),
        "response_residual_variance": output_interval(residual_squared),
        "density_lower": output_interval(density_lower),
        "density_upper": output_interval(density_upper),
        "feature_correlation_lower": output_interval(feature_correlation_lower),
        "paired_loss_upper": output_interval(paired_loss_upper),
        "gain": output_interval(gain),
        "gain_target": "1/1000000",
        "base_lower_bound": output_interval(baseline),
        "conference_lower_bound": output_interval(baseline+gain),
        "verified": True,
    }


if __name__ == "__main__":
    print(json.dumps(certify(), indent=2))
