"""Exact outward constant conversion; phase proof has a separate certificate."""
import json
from fractions import Fraction as F
from continued_feedback_ternary_latent_interval_certificate_2026_09_06 import sqrt_interval, log_interval

p = F(31, 32)
slo, shi = sqrt_interval(257)
rlo, rhi = (slo-1)/16, (shi-1)/16
glo = -4*(1-rlo) + log_interval(1-rhi*rhi)[0]/4
ghi = -4*(1-rhi) + log_interval(1-rlo*rlo)[1]/4
plo, phi = sqrt_interval(p)
llo, lhi = log_interval(2)
# C=(4+p log(2)+g_4(1))/(8 sqrt(p)). All numerators positive.
clo = (4+p*llo+glo)/(8*phi)
chi = (4+p*lhi+ghi)/(8*plo)
assert chi < F(494515125, 1000000000)
print(json.dumps(dict(status="exact rational outward bounds", lower=str(clo),
                     upper=str(chi), lower_float=float(clo), upper_float=float(chi),
                     g_lower=str(glo), g_upper=str(ghi)), indent=2))
