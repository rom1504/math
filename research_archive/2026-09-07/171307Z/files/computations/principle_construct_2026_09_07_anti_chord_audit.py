"""Independent entropy-tangent enclosure of the anti-weave density chords.

Conditional on the listed rational E bounds; does not certify those bounds.
No optimized critical-point formula is used here.
"""
from fractions import Fraction as F
import json
import mpmath as mp
from flatify_independent_2026_09_07_ternary_interval import rational_iv, endpoints

mp.iv.dps = 55
p, t = F(24, 25), F(7, 2)
raw = [0, -224, -378, -465, -529, -581, -623, -660, -691,
       -720, -745, -768, -789, -808, -826, -843, -859,
       -874, -888, -901, -832]
values = [F(0)] + [F(x, 1000) + F(1, 2000) for x in raw[1:]]
alpha = rational_iv(p / 2)
best = None
records = []
for j in range(10):
    left, right = F(j, 20), F(j + 1, 20)
    el = (values[j] + values[20-j]) / 2
    er = (values[j+1] + values[19-j]) / 2
    slope = (er-el)/(right-left)
    intercept = el-slope*left
    interval_best = None
    for k in range(200):
        lo = left+(right-left)*F(k, 200)
        hi = left+(right-left)*F(k+1, 200)
        mid = (lo+hi)/2
        z = rational_iv(mid)
        entropy = -z*mp.iv.ln(z)-(1-z)*mp.iv.ln(1-z)
        derivative = mp.iv.ln((1-z)/z)
        # Concavity gives h(x)<=h(mid)+h'(mid)(x-mid).
        tangent_values = [alpha*(entropy+derivative*rational_iv(x-mid))
                          +rational_iv(slope*x+intercept)
                          for x in (lo, hi)]
        upper = max(endpoints(v)[1] for v in tangent_values)
        interval_best = upper if interval_best is None else max(interval_best, upper)
    exponent = rational_iv(interval_best)+alpha*mp.iv.ln(2)
    cap = (rational_iv(t)+exponent)/(2*rational_iv(t)*mp.iv.sqrt(rational_iv(p)))
    cap_upper = endpoints(cap)[1]
    assert cap_upper < F(499, 1000)
    best = cap_upper if best is None else max(best, cap_upper)
    records.append({"theta_cell": j, "cap_upper": str(cap_upper),
                    "display_only": float(cap_upper)})
print(json.dumps({"status": "CONDITIONAL_CHORD_TANGENT_AUDIT_PASS",
                  "requires": "all 20 rational E point bounds",
                  "max_cap_upper": str(best),
                  "display_only": float(best), "cells": records}, indent=2))
