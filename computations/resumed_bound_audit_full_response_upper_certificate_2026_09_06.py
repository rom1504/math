#!/usr/bin/env python3
"""Exact rational upper certificate for the full ONE-response functional.

This does not upper-bound M_n. It bounds the variational lower-bound class.
All arithmetic is outward Fraction interval arithmetic; adaptive alpha
rectangles use proved separate monotonicity, not sampled interpolation.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import json
from pathlib import Path

from fresh_limit_rooted_lower_certificate import I
from fresh_limit_mask_ascent_certificate import cdf_point, phi_large


@lru_cache(None)
def endpoint(alpha):
    mass = 2*(1-cdf_point(alpha))
    first = 2*phi_large(alpha)
    return mass, first, first*first


def envelope(left, right, target):
    # h=1-p and u=r^2. F(h,p,u)=h Psi(sqrt(u/h),p-u).
    # For alpha in [left,right], h<=h(right), p<=p(left),
    # u<=m(left)^2. F is increasing in all three variables.
    p, m, u = endpoint(left)
    h = 1-endpoint(right)[0]
    t = p-u
    assert h.lo > 0 and h.hi <= 1 and t.lo > 0
    cauchy = (h*u+h*h*t).sqrt()
    if cauchy.hi < target:
        return cauchy, "Cauchy"
    k = m/h.sqrt()
    tau = t.sqrt()
    z = k/tau
    if z.hi > 8:
        return cauchy, "Cauchy"
    assert z.lo >= 0
    value = h*(k*(2*cdf_point(z)-1)+2*tau*phi_large(z))
    return value, "Gaussian"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--target", type=F, default=F(9,20))
    args = parser.parse_args()
    pending = [(F(0),F(3),0)]
    accepted = []
    evaluations = 0
    while pending:
        left,right,depth = pending.pop()
        assert depth < 30
        bound,method = envelope(left,right,args.target)
        evaluations += 1
        if bound.hi < args.target:
            accepted.append({"left":str(left),"right":str(right),
                             "method":method,"upper_bound":bound.json()})
        else:
            mid = (left+right)/2
            pending.extend([(mid,right,depth+1),(left,mid,depth+1)])
    accepted.sort(key=lambda row:F(row["left"]))
    assert accepted[0]["left"] == "0" and accepted[-1]["right"] == "3"
    assert all(a["right"] == b["left"] for a,b in zip(accepted,accepted[1:]))
    # At alpha>=3, h<=1 and the same monotone envelope is at most
    # E|m(3)+sqrt(p(3)-m(3)^2)N|<=sqrt(p(3)).
    tail = endpoint(F(3))[0].sqrt()
    assert tail.hi < args.target
    upper = max([F(row["upper_bound"]["upper"]) for row in accepted]+[tail.hi])
    assert upper < args.target
    result = {"method":"exact_fraction_monotone_alpha_envelopes",
              "scope":"full_final_one_response_functional_only_not_actual_signings",
              "verified":True,"target":str(args.target),
              "finite_intervals":len(accepted),"evaluations":evaluations,
              "uniform_upper_bound":I(upper).json(),
              "margin_below_target":(args.target-I(upper)).json(),
              "alpha_at_least_3_upper_bound":tail.json(),"intervals":accepted}
    report = {key:value for key,value in result.items() if key!="intervals"}
    print(json.dumps(report,indent=2),flush=True)
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+"\n")


if __name__ == "__main__":
    main()
