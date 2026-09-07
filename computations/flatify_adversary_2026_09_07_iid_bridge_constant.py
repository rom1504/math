"""Rational certificate for constants in the analytic iid-bridge floor.

Uses only pi<22/7 and integer/Fraction arithmetic. It does not numerically
solve a Parisi equation or estimate a spin-glass ground-state constant.
"""
from fractions import Fraction as F
from math import isqrt
from decimal import Decimal,localcontext
from pathlib import Path
import json


def sqrt_interval(x,digits=40):
    scale=10**digits
    k=isqrt(x.numerator*scale*scale//x.denominator)
    lo=F(k,scale); hi=F(k+1,scale)
    assert lo*lo<=x<hi*hi
    return lo,hi


def dec(x):
    with localcontext() as c:
        c.prec=45
        return str(Decimal(x.numerator)/Decimal(x.denominator))


def run():
    sk_lower_squared=F(2240**2*7**7,57**2*22**7)
    certified_sk=F(71,100)
    assert sk_lower_squared>certified_sk**2
    exact_integer_gap=2240**2*7**7*100**2-71**2*57**2*22**7
    assert exact_integer_gap==468978757537408
    bridge_lower_squared=sk_lower_squared/2
    safe_bridge_squared=certified_sk**2/2
    target=F(493608094,10**9)
    assert safe_bridge_squared>F(1,2)**2>target**2
    result=dict(status='PASS: EXACT RATIONAL CONSTANT CERTIFICATE',
        pi_upper='22/7',sk_lower_squared=str(sk_lower_squared),
        certified_sk_lower='71/100',exact_integer_gap=exact_integer_gap,
        sk_lower_interval=[dec(t) for t in sqrt_interval(sk_lower_squared)],
        bridge_lower_interval=[dec(t) for t in sqrt_interval(bridge_lower_squared)],
        safe_bridge_interval=[dec(t) for t in sqrt_interval(safe_bridge_squared)],
        safe_bridge_squared_gap_above_one_half_squared=str(safe_bridge_squared-F(1,4)),
        favorable_upper_target=str(target),
        safe_bridge_squared_gap_above_target_squared=str(safe_bridge_squared-target**2))
    print(json.dumps(result,indent=2))
    return result


if __name__=='__main__':
    result=run()
    target=Path(__file__).resolve().parent/'results'/'flatify_adversary_2026_09_07_iid_bridge_constant.json'
    target.write_text(json.dumps(result,indent=2)+'\n')
