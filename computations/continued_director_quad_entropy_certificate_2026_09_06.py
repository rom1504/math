"""Rational interval certificate for the four-point stable-profile constant.

All inequalities establishing negativity use fractions, not floating point.
The script also checks the projection fibres directly on all 81 patterns.
"""
from collections import defaultdict
from fractions import Fraction as F
from itertools import product
from math import factorial, isqrt
from pathlib import Path
import argparse
import json

D = 10 ** 25


def enclose(x):
    a = x.numerator * D // x.denominator
    return F(a, D), F(a + 1, D)


def unit_log(x, terms=35):
    assert 1 <= x <= 2
    y = (x - 1) / (x + 1)
    lower = 2 * sum((y ** (2*j+1) / (2*j+1) for j in range(terms)), F(0))
    tail = 2*y**(2*terms+1)/((2*terms+1)*(1-y*y))
    return lower, lower+tail


def log_bound(x):
    assert x > 0
    power = 0
    while x < 1:
        x *= 2
        power -= 1
    while x >= 2:
        x /= 2
        power += 1
    lo, hi = unit_log(x)
    l2, u2 = unit_log(F(2))
    if power >= 0:
        lo, hi = lo+power*l2, hi+power*u2
    else:
        lo, hi = lo+power*u2, hi+power*l2
    return enclose(lo)[0], enclose(hi)[1]


def exp_bound(x, terms=90):
    assert x >= 0
    lower = sum((x**j/factorial(j) for j in range(terms+1)), F(0))
    ratio = x/(terms+2)
    assert ratio < 1
    upper = lower+x**(terms+1)/factorial(terms+1)/(1-ratio)
    return enclose(lower)[0], enclose(upper)[1]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', default='computations/results/continued_director_quad_entropy_certificate_2026_09_06.json')
    args = ap.parse_args()
    p, q, t = F(15,16), F(1,16), F(33,32)
    probs = {-1:p/2, 0:q, 1:p/2}
    law, fibres = defaultdict(F), defaultdict(list)
    for x in product((-1,0,1), repeat=4):
        y = (sum(x), x[0]+x[1]-x[2]-x[3], x[0]-x[1]+x[2]-x[3])
        prob = F(1)
        for a in x:
            prob *= probs[a]
        law[y] += prob
        fibres[y].append(x)
    assert sum(law.values()) == 1
    chi = (1,-1,-1,1)
    for xs in fibres.values():
        for x in xs:
            diff = tuple(x[i]-xs[0][i] for i in range(4))
            assert all(diff[i] == diff[0]*chi[i] for i in range(4))
    hlo = hhi = F(0)
    for a in law.values():
        lo, hi = log_bound(a)
        hlo -= a*hi
        hhi -= a*lo
    l2,u2 = log_bound(F(2))
    clo,chi_bound = (p+F(3,4))*l2-hhi/4, (p+F(3,4))*u2-hlo/4
    # e^(-4t), then log[(1+e^(-4t))/2].
    elo,ehi = exp_bound(4*t)
    klo,khi = enclose((1+1/ehi)/2)[0], enclose((1+1/elo)/2)[1]
    lklo,lkhi = log_bound(klo)[0], log_bound(khi)[1]
    root_floor = isqrt(p.numerator*D*D//p.denominator)
    rlo,rhi = F(root_floor,D), F(root_floor+1,D)
    assert rlo*rlo <= p <= rhi*rhi
    lower = clo+lklo/2+t*(1-rhi)
    upper = chi_bound+lkhi/2+t*(1-rlo)
    assert upper < F(-3,100000)
    def interval(lo,hi):
        return {'lower_rational':str(enclose(lo)[0]),
                'upper_rational':str(enclose(hi)[1]),
                'display_lower':float(lo), 'display_upper':float(hi)}
    result = {'status':'proved rational interval inequalities; no full cap claim',
              'p':str(p),'t':str(t),'enumerated_patterns':81,
              'projection_values':len(law),'entropy':interval(hlo,hhi),
              'count_rate':interval(clo,chi_bound),'tilted_class_exponent':interval(lower,upper),
              'exact_assertion':'tilted class exponent < -3/100000',
              'scope':'epsilon->0 nearly-flat class only; this tilt fails the Gaussian full-sum test'}
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
