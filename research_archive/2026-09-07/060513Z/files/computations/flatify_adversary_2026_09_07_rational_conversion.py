"""Rational-only logarithm/square-root endpoint reconstruction; no mpmath."""
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import json

def positive_log_bounds(x,terms=80):
    assert x>0
    k=0
    while x<1:x*=2;k-=1
    while x>=2:x/=2;k+=1
    def series(y):
        z=(y-1)/(y+1)
        s=2*sum((z**(2*j+1))/F(2*j+1) for j in range(terms))
        remainder=2*z**(2*terms+1)/(F(2*terms+1)*(1-z*z))
        return s,s+remainder
    lo,hi=series(x)
    l2,h2=series(F(2))
    return (lo+k*l2,hi+k*h2) if k>=0 else (lo+k*h2,hi+k*l2)

def sqrt_bounds(x,digits=40):
    scale=10**digits
    k=isqrt(x.numerator*scale*scale//x.denominator)
    lo,hi=F(k,scale),F(k+1,scale)
    assert lo*lo<=x<=hi*hi
    return lo,hi

p,t,E=F(24,25),F(97,20),-F(5151,6250)
l2,h2=positive_log_bounds(F(2))
pl,pu=sqrt_bounds(p)
numerator_lo,numerator_hi=t+p*l2+E,t+p*h2+E
assert numerator_lo>0
assert p*h2+E<0
caplo,caphi=numerator_lo/(2*t*pu),numerator_hi/(2*t*pl)
endpoint=F(493608094,10**9)
assert caphi<endpoint
sl,su=sqrt_bounds(1+16*t*t)
rl,ru=(sl-1)/(4*t),(su-1)/(4*t)
assert 0<rl<ru<1
loglo=positive_log_bounds(1-ru*ru)[0]
loghi=positive_log_bounds(1-rl*rl)[1]
glo,ghi=-t*(1-rl)+loglo/4,-t*(1-ru)+loghi/4
assert ghi<E
assert t*(1-rl)<F(1,2)
# Store concise decimal-rational enclosing endpoints plus rigorous gaps.
scale=10**30
floor=lambda x:F(x.numerator*scale//x.denominator,scale)
ceil=lambda x:-floor(-x)
out=dict(status='exact rational PASS',p=str(p),t=str(t),E_upper=str(E),
         cap_enclosure=list(map(str,(floor(caplo),ceil(caphi)))),
         published_endpoint=str(endpoint),gap_lower=str(floor(endpoint-caphi)),
         gaussian_enclosure=list(map(str,(floor(glo),ceil(ghi)))),
         gaussian_target_gap_lower=str(floor(E-ghi)),
         lambda_star_enclosure=list(map(str,(floor(t*(1-ru)),ceil(t*(1-rl))))),
         logarithm_terms=80,sqrt_decimal_digits=40)
(Path(__file__).resolve().parent/'results'/(Path(__file__).stem+'.json')).write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
