"""Focused adaptive monotone rectangle phase exclusion.

Default is a floating pilot ONLY; --verify uses directed mpmath intervals
on every accepted rectangle. The target is fixed strictly BELOW the exact
Gaussian value, so this proves a phase gap, not a decimal optimization.
"""
import argparse
import json
import math
import time
from fractions import Fraction
import mpmath as mp
from continued_feedback_ternary_latent_interval_certificate_2026_09_06 import sqrt_interval,log_interval

p=31/32
target=-.777671  # Strictly below g_4(1)=-.777670345...

def gain(A,B):
    if B==0 or p*B<=A:
        return 0.0
    if A==0 or p*B-A >=(1-p)*A*B:
        return p*math.log1p(B)-math.log1p(A)
    return p*math.log(p)+(1-p)*math.log1p(-p)-p*math.log(A)+math.log(B)-(1-p)*math.log(B-A)

def upper(rect):
    l,u,a,b=rect
    A=math.expm1(l*a*a/p)
    B=math.cosh(2*u*b/p)-1
    return .25*math.log(l*(8-l)/16)-l+gain(A,B)

def interval_upper(rect):
    iv=mp.iv
    l,u,a,b=map(iv.mpf,rect)
    pp=iv.mpf(31)/32
    A=(iv.exp(l*a*a/pp)-1).a
    z=2*u*b/pp
    B=((iv.exp(z)+iv.exp(-z))/2-1).b
    # Endpoints are exact dyadic numbers. Rational comparisons select the
    # constrained-v branch without any rounded sign decision.
    def rational(raw):
        sign,mantissa,exponent,_=raw
        value=Fraction(-mantissa if sign else mantissa)
        return value*(2**exponent) if exponent>=0 else value/Fraction(2**(-exponent))
    aa=rational(A._mpi_[0]);bb=rational(B._mpi_[1])
    if bb==0:
        G=iv.mpf(0)
    elif 31*bb<=32*aa:
        G=iv.mpf(0)
    elif aa==0 or 31*bb-32*aa>=aa*bb:
        G=pp*iv.ln(1+B)-iv.ln(1+A)
    else:
        G=pp*iv.ln(pp)+(1-pp)*iv.ln(1-pp)-pp*iv.ln(A)+iv.ln(B)-(1-pp)*iv.ln(B-A)
    return (.25*iv.ln(l*(8-l)/16)-l+G).b

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--verify',action='store_true')
    parser.add_argument('--max-boxes',type=int,default=2000000)
    args=parser.parse_args()
    # Interval endpoints have about 203 bits; 200 decimal digits make the
    # endpoint products/sign tests exact (well below 665 binary bits).
    mp.mp.dps=200
    mp.iv.dps=60
    slo,shi=sqrt_interval(257)
    rlo,rhi=(slo-1)/16,(shi-1)/16
    gaussian_lower=-4*(1-rlo)+log_interval(1-rhi*rhi)[0]/4
    assert gaussian_lower>Fraction(-777671,1000000)
    stack=[(.5,4.,0.,1.)]
    checked=accepted=0
    max_accepted=-math.inf
    start=time.monotonic()
    while stack:
        rect=stack.pop()
        bound=upper(rect)
        checked+=1
        if bound < target-1e-10:
            if args.verify:
                certified=interval_upper(rect)
                assert certified < mp.iv.mpf('-0.777671'), (rect,bound,certified)
            accepted+=1
            max_accepted=max(max_accepted,bound)
        else:
            l,u,a,b=rect
            if 4*(u-l)>=16*(b-a):
                mid=(l+u)/2
                stack.extend([(l,mid,a,b),(mid,u,a,b)])
            else:
                mid=(a+b)/2
                stack.extend([(l,u,a,mid),(l,u,mid,b)])
        if checked%100000==0:
            print(json.dumps(dict(checked=checked,accepted=accepted,pending=len(stack),
                                  seconds=time.monotonic()-start)),flush=True)
        if checked>=args.max_boxes:
            print(json.dumps(dict(status='budget stop, not certified',checked=checked,
                                  accepted=accepted,pending=len(stack))))
            return
    print(json.dumps(dict(status='directed interval exclusion' if args.verify else 'FLOATING PILOT ONLY',
                          checked=checked,accepted=accepted,max_float_upper=max_accepted,
                          threshold=target,gaussian_lower_rational=str(gaussian_lower),
                          seconds=time.monotonic()-start),indent=2))

if __name__=='__main__':
    main()
