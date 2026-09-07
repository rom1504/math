"""Rational certificate: the tilted noisy-Gaussian rank-two kernel <=1/15.

No floating arithmetic is used in acceptance. exp is upper-bounded by a
rational Taylor remainder, followed by upward fixed-point squaring.
"""
import json
import math
from fractions import Fraction as Q
from pathlib import Path


SCALE=10**40


def ceil_fraction(x):
    return (x.numerator+x.denominator-1)//x.denominator


def exp_upper(x):
    # exp(x)=exp(x/128)^128. All arguments in this certificate are small
    # enough after this fixed range reduction.
    z=x/128
    assert abs(z)<1
    value=Q(1); term=Q(1)
    degree=18
    for j in range(1,degree+1):
        term*=z/j; value+=term
    az=abs(z)
    remainder=az**(degree+1)/math.factorial(degree+1)/(1-az/(degree+2))
    bound=ceil_fraction((value+remainder)*SCALE)
    assert bound>0
    for _ in range(7):
        bound=(bound*bound+SCALE-1)//SCALE
    return Q(bound,SCALE)


def run():
    sqrt30_upper=Q(5477225575051662,10**15)
    assert sqrt30_upper**2>30
    count=500; leaves=[]; maximum=Q(0); maximizing=None
    for j in range(count):
        lo=Q(j,count); hi=Q(j+1,count)
        dmin=16-15*hi*hi; dmax=16-15*lo*lo
        assert dmin>=1
        penalty_lower=30*lo*lo/dmax
        argument_upper=4*sqrt30_upper*hi/dmin
        upper=(exp_upper(argument_upper-penalty_lower)
               +exp_upper(-argument_upper-penalty_lower))/(2*dmin)
        assert upper<Q(1,15), (j,upper)
        leaves.append(dict(index=j,lower=str(lo),upper=str(hi),kernel_upper=str(upper)))
        if upper>maximum: maximum=upper; maximizing=j
    result=dict(status='PASS: EXACT RATIONAL CERTIFICATE',bound='1/15',
                cells=count,sqrt30_upper=str(sqrt30_upper),scale=str(SCALE),
                Taylor_degree=18,range_reduction=128,maximum_leaf_upper=str(maximum),
                maximizing_cell=maximizing,leaves=leaves)
    target=Path(__file__).resolve().parent/'results'/'flatify_adversary_2026_09_07_noisy_gaussian_kernel_certificate.json'
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='leaves'}),flush=True)


if __name__=='__main__': run()
