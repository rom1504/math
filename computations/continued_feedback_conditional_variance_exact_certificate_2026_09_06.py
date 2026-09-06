"""Exact integer/rational upper certificate for the conditional-variance envelope."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import numpy as np

from continued_feedback_ternary_latent_interval_certificate_2026_09_06 import (
    SCALE,entropy_interval,log_interval,sqrt_interval,hull_at)


def ceil_scaled(x):
    return -((-x.numerator*SCALE)//x.denominator)


def gaussian_upper(variance,t):
    if variance==0:return F(0)
    x=t*variance
    lo,hi=sqrt_interval(1+16*x*x)
    rlo=(lo-1)/(4*x);rhi=(hi-1)/(4*x)
    assert 0<=rlo<=rhi<1
    return -x*(1-rhi)+log_interval(1-rlo*rlo)[1]/4


def certificate(n=2500,variance_grid=1000):
    p=F(31,32);t=F(4)
    assert n>=2 and variance_grid>=1
    assert 32*variance_grid*n**4<2**63
    assert n*SCALE<2**63
    h_z=np.array([ceil_scaled(entropy_interval(F(i,n))[1]) for i in range(n+1)],dtype=np.int64)
    h_s=np.array([ceil_scaled(entropy_interval(F(n+j,2*n))[1]) for j in range(n+1)],dtype=np.int64)
    print(json.dumps({'stage':'entropy_tables_ready','n':n}),flush=True)
    last=(variance_grid*p.denominator+p.numerator-1)//p.numerator
    gs=np.array([ceil_scaled(gaussian_upper(F(i,variance_grid),t)) for i in range(last+1)],dtype=np.int64)
    print(json.dumps({'stage':'gaussian_table_ready','size':len(gs)}),flush=True)
    js=np.arange(n+1,dtype=np.int64);jsquare=js*js
    denominator=p.numerator*n**4
    values=[];maxbin=0
    for i in range(n+1):
        numerator=(i*n**3-i*i*jsquare)*(p.denominator*variance_grid)
        bins=numerator//denominator
        assert int(bins.min())>=0 and int(bins.max())<=last
        maxbin=max(maxbin,int(bins.max()))
        value=h_z[i]+(i*h_s+n-1)//n+gs[bins]
        values.append(F(int(value.max()),SCALE))
    hull=hull_at(values,p)
    delta=F(1,n)
    modulus=entropy_interval(delta)[1]+delta*log_interval(2)[1]
    modulus+=entropy_interval(delta/4)[1]+2*t*delta/p
    offset=-entropy_interval(p)[0]+t*(1-sqrt_interval(p)[0])
    point=hull+offset;answer=point+modulus
    assert answer<0
    return {'target':'conditional-variance latent envelope offset at p31/32,t4',
        'p':str(p),'t':str(t),'z_and_s_denominator':n,'variance_grid':variance_grid,
        'maximum_variance_bin':maxbin,'exact_integer_grid':True,
        'grid_hull_offset_upper':str(point),'grid_hull_offset_display':float(point),
        'continuous_modulus_upper':str(modulus),'continuous_modulus_display':float(modulus),
        'certified_offset_upper':str(answer),'certified_offset_display':float(answer),
        'strictly_negative':True,
        'scope':'exact envelope upper; its signing consequence requires the separately proved supersolution and stopped-tree theorems'}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--n',type=int,default=2500);parser.add_argument('--variance-grid',type=int,default=1000)
    args=parser.parse_args();result=certificate(args.n,args.variance_grid)
    print(json.dumps(result,indent=2),flush=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
