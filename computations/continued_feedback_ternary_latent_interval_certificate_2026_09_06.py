"""Exact rational upper certificate for the TERNARY LATENT LOWER ENVELOPE.

No optimizer and no floating transcendental value enter the certificate.
Its negative result is NOT an upper bound on the recursive Bellman value.
"""
import argparse
from bisect import bisect_left, bisect_right
from fractions import Fraction as F
import json
from math import isqrt
from pathlib import Path


SCALE=10**12


def log_interval(x):
    x=F(x)
    if x<=0:raise ValueError('positive argument required')
    power=0
    while x<1:x*=2;power-=1
    while x>=2:x/=2;power+=1
    def unit(y):
        z=(y-1)/(y+1)
        partial=2*sum(z**(2*j+1)/F(2*j+1) for j in range(25))
        remainder=2*z**51/(51*(1-z*z))
        return partial,partial+remainder
    low,high=unit(x);l2,u2=unit(F(2))
    low+=power*(l2 if power>=0 else u2)
    high+=power*(u2 if power>=0 else l2)
    lower=(low.numerator*SCALE)//low.denominator
    upper=-((-high.numerator*SCALE)//high.denominator)
    return F(lower,SCALE),F(upper,SCALE)


def entropy_interval(x):
    lo=F(0);hi=F(0)
    for probability in (F(x),1-F(x)):
        if probability:
            a,b=log_interval(probability)
            lo-=probability*b;hi-=probability*a
    return lo,hi


def sqrt_interval(x):
    x=F(x);a=isqrt((x.numerator*SCALE*SCALE)//x.denominator)
    lo=F(a,SCALE);hi=F(a+1,SCALE)
    assert lo*lo<=x<=hi*hi
    return lo,hi


class InnerCertificate:
    def __init__(self,grid):
        self.grid=grid;self.log2=log_interval(2)
        self.threshold_low=[F(1,2)];self.threshold_high=[F(1,2)]
        for j in range(1,grid):
            s=F(j,grid)
            lo,hi=log_interval((1+s)/(1-s))
            self.threshold_low.append(lo/(4*s));self.threshold_high.append(hi/(4*s))
        assert all(a<=b for a,b in zip(self.threshold_low,self.threshold_low[1:]))
        assert all(a<=b for a,b in zip(self.threshold_high,self.threshold_high[1:]))
        # Midpoint entropy tangents; index is twice the midpoint's grid index.
        self.tangents={}
        for j in range(1,2*grid):
            midpoint=F(j,2*grid)
            hlo,hhi=entropy_interval((1+midpoint)/2)
            alo,ahi=log_interval((1+midpoint)/(1-midpoint))
            self.tangents[j]=(hhi,-ahi/2,-alo/2)
        self.maximum_bracket_width=F(0)

    def upper(self,alpha):
        alpha=F(alpha)
        if alpha<=F(1,2):return self.log2[1]
        lower=bisect_right(self.threshold_high,alpha)-1
        upper=bisect_left(self.threshold_low,alpha)
        if upper>=self.grid:upper=self.grid
        upper=max(upper,lower+1)
        assert self.threshold_high[lower]<=alpha
        assert upper==self.grid or alpha<=self.threshold_low[upper]
        l=F(lower,self.grid);u=F(upper,self.grid);mid=(l+u)/2
        self.maximum_bracket_width=max(self.maximum_bracket_width,u-l)
        entropy_upper,derivative_lower,derivative_upper=self.tangents[lower+upper]
        left=entropy_upper+derivative_lower*(l-mid)+alpha*l*l
        right=entropy_upper+derivative_upper*(u-mid)+alpha*u*u
        return max(left,right)


def hull_at(values,p):
    count=len(values)-1;stack=[]
    for index in range(count+1):
        while len(stack)>=2:
            a,b=stack[-2:]
            if (values[b]-values[a])*(index-b) <= (values[index]-values[b])*(b-a):
                stack.pop()
            else:break
        stack.append(index)
    target=p*count
    for a,b in zip(stack,stack[1:]):
        if a<=target<=b:
            return ((b-target)*values[a]+(target-a)*values[b])/(b-a)
    raise AssertionError('mean not covered')


def certificate(z_grid=1000,kappa_denominator=100,spin_grid=100):
    assert z_grid>=2 and spin_grid>=2
    assert kappa_denominator>=2 and kappa_denominator%2==0
    p=F(31,32);t=F(4);maximum_kappa=t/p
    inner=InnerCertificate(spin_grid)
    entropies=[entropy_interval(F(j,z_grid))[1] for j in range(z_grid+1)]
    hp_low,hp_high=entropy_interval(p)
    sqrtp_low,sqrtp_high=sqrt_interval(p)
    constant_upper=-hp_low+t*(1-sqrtp_low)
    delta_z=F(1,z_grid)
    entropy_modulus=entropy_interval(delta_z)[1]+delta_z*inner.log2[1]
    first=kappa_denominator//2
    kappas=[F(j,kappa_denominator) for j in range(first,int(maximum_kappa*kappa_denominator)+1)]
    if kappas[-1]!=maximum_kappa:kappas.append(maximum_kappa)
    point_upper=[];best=None
    for index,kappa in enumerate(kappas):
        values=[]
        for j in range(z_grid+1):
            z=F(j,z_grid)
            values.append(entropies[j]+z*inner.upper(kappa*z))
        hull=hull_at(values,p)
        r=1-p*kappa/t
        logterm_upper=log_interval(1-r*r)[1]/4
        upper=constant_upper-p*kappa+logterm_upper+hull+entropy_modulus+2*kappa*delta_z
        point_upper.append(upper)
        if best is None or upper>best[0]:best=(upper,kappa)
        if index%40==0:
            print(json.dumps({'completed_kappa_points':index+1,'total':len(kappas),
                'largest_point_upper_so_far':float(best[0])}),flush=True)
    maximum_gap=max(b-a for a,b in zip(kappas,kappas[1:]))
    # On kappa>=1/2 the exact outer expression is p-Lipschitz for p>=1/2.
    outer_upper=max(point_upper)+p*maximum_gap/2
    discr_low,discr_high=sqrt_interval(1+16*t*t)
    rho_low=(discr_low-1)/(4*t);rho_high=(discr_high-1)/(4*t)
    gaussian_upper=p*inner.log2[1]+t*(1-sqrtp_low)-t*(1-rho_high)+log_interval(1-rho_low*rho_low)[1]/4
    final=max(outer_upper,gaussian_upper)
    assert final<0
    return {'target':'ternary finite-channel latent LOWER envelope at p31/32,t4',
       'p':str(p),'t':str(t),'z_grid':z_grid,'spin_grid':spin_grid,
       'kappa_denominator':kappa_denominator,'kappa_points':len(kappas),
       'maximum_spin_bracket_width':str(inner.maximum_bracket_width),
       'largest_grid_upper_rational':str(best[0]),'largest_grid_upper':float(best[0]),
       'largest_grid_upper_kappa':str(best[1]),
       'outer_mesh_upper_rational':str(outer_upper),'outer_mesh_upper':float(outer_upper),
       'gaussian_upper_rational':str(gaussian_upper),'gaussian_upper':float(gaussian_upper),
       'certified_envelope_upper_rational':str(final),'certified_envelope_upper':float(final),
       'strictly_negative':True,
       'qualification':'exact rational upper certificate for the latent obstruction class ONLY; no Bellman upper or original cap follows'}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--z-grid',type=int,default=1000)
    parser.add_argument('--kappa-denominator',type=int,default=100)
    parser.add_argument('--spin-grid',type=int,default=100)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    result=certificate(args.z_grid,args.kappa_denominator,args.spin_grid)
    print(json.dumps(result,indent=2),flush=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
