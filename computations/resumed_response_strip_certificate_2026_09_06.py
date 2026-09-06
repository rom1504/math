"""Exact rational certificate for a high-value creation mask with zero-strip slack.

The proof is in artifacts/resumed_response_high_value_zero_strip_2026_09_06.md.
No optimizer, floating point arithmetic, or unenclosed quadrature is used.
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import product
from math import comb, factorial
import json
import argparse
from pathlib import Path

from fresh_limit_rooted_lower_certificate import I, phi, Phi, INV_SQRT_2PI
from fresh_limit_mask_ascent_certificate import phi_large, cdf_point
from fresh_finite_anchor_fixed_point_certificate import CHILDREN


def square(x):
    x=I.get(x)
    return I(0,max(x.lo*x.lo,x.hi*x.hi)) if x.lo<=0<=x.hi else I(min(x.lo*x.lo,x.hi*x.hi),max(x.lo*x.lo,x.hi*x.hi))


def inverse_negative_part():
    alpha,a,degree=F(361,500),F(17,5),200
    rho=[F(t,10000) for t in [8108,-3110,1441,1662,-691,-1088,-758,-874,-444,331,638,496,572,502,357,563,392,452,230,286,330]]
    total=sum((r*r for r in rho),F(0)); innovation=1-total
    density=phi(alpha); mass=2*Phi(alpha)-1
    he=[F(1),alpha]
    for d in range(1,degree):he.append(alpha*he[-1]-d*he[-2])
    grouped=[F(0) for _ in range(degree+1)]
    for d in range(2,degree+1,2):
        raw=he[d-1]**2/factorial(d)
        for ell in range(d+1):grouped[ell]+=raw*comb(d,ell)*total**(d-ell)*innovation**ell
    ac=[mass];removed=F(0)
    for ch in CHILDREN[1:]:
        mult=Counter(ch);denom=1;monomial=F(1)
        for j,k in mult.items():denom*=factorial(k);monomial*=rho[j]**k
        raw=he[len(ch)-1]*monomial
        removed+=raw*raw/denom
        ac.append(-2*density*raw/I(denom).sqrt())
    grouped[0]-=removed
    assert all(v>=0 for v in grouped)
    norm=sum((v/(a+ell)**2 for ell,v in enumerate(grouped)),F(0))
    scale=I(innovation).sqrt()/(2*density*I(norm).sqrt())
    delta=[rho[j]-scale/a*ac[j] for j in range(len(rho))]
    assert delta[0].lo>0
    indices=[]
    for ch in CHILDREN:
        ct=Counter(ch);indices.append(tuple(ct[j] for j in range(9)))
    coefficients=defaultdict(lambda:I(0))
    for aa in range(1,len(delta)):
        for bb in range(1,len(delta)):
            options=[]
            for m,n in zip(indices[aa],indices[bb]):
                entries=[]
                for k in range(min(m,n)+1):
                    j=m+n-2*k
                    cc=factorial(k)*comb(m,k)*comb(n,k)*I(F(factorial(j),factorial(m)*factorial(n))).sqrt()
                    entries.append((j,cc))
                options.append(entries)
            for items in product(*options):
                idx=tuple(v[0] for v in items)
                cc=delta[aa]*delta[bb]
                for _,v in items:cc*=v
                coefficients[idx]+=cc
    fourth=sum((square(v) for v in coefficients.values()),I(0))
    correction_negative=fourth.sqrt()/(4*delta[0])
    retained=square(mass)+4*square(density)*sum((he[d-1]**2/F(factorial(d)) for d in range(2,degree+1,2)),F(0))
    hermite_tail=mass-retained
    assert hermite_tail.lo>0
    d=degree+2
    binomial_resolvent=sum((F(comb(d,ell))*total**(d-ell)*innovation**ell/(a+ell)**2 for ell in range(d+1)),F(0))
    tail_bound=scale*(hermite_tail*binomial_resolvent).sqrt()
    negative_bound=correction_negative+tail_bound
    assert negative_bound.hi<F(3,100)
    covariance=sum((rho[j]*ac[j] for j in range(len(rho))),I(0))+I(innovation).sqrt()*2*density*sum((v/(a+ell) for ell,v in enumerate(grouped)),F(0))/I(norm).sqrt()
    assert F(5297,10000)<mass.lo<=mass.hi<F(5298,10000)
    assert F(7004,10000)<covariance.lo
    assert 2*density.hi<F(62,100)
    return {'mask_mass':mass.json(),'covariance':covariance.json(),'correction_constant':delta[0].json(),'centered_correction_fourth_moment':fourth.json(),'correction_negative_norm_bound':correction_negative.json(),'resolvent_tail_norm_bound':tail_bound.json(),'inverse_negative_norm_bound':negative_bound.json(),'fourth_moment_terms':len(coefficients)}


def strip_bounds(bins=128):
    alpha=F(361,500)
    e,r=F(3,100),F(23,1000)
    inner,outer=F(1,1000),F(3,25)
    pmin,pmax=F(5297,10000),F(5298,10000)
    cmin=F(7004,10000)-e*r
    vmax=pmax+r*r
    sd=I(1-cmin*cmin/vmax).sqrt()
    rootp=I(pmin).sqrt()
    # Monotonicity hypotheses, stated explicitly in the proof.
    assert outer<cmin*alpha
    assert 2*outer*vmax<alpha*pmin*cmin
    assert outer*outer<pmin
    assert ((I(alpha)-cmin/vmax*outer)/sd).lo**2>2
    def point(w):
        u=(I(alpha)-cmin/vmax*w)/sd
        b=(I(alpha)+cmin/vmax*w)/sd
        t0=2-cdf_point(u)-cdf_point(b)
        t2=t0+u*phi_large(u)+b*phi_large(b)
        assert t0.lo>0 and t2.lo>0
        dens=phi(I(w)/rootp)/rootp
        return dens,t0,t2
    mass_upper=I(0);cost_upper=I(0);trace_upper=I(0)
    # On each bin density decreases, t0,t2 increase, and q decreases.
    # Use endpoint products, giving literal Darboux upper sums.
    for start,stop in ((F(0),inner),(inner,outer)):
        count=1 if stop==inner else bins
        width=(stop-start)/count
        left=point(start)
        for k in range(count):
            wl=start+k*width;wr=wl+width;right=point(wr)
            qleft=F(1) if wl<=inner else (outer-wl)/(outer-inner)
            density=I(left[0].hi)
            t0=I(right[1].hi);t2=I(right[2].hi)
            mass_upper+=2*width*density*t0*qleft*qleft
            cost_upper+=2*width*density*t0*wr*qleft
            if start==inner:
                trace_upper+=2*width/(outer-inner)**2*density*(wr*wr/pmin*t0+t2)
            left=right
    assert mass_upper.hi<r*r
    assert trace_upper.hi<1
    assert cost_upper.hi<F(4,100000)
    base_lower=F(43065,100000)
    final_lower=I(base_lower-F(62,100)*e*r)-cost_upper
    assert final_lower.lo>F(43,100)
    return {'bins':bins,'invariant_radius':str(r),'negative_part_budget':str(e),'plateau_width':str(inner),'support_width':str(outer),'worst_covariance':str(cmin),'worst_variance':str(vmax),'feedback_L2_squared_upper':mass_upper.json(),'radius_squared':str(r*r),'full_direction_derivative_squared_upper':trace_upper.json(),'removed_absolute_energy_upper':cost_upper.json(),'new_mask_J_lower':final_lower.json()}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    negative=inverse_negative_part()
    strips=strip_bounds()
    rendered=json.dumps({'negative_part':negative,'strip_feedback':strips,'verified':True},indent=2)
    if args.output is not None:
        args.output.write_text(rendered+'\n',encoding='utf-8')
    print(rendered,flush=True)


if __name__=='__main__':main()
