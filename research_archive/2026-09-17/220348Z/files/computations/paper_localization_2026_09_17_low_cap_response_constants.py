"""Replay of explicit margin-to-response constants; no finite-n CLT claim."""

import json
import math
from pathlib import Path


def main():
    kappa=math.sqrt(2/math.pi)
    b0=0.5-math.pi/12
    c1=math.sqrt(8/(kappa**5*b0*b0))
    c2=(373248/(25*kappa**5))**(1/6)
    rows=[]
    checks=0
    for groth,rstar,cap,energy in ((2,0.429,0.5,0.30),
                                 (2,0.43,0.5,0.30),
                                 (math.pi/(2*math.asinh(1)),0.429,0.5,0.433),
                                 (math.pi/(2*math.asinh(1)),0.433,0.494,0.433)):
        gap=2*energy-4*groth*(cap-rstar)
        assert gap>0
        ell=c2+(5/math.sqrt(2))*math.sqrt(c1)+6*groth*(cap+rstar)*c1**(1/3)
        d0=min(1,(8*c1)**(-2),(gap/(2*ell))**6)
        assert c1*math.sqrt(d0)<=1/8+1e-15
        assert ell*d0**(1/6)<=gap/2+1e-15
        for j in range(1,301):
            d=10**(-1-j/8)
            a=min(c1*math.sqrt(d),1/8)
            b=min(c2*c2*d**(1/3),9)
            t=math.sqrt(b)/36
            assert t<=1/12+1e-15
            signed_variance=kappa*kappa*(t*b-216*t**3)
            assert abs(signed_variance-(5*kappa*kappa/216)*b**1.5)<1e-13
            assert abs(kappa*signed_variance**2/8-25*kappa**5*b**3/373248)<1e-13
            exact=(math.sqrt(b)+(5/math.sqrt(2))*math.sqrt(a)
                   +4*groth*rstar*(1-(1-a**(1/3))**1.5)
                   +groth*cap*(4*a**(2/3)+2*a**(1/3)))
            upper=ell*d**(1/6)
            assert exact<=upper+1e-12
            checks+=1
        rows.append(dict(K=groth,half_range=rstar,cap=cap,energy=energy,
                         gap=gap,L=ell,discount=d0/2))
    floor=kappa*(math.sqrt(4/3)+math.sqrt(2/3))/2
    assert floor>0.75
    result=dict(status='PASS',algebra_checks=checks,C1=c1,C2=c2,
                rows=rows,paired_half_band_floor=floor,
                scope='Algebraic constants and exponent replay; asymptotic source theorem is separately audited.')
    path=Path('tmp/paper_portfolio_2026_09_17/localization/low_cap_response_constants.json')
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
