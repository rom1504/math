"""Directed interval upper bound for the actual ternary E construction.

Every accepted rectangle is verified; floating point never prunes a box.
The optional mean-value bound uses the exact envelope derivative and the
interval range of its clipped maximizing mixture weight.
"""
import argparse
from fractions import Fraction as F
import json
import math
from pathlib import Path
import time
import mpmath as mp


def raw_fraction(raw):
    sign, mantissa, exponent, _ = raw
    answer = F(-mantissa if sign else mantissa)
    return answer * (2**exponent) if exponent >= 0 else answer / (2**(-exponent))


def endpoints(interval):
    return tuple(map(raw_fraction, interval._mpi_))


def rational_iv(value):
    value = F(value)
    return mp.iv.mpf(value.numerator) / value.denominator


def range_iv(lower, upper):
    return mp.iv.mpf([rational_iv(lower).a, rational_iv(upper).b])


class Certificate:
    def __init__(self, p, t, target, precision=50, monotone_only=False):
        mp.iv.dps = precision
        self.p, self.t, self.target = p, t, target
        self.pp, self.tt = rational_iv(p), rational_iv(t)
        self.monotone_only = monotone_only

    def optimized_gain(self, A, B):
        # Gain is increasing in B and decreasing in A. Outward endpoints
        # therefore give a rigorous bound even when a branch boundary is near.
        aa = endpoints(A)[0]
        bb = endpoints(B)[1]
        A, B, p = rational_iv(aa), rational_iv(bb), self.pp
        if bb == 0 or self.p*bb <= aa:
            return mp.iv.mpf(0)
        if aa == 0 or self.p*bb-aa >= (1-self.p)*aa*bb:
            return p*mp.iv.ln(1+B)-mp.iv.ln(1+A)
        return (p*mp.iv.ln(p)+(1-p)*mp.iv.ln(1-p)-p*mp.iv.ln(A)
                +mp.iv.ln(B)-(1-p)*mp.iv.ln(B-A))

    def gain_arguments(self, lam, z):
        A = mp.iv.exp(lam*z*z/self.pp)-1
        h = lam*z/self.pp
        # 2 sinh(h)^2, equivalent to cosh(2h)-1, keeps B nonnegative.
        B = (mp.iv.exp(h)-mp.iv.exp(-h))**2/2
        return A, B

    def point_upper(self, lam, z):
        ll, zz = rational_iv(lam), rational_iv(z)
        A, B = self.gain_arguments(ll, zz)
        cminus = mp.iv.ln(ll*(2*self.tt-ll)/(self.tt*self.tt))/4-ll
        return cminus+self.optimized_gain(A, B)

    def monotone_upper(self, rect):
        l,u,a,b = rect
        ll, uu, aa, bb = map(rational_iv, rect)
        A = mp.iv.exp(ll*aa*aa/self.pp)-1
        h = uu*bb/self.pp
        B = (mp.iv.exp(h)-mp.iv.exp(-h))**2/2
        cminus = mp.iv.ln(ll*(2*self.tt-ll)/(self.tt*self.tt))/4-ll
        return endpoints(cminus+self.optimized_gain(A,B))[1]

    def meanvalue_upper(self, rect):
        l,u,a,b = rect
        if a == 0:
            return None
        lam, z = range_iv(l,u), range_iv(a,b)
        p, t = self.pp, self.tt
        A, B = self.gain_arguments(lam,z)
        if endpoints(A)[0] <= 0 or endpoints(B)[0] <= 0:
            return None
        unclipped = (p*B-A)/((1-p)*A*B)
        vlo,vhi = endpoints(unclipped)
        vlo,vhi = max(F(0),min(F(1),vlo)), max(F(0),min(F(1),vhi))
        v = range_iv(vlo,vhi)
        expA = mp.iv.exp(lam*z*z/p)
        q = 2*lam*z/p
        sinhq = (mp.iv.exp(q)-mp.iv.exp(-q))/2
        dl = (1/lam-1/(2*t-lam))/4-1
        dl += 2*v*z*sinhq/(1+v*B)-v*z*z*expA/(p*(1+v*A))
        dz = 2*v*lam*sinhq/(1+v*B)-2*v*lam*z*expA/(p*(1+v*A))
        dlmax = max(map(abs,endpoints(dl)))
        dzmax = max(map(abs,endpoints(dz)))
        center = endpoints(self.point_upper((l+u)/2,(a+b)/2))[1]
        return center+dlmax*(u-l)/2+dzmax*(b-a)/2

    def run(self, max_boxes):
        stack = [(F(1,2),self.t,F(0),F(1))]
        checked=accepted=meanvalue_accepted=0
        worst=None
        start=time.monotonic()
        while stack:
            rect=stack.pop()
            checked+=1
            upper=self.monotone_upper(rect)
            use_meanvalue=False
            if upper>=self.target and not self.monotone_only:
                meanvalue=self.meanvalue_upper(rect)
                if meanvalue is not None and meanvalue<upper:
                    upper=meanvalue
                    use_meanvalue=True
            if upper<self.target:
                accepted+=1
                meanvalue_accepted+=use_meanvalue
                worst=upper if worst is None else max(worst,upper)
            else:
                l,u,a,b=rect
                if (u-l)>=self.t*(b-a):
                    mid=(l+u)/2
                    stack.extend([(l,mid,a,b),(mid,u,a,b)])
                else:
                    mid=(a+b)/2
                    stack.extend([(l,u,a,mid),(l,u,mid,b)])
            if checked%10000==0:
                print(json.dumps({'checked':checked,'accepted':accepted,'pending':len(stack),
                                  'seconds':time.monotonic()-start}),flush=True)
            if checked>=max_boxes:
                return {'status':'budget stop; NOT CERTIFIED','checked':checked,
                        'accepted':accepted,'pending':len(stack)}
        rho=(mp.iv.sqrt(1+16*self.tt*self.tt)-1)/(4*self.tt)
        gaussian=-self.tt*(1-rho)+mp.iv.ln(1-rho*rho)/4
        assert endpoints(gaussian)[1]<self.target
        cap=(self.tt+self.pp*mp.iv.ln(2)+rational_iv(self.target))/(2*self.tt*mp.iv.sqrt(self.pp))
        return {'status':'directed interval certificate','p':str(self.p),'t':str(self.t),
                'E_upper':str(self.target),'cap_interval':list(map(str,endpoints(cap))),
                'cap_float_upper_DISPLAY_ONLY':float(endpoints(cap)[1]),
                'gaussian_interval':list(map(str,endpoints(gaussian))),
                'checked':checked,'accepted':accepted,'meanvalue_accepted':meanvalue_accepted,
                'monotone_only':self.monotone_only,
                'worst_accepted':str(worst),'seconds':time.monotonic()-start}


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--p',default='24/25')
    parser.add_argument('--t',default='97/20')
    parser.add_argument('--target',default='-10302/12500')
    parser.add_argument('--max-boxes',type=int,default=2000000)
    parser.add_argument('--precision',type=int,default=50)
    parser.add_argument('--monotone-only',action='store_true')
    parser.add_argument('--output')
    args=parser.parse_args()
    result=Certificate(F(args.p),F(args.t),F(args.target),args.precision,args.monotone_only).run(args.max_boxes)
    print(json.dumps(result,indent=2),flush=True)
    if args.output:
        Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
