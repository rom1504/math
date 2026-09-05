"""Rigorous rational midpoint certificate for a two-field Gaussian mask.

Every arithmetic operation uses the audited outward Fraction interval class.
Second-derivative midpoint errors and threshold-crossing range errors are
explicit. No floating-point arithmetic or numerical integration is used.
"""

from fractions import Fraction as F
from math import factorial
import json
import sys

from fresh_limit_rooted_lower_certificate import (
    I, DIGITS, INV_SQRT_2PI, phi, Phi, exp_negative,
)


def magnitude(x):
    return max(abs(x.lo), abs(x.hi))


def phi_large(x):
    x = I.get(x)
    if x.lo <= 0 <= x.hi:
        square = I(0,max(x.lo*x.lo,x.hi*x.hi))
    else:
        square = I(min(x.lo*x.lo,x.hi*x.hi),max(x.lo*x.lo,x.hi*x.hi))
    y = exp_negative(square / 64)
    for _ in range(5):
        y = y * y
    return INV_SQRT_2PI * y


CDF_TAIL = F(1, 10 ** 24)
assert 8 * F(32) ** 129 / (factorial(129) * 259) < CDF_TAIL


def cdf_point(x):
    """Taylor interval for an essentially point input in [-8,8]."""
    x = I.get(x)
    assert -8 <= x.lo <= x.hi <= 8
    if x.hi < 0:
        return 1 - cdf_point(-x)
    if x.lo < 0:
        return I(F(1, 2)) + INV_SQRT_2PI * I(x.lo, x.hi)
    u = x * x / 2
    term = I(1)
    total = term
    for k in range(1, 129):
        term = term * (-u) / k
        total += term / (2 * k + 1)
    return I(F(1, 2)) + INV_SQRT_2PI * (x * total + I(-CDF_TAIL, CDF_TAIL))


class Jet:
    def __init__(self, value, first=0, second=0):
        self.v, self.d, self.dd = I.get(value), I.get(first), I.get(second)

    @staticmethod
    def get(x):
        return x if isinstance(x, Jet) else Jet(x)

    def __add__(self, x):
        x = Jet.get(x)
        return Jet(self.v + x.v, self.d + x.d, self.dd + x.dd)

    __radd__ = __add__

    def __neg__(self):
        return Jet(-self.v, -self.d, -self.dd)

    def __sub__(self, x):
        return self + (-Jet.get(x))

    def __rsub__(self, x):
        return Jet.get(x) - self

    def __mul__(self, x):
        x = Jet.get(x)
        return Jet(self.v*x.v,
                   self.d*x.v+self.v*x.d,
                   self.dd*x.v+2*self.d*x.d+self.v*x.dd)

    __rmul__ = __mul__

    def __truediv__(self, x):
        # All divisions in this script are by constants, not variable jets.
        assert not isinstance(x, Jet)
        return self * (1 / I.get(x))


def density_jet(x):
    p = phi_large(x.v)
    return Jet(p, -x.v*p*x.d,
               (x.v*x.v-1)*p*x.d*x.d-x.v*p*x.dd)


def cdf_jet(x, coarse=False):
    p = phi_large(x.v)
    value = I(0, 1) if coarse else cdf_point(x.v)
    return Jet(value, p*x.d, p*x.dd-x.v*p*x.d*x.d)


class Certificate:
    def __init__(self):
        self.alpha, self.resolvent = F(37, 50), F(97, 10)
        self.degree = 200
        self.sqrt_int = [I(r).sqrt() for r in range(self.degree + 1)]
        self.mass = 2*Phi(I(self.alpha))-1
        self.density = phi(I(self.alpha))
        basis = self.basis(I(self.alpha))
        raw = [I(0) for _ in range(self.degree + 1)]
        coeff = [I(0) for _ in range(self.degree + 1)]
        coeff[0] = self.mass
        raw[0] = self.mass/self.resolvent
        for r in range(2, self.degree + 1, 2):
            coeff[r] = -2*self.density*basis[r-1]/self.sqrt_int[r]
            raw[r] = coeff[r]/(self.resolvent+r)
        norm = sum((x*x for x in raw), I(0)).sqrt()
        self.u = [x/norm for x in raw]
        self.derivative_energy = sum((r*x*x for r,x in enumerate(self.u)), I(0))
        assert self.derivative_energy.hi < 1
        self.w = sum((x*y for x,y in zip(self.u,coeff)), I(0))
        self.s = (self.mass-self.w*self.w).sqrt()
        self.a = 2*self.density*(2*cdf_point(self.w*self.alpha/self.s)-1)
        self.b = 4*INV_SQRT_2PI/self.mass.sqrt()*(1-cdf_point(self.alpha*self.mass.sqrt()/self.s))
        self.baseline = self.a*self.w+self.b*self.mass
        self.old_absolute = 2*INV_SQRT_2PI*self.mass.sqrt()-self.baseline
        self.old_signed = 2*INV_SQRT_2PI*self.w/self.mass.sqrt()-self.a-self.b*self.w
        self.derivative_coeff = [self.u[r]*self.sqrt_int[r]
                                 for r in range(1,self.degree+1)]
        self.second_coeff = [self.u[r]*self.sqrt_int[r]*self.sqrt_int[r-1]
                             for r in range(2,self.degree+1)]

        # Mehler: sum_{j<=197} He_j(v)^2/j! <= 128 for |v|<=1.
        # t=99/100, t^-197<8, (1-t²)^-1/2<8, exp(t/(1+t))<2.
        t = F(99,100)
        assert t**197 > F(1,8)
        assert 1-t*t > F(1,64)
        # ||g'''||² <=199*198*D<39402.  Cauchy then gives this bound.
        self.third_bound = F(2300)
        assert 199*198*128 < self.third_bound**2

    def basis(self, x):
        result = [I(1), x]
        for r in range(1,self.degree):
            result.append((x*result[-1]-self.sqrt_int[r]*result[-2])/self.sqrt_int[r+1])
        return result

    def g_at(self, center):
        basis = self.basis(I(center))
        g = sum((x*y for x,y in zip(self.u,basis)),I(0))
        d = sum((x*y for x,y in zip(self.derivative_coeff,basis)),I(0))
        dd = sum((x*y for x,y in zip(self.second_coeff,basis)),I(0))
        return Jet(g,d,dd)

    def g_range(self, g, radius):
        delta = I(-radius,radius)
        delta_squared = I(0,radius*radius)
        rem0 = self.third_bound*radius**3/6
        rem1 = self.third_bound*radius**2/2
        rem2 = self.third_bound*radius
        return Jet(g.v+g.d*delta+g.dd*delta_squared/2+I(-rem0,rem0),
                   g.d+g.dd*delta+I(-rem1,rem1),
                   g.dd+I(-rem2,rem2))

    def integrands(self, v, g, old_mask, coarse=False):
        threshold = g*(F(5,2)*self.a)+self.b*old_mask-F(4,5)
        mean = v*self.w
        upper, lower, zero = (threshold-mean)/self.s, (-threshold-mean)/self.s, -mean/self.s
        pu, pl, pz = cdf_jet(upper,coarse), cdf_jet(lower,coarse), cdf_jet(zero,coarse)
        probability = pu-pl
        signed = pu+pl-2*pz
        absolute = mean*signed+(2*density_jet(zero)-density_jet(upper)-density_jet(lower))*self.s
        density = 2*density_jet(v)
        return [density*probability, density*g*probability,
                old_mask*density*probability, density*v*signed, density*absolute]

    def integrate(self, bins_per_unit):
        totals = [I(0) for _ in range(5)]
        crossing_bins = 0
        active_bins = 0
        # Both discontinuities (alpha and the imposed support endpoint 1)
        # are exact partition endpoints.
        for left,right,old_mask in [(F(0),self.alpha,1),(self.alpha,F(1),0)]:
            count = (int((right-left)*bins_per_unit)+1)
            width = (right-left)/count
            for j in range(count):
                lo,hi = left+j*width,left+(j+1)*width
                center,radius = (lo+hi)/2,width/2
                g0 = self.g_at(center)
                gr = self.g_range(g0,radius)
                threshold = F(5,2)*self.a*gr.v+self.b*old_mask-F(4,5)
                if threshold.hi <= 0:
                    continue
                if threshold.lo <= 0:
                    crossing_bins += 1
                    # P(|N(mean,s²)|<=T) <=2*T/(s*sqrt(2*pi)).
                    prob = 2*I(0,threshold.hi)*INV_SQRT_2PI/self.s
                    weight = 2*INV_SQRT_2PI*width
                    bounds = [prob, I(-magnitude(gr.v),magnitude(gr.v))*prob,
                              old_mask*prob, I(-hi,hi)*prob,
                              I(0,threshold.hi)*prob]
                    for k,bound in enumerate(bounds):
                        totals[k] += weight*bound
                    continue
                active_bins += 1
                central = self.integrands(Jet(center,1),g0,old_mask)
                ranged = self.integrands(Jet(I(lo,hi),1),gr,old_mask,coarse=True)
                for k in range(5):
                    error = magnitude(ranged[k].dd)*width**3/24
                    totals[k] += width*central[k].v+I(-error,error)
                if active_bins % 250 == 0:
                    print("certified active bins",active_bins,file=sys.stderr,flush=True)
        return totals,{"active_bins":active_bins,"threshold_crossing_bins":crossing_bins}

    def finish(self, integrals, metadata):
        p1,r1,r0,new_signed,new_absolute = integrals
        det = self.mass-self.w*self.w
        c = (self.mass*r1-self.w*r0)/det
        d = (r0-self.w*r1)/det
        absolute_difference = new_absolute-self.old_absolute
        signed_difference = new_signed-self.old_signed
        gamma = self.a*(r1-self.w)+self.b*(r0-self.mass)-absolute_difference
        quadratic = c*signed_difference+(d-1)*absolute_difference
        value = self.baseline+gamma-quadratic  # epsilon=1 exactly.
        target = F(429,1000)
        assert value.lo > target, value.json()
        return {
            "method":"exact_fraction_second_derivative_midpoint_intervals",
            "grid_decimal_digits":DIGITS,
            "alpha":"37/50","resolvent":"97/10","Hermite_degree":200,
            "mask":"1{|V|<=1, |W|<max(0,(5/2)*a*g(V)+b*H0-4/5)}",
            "step_epsilon":"1","uniform_g_third_derivative_bound":"2300",
            "quadrature":metadata,
            "baseline":self.baseline.json(),
            "conditional_derivative_energy":self.derivative_energy.json(),
            "new_mask_mass":p1.json(),"covariance_V_new":r1.json(),
            "covariance_W_new":r0.json(),"new_V_signed_moment":new_signed.json(),
            "new_absolute_moment":new_absolute.json(),
            "gradient":gamma.json(),"quadratic_loss":quadratic.json(),
            "lower_bound":value.json(),"target":"429/1000",
            "margin_above_target":(value-target).json(),"verified":True,
        }


def main():
    certificate = Certificate()
    totals,metadata = certificate.integrate(2000)
    metadata["nominal_bins_per_unit"] = 2000
    print(json.dumps(certificate.finish(totals,metadata),indent=2))


if __name__ == "__main__":
    main()
