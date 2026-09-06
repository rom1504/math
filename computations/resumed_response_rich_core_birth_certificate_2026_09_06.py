#!/usr/bin/env python3
"""Exact rational rectangle certificate for one rich-core response birth.

The policy file contributes ONLY rational rectangle endpoints. All inverse
coefficients, variances, Gaussian moments and Jensen terms are reconstructed
with outward rational intervals. The omitted conditional Z-Hermite levels
are charged by the Gaussian-regression attenuation bound.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from math import comb, factorial
import json
from pathlib import Path

from fresh_finite_anchor_fixed_point_certificate import CHILDREN
from fresh_limit_rooted_lower_certificate import I, INV_SQRT_2PI
from fresh_limit_mask_ascent_certificate import cdf_point, phi_large
from resumed_response_full_center_certificate_2026_09_06 import canonical_pair
from resumed_response_conditional_v_certificate_2026_09_06 import conditional_coefficients, hermite_endpoint
from resumed_response_feedforward_inverse_certificate_2026_09_06 import central_second_moment, to_power_coefficients, gaussian_power_moments


def power(x, k):
    result = I(1)
    for _ in range(k):
        result *= x
    return result


def magnitude(x):
    return max(abs(x.lo), abs(x.hi))


def hermite_powers(degree):
    answer = [I(0) for _ in range(degree+1)]
    for j in range(degree//2+1):
        answer[degree-2*j] = I(factorial(degree)).sqrt() * F((-1)**j, 2**j*factorial(j)*factorial(degree-2*j))
    return answer


def endpoint(x, degree):
    values = [I(1)]
    if degree:
        values.append(I(x))
    for k in range(1, degree):
        values.append(x/I(k+1).sqrt()*values[-1]-I(F(k,k+1)).sqrt()*values[-2])
    return cdf_point(x), phi_large(x), values


def moments_between(left, right, degree, cache):
    if (left,degree) not in cache:
        cache[left,degree] = endpoint(left, degree)
    if (right,degree) not in cache:
        cache[right,degree] = endpoint(right, degree)
    a, b = cache[left,degree], cache[right,degree]
    values = [b[0]-a[0]]
    for k in range(1, degree+1):
        values.append((a[1]*a[2][k-1]-b[1]*b[2][k-1])/I(k).sqrt())
    return values


def positive_resolvent_moment(j, k, R2, s2, a):
    # Integral r^(a-1) (R2+s2*r)^j (1-r)^k dr. Every summand is positive.
    total = I(0)
    for ell in range(j+1):
        beta = F(factorial(k))
        for h in range(k+1):
            beta /= a+ell+h
        term = F(comb(j,ell))*R2**(j-ell)*s2**ell*beta
        total += I(term)
    assert total.lo > 0
    return total


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--policy", type=Path, default=Path("computations/results/resumed_response_rich_core_rectangle_birth_policy_2026_09_06.json"))
    parser.add_argument("--output", type=Path)
    parser.add_argument("--u-levels", type=int, default=20)
    parser.add_argument("--theta", type=F, default=F(437,1000))
    parser.add_argument("--target", type=F, default=F(43204,100000))
    args = parser.parse_args()
    K = args.u_levels
    assert 2 <= K <= 40 and K % 2 == 0
    raw_policy = json.loads(args.policy.read_text())["first"]["rectangle_policy"]
    policy = [{key:F(row[key]) for key in ["left","right","lower_Z","upper_Z"]} for row in raw_policy]
    assert policy[0]["left"] == 0 and policy[-1]["right"] == 2
    assert all(row["left"] < row["right"] and -8 <= row["lower_Z"] <= row["upper_Z"] <= 8 for row in policy)
    assert all(a["right"] == b["left"] for a,b in zip(policy,policy[1:]))
    alpha, core_alpha, resolvent = F(91,125), F(3623,5000), F(3479,1000)
    assert all(row["right"] <= alpha or row["left"] >= alpha for row in policy)
    theta, rr = args.theta, F(9,250)
    assert 0 <= theta <= 1
    polynomial = [F(9969843,10**7),-F(728676,10**7),F(5715628,10**7),-F(14557984,10**7)]
    _, core_density, core_mass, core_covariance, _, derivative = canonical_pair(core_alpha,resolvent)
    projection, reconstructed, norm_stripped = conditional_coefficients(core_alpha,resolvent,core_density,core_mass)
    assert derivative.hi < 1
    check = reconstructed-core_covariance
    assert check.lo <= 0 <= check.hi
    degree = len(projection)-1
    norm = 2*core_density*norm_stripped.sqrt()
    rho = [F(z,10000) for z in [8108,-3110,1441,1662,-691,-1088,-758,-874,-444,331,638,496,572,502,357,563,392,452,230,286,330]]
    R2 = sum((r*r for r in rho),F(0))
    s2, s = 1-R2, I(1-R2).sqrt()
    density, mass = phi_large(alpha), 2*cdf_point(alpha)-1
    ep = hermite_endpoint(alpha,degree)
    Hcoef = [I(0) for _ in range(degree+1)]
    Hcoef[0] = mass
    for d in range(2,degree+1,2):
        Hcoef[d] = -2*density*ep[2][d-1]/I(d).sqrt()
    mean = sum((projection[d]*Hcoef[d] for d in range(0,degree+1,2)),I(0))
    rho_m = sum((x*x for x in projection),I(0))
    center_second, check_mean, check_mass = central_second_moment(projection,alpha)
    for error in [check_mean-mean,check_mass-mass]:
        assert error.lo <= 0 <= error.hi and magnitude(error) < F(1,10**22)
    v1, v2 = center_second-mean*mean/mass, rho_m-center_second
    x1,x2 = F(229,20000),-F(647,50000)
    x0 = I(1-x1*x1-x2*x2).sqrt()
    tail_scale = x2/v2.sqrt()
    center_linear = x1/v1.sqrt()-tail_scale
    center_constant = x0/mass.sqrt()-x1*mean/(mass*v1.sqrt())
    c = x0*mean/mass.sqrt()+x1*v1.sqrt()+x2*v2.sqrt()
    d = (1-c*c).sqrt()
    lam = 2*density*(2*cdf_point(c*alpha/d)-1)
    gam = 4*INV_SQRT_2PI*(1-cdf_point(alpha/d))
    old_norm = lam*lam+2*c*lam*gam+gam*gam
    old_variance = 1-mass-old_norm
    assert old_variance.lo > 0
    powers_m = to_power_coefficients(projection)
    gm = gaussian_power_moments(alpha,degree+12)
    mp = sum((powers_m[j]*polynomial[k]*gm[(j+2*k)//2]
              for j in range(0,degree+1,2) for k in range(4)),I(0))
    pp = sum((polynomial[j]*polynomial[k]*gm[j+k] for j in range(4) for k in range(4)),I(0))
    Hp = sum((polynomial[j]*gm[j] for j in range(4)),I(0))
    A = mp-rr*rho_m
    norm_f = pp-2*rr*mp+rr*rr*rho_m
    nu2 = norm_f-A*A
    assert nu2.lo > 0
    nu = nu2.sqrt()
    uf = (tail_scale*(mp-rr*rho_m)+center_linear*(mp-rr*center_second)
          +center_constant*(Hp-rr*mean))
    core_endpoint = hermite_endpoint(core_alpha,degree)
    beta = [I(0) for _ in range(degree+1)]
    beta[0] = core_mass
    for d in range(2,degree+1,2):
        beta[d] = -2*core_density*core_endpoint[2][d-1]/I(d).sqrt()
    t_anchor, anchor_scales = [], []
    for j,children in enumerate(CHILDREN):
        d = len(children)
        mult = Counter(children)
        denom, monomial = 1,F(1)
        for child,count in mult.items():
            denom *= factorial(count)
            monomial *= rho[child]**count
        vj = I(F(factorial(d),denom)).sqrt()*monomial
        hp = hermite_powers(d)
        integral = sum((hp[a]*polynomial[b]*gm[(a+2*b)//2]
                        for a in range(0,d+1,2) for b in range(4)),I(0))
        t_anchor.append(vj*(integral-rr*projection[d]))
        anchor_scales.append(rho[j]-s*beta[d]*vj/(resolvent*norm))
    t_star = (A-sum((rho[j]*t_anchor[j] for j in range(len(rho))),I(0)))/s
    b_core = [(t_anchor[j]-A*rho[j])/nu for j in range(len(rho))]
    b_core.append((t_star-A*s)/nu)
    b_norm = sum((z*z for z in b_core),I(0))
    orthogonality = sum((rho[j]*b_core[j] for j in range(len(rho))),I(0))+s*b_core[-1]
    assert orthogonality.lo <= 0 <= orthogonality.hi
    assert b_norm.hi < F(31,100)
    print(json.dumps({"stage":"exact_inverse_moments","A":A.json(),"nu2":nu2.json(),"b_norm_squared":b_norm.json()}),flush=True)

    coeff = [[I(0) for _ in range(K)] for _ in range(degree+1)]
    b_powers = [power(-s*b_core[-1],k) for k in range(K)]
    for d in range(0,degree+1,2):
        for k in range(min(d+1,K)):
            integral = positive_resolvent_moment(d-k,k,R2,s2,resolvent)*b_powers[k]
            coeff[d-k][k] += s/norm*beta[d]*I(comb(d,k)).sqrt()*integral
    for j,children in enumerate(CHILDREN):
        mult = Counter(children)
        d = len(children)
        pol = [I(1)]
        denominator = 1
        for child,count in mult.items():
            denominator *= factorial(count)
            for _ in range(count):
                following = [I(0) for _ in range(len(pol)+1)]
                for ell,z in enumerate(pol):
                    following[ell] += rho[child]*z
                    following[ell+1] += b_core[child]*z
                pol = following
        for k in range(min(d+1,K)):
            coeff[d-k][k] += anchor_scales[j]*pol[k]*I(F(factorial(d-k)*factorial(k),denominator)).sqrt()
    column_error = I(0)
    for j in range(degree+1):
        error = coeff[j][0]-projection[j]
        assert error.lo <= 0 <= error.hi and magnitude(error)<F(1,10**20)
        column_error += I(0,magnitude(error))
    print(json.dumps({"stage":"exact_bivariate_coefficients","first_column_total_error":column_error.json()}),flush=True)
    cache = {}
    records = []
    support, av, az = 2*(1-cdf_point(F(2))), 2*phi_large(F(2)), I(0)
    for index,row in enumerate(policy):
        left,right,lower,upper = [row[key] for key in ["left","right","lower_Z","upper_Z"]]
        vm = [2*z for z in moments_between(left,right,degree,cache)]
        zm = moments_between(lower,upper,K-1,cache)
        q = 1-zm[0]
        signmean = 1-cdf_point(upper)-cdf_point(lower)
        support += vm[0]*q
        av += vm[1]*signmean
        az += vm[0]*(phi_large(lower)+phi_large(upper))
        mean_m = sum((vm[j]*projection[j] for j in range(0,degree+1,2)),I(0))
        mean_g = I(0)
        for k in range(K):
            conditional = sum((vm[j]*coeff[j][k] for j in range(k%2,degree-k+1,2)),I(0))
            mean_g += conditional*zm[k]
        oldH = int(right<=alpha)
        # P(v)=P0+P1*v²+P2*v⁴+P3*v⁶ in normalized Hermites.
        poly_h = [I(polynomial[0]+polynomial[1]+3*polynomial[2]+15*polynomial[3]),
                  I(2).sqrt()*(polynomial[1]+6*polynomial[2]+45*polynomial[3]),
                  I(24).sqrt()*(polynomial[2]+15*polynomial[3]),
                  I(720).sqrt()*polynomial[3]]
        mean_poly = sum((poly_h[j]*vm[2*j] for j in range(4)),I(0))
        mean_f = -rr*mean_m+oldH*mean_poly
        mean_u = tail_scale*mean_m+oldH*(center_linear*mean_m+center_constant*vm[0])
        records.append((row,vm[0],zm[0],oldH,mean_m,mean_g,mean_f,mean_u))
        if index % 64 == 63:
            print(json.dumps({"stage":"rectangle_moments","completed":index+1}),flush=True)
    new_a,new_b = av-az*A/nu,az/nu
    norm_new = av*av+az*az
    cross = lam*new_a+lam*new_b*A+gam*new_a*c+gam*new_b*uf
    support_theta = (1-theta)*(1-mass)+theta*support
    norm_theta = (1-theta)**2*old_norm+2*theta*(1-theta)*cross+theta*theta*norm_new
    variance = support_theta-norm_theta
    assert variance.lo > 0
    tau = variance.sqrt()
    g_scale = (1-theta)*lam+theta*new_a
    total = I(0)
    report_bins = []
    for row,vmass,zmass,oldH,mean_m,mean_g,mean_f,mean_u in records:
        factor = (1-theta)*oldH+theta*zmass
        bin_mass = vmass*factor
        if bin_mass.lo <= F(1,10**18):
            # The true contribution is nonnegative. Avoid dividing by a tiny
            # interval mass whose enclosure could straddle zero.
            report_bins.append({**{k:str(v) for k,v in row.items()},"discarded_nonnegative_bin":True})
            continue
        bin_g = (1-theta)*oldH*mean_m+theta*mean_g
        bin_scalar = factor*((1-theta)*gam*mean_u+theta*new_b*mean_f)
        bin_K = g_scale*bin_g+bin_scalar
        argument = bin_K/(bin_mass*tau)
        assert -8 < argument.lo <= argument.hi < 8
        value = bin_K*(2*cdf_point(argument)-1)+2*bin_mass*tau*phi_large(argument)
        total += value
        report_bins.append({**{k:str(v) for k,v in row.items()},"mass":bin_mass.json(),"integral_K_truncated":bin_K.json(),"Jensen_term":value.json()})
    truncation_penalty = I(magnitude(g_scale))*theta*power(b_norm,K//2)
    lower = total-truncation_penalty
    assert lower.lo > args.target
    result = {"method":"exact_fraction_rich_core_rectangle_response_birth",
              "verified":True,"target":str(args.target),"theta":str(theta),
              "surrogate_P_even_power_coefficients":[str(x) for x in polynomial],
              "surrogate_tail_m_coefficient":str(-rr),"retained_Z_Hermite_levels":K,
              "core_derivative_energy":derivative.json(),"inverse_A":A.json(),
              "inverse_nu_squared":nu2.json(),"core_Z_covariance_norm_squared":b_norm.json(),
              "new_support":support.json(),"new_aV":av.json(),"new_aZ":az.json(),
              "old_new_inverse_inner_product":cross.json(),"line_residual_variance":variance.json(),
              "line_g_coefficient":g_scale.json(),"bin_sum_before_truncation":total.json(),
              "conditional_projection_tail_penalty":truncation_penalty.json(),
              "lower_bound":lower.json(),"margin_above_target":(lower-args.target).json(),
              "policy_source":str(args.policy),"bin_certificates":report_bins}
    rendered = json.dumps(result,indent=2)+"\n"
    if args.output:
        args.output.write_text(rendered)
        print(json.dumps({k:v for k,v in result.items() if k!="bin_certificates"},indent=2))
    else:
        print(rendered,end="")


if __name__ == "__main__":
    main()
