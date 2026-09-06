#!/usr/bin/env python3
"""Floating diagnostic for an actual rich-core response birth.

The new inverse is a scalar surrogate f(V), but the old inverse g is NOT
replaced by E[g|V]. Its exact bivariate polynomial projection onto (V,Z)
is integrated against the new center before conditional Jensen.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import json
from math import comb, factorial
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import ndtr, roots_jacobi

from fresh_finite_anchor_fixed_point_certificate import CHILDREN
from resumed_response_scalar_full_optimizer_2026_09_06 import hermites, phi, PHI0


def psi(k, variance):
    tau = np.sqrt(variance)
    return 2 * tau * phi(k / tau) + k * (2 * ndtr(k / tau) - 1)


class RichCore:
    def __init__(self, source):
        self.mcoef = np.array([float(z["lower"]) for z in source["conditional_projection_coefficients"]])
        self.degree = len(self.mcoef) - 1
        self.rho = np.array([8108, -3110, 1441, 1662, -691, -1088,
                             -758, -874, -444, 331, 638, 496, 572, 502,
                             357, 563, 392, 452, 230, 286, 330]) / 10000
        self.R2 = float(self.rho @ self.rho)
        self.innovation = np.sqrt(1 - self.R2)
        self.resolvent = 3.479
        self.alpha_core = .7246
        self.norm = 2 * phi(self.alpha_core) * np.sqrt(float(source["projection_norm_stripped"]["lower"]))
        self.beta = np.zeros(self.degree + 1)
        self.beta[0] = 2 * ndtr(self.alpha_core) - 1
        endpoint = hermites(np.array([self.alpha_core]), self.degree)[0]
        for d in range(2, self.degree + 1, 2):
            self.beta[d] = -2 * phi(self.alpha_core) * endpoint[d - 1] / np.sqrt(d)
        self.anchor_projection = []
        self.anchor_scale = []
        for j, children in enumerate(CHILDREN):
            mult = Counter(children)
            denominator = np.prod([factorial(k) for k in mult.values()], dtype=float)
            v = np.sqrt(factorial(len(children)) / denominator)
            for child, count in mult.items():
                v *= self.rho[child] ** count
            self.anchor_projection.append(v)
            self.anchor_scale.append(self.rho[j] - self.innovation * self.beta[len(children)] * v / (self.resolvent * self.norm))
        self.anchor_projection = np.array(self.anchor_projection)
        self.anchor_scale = np.array(self.anchor_scale)

    def conditional_polynomial(self, b, resolvent_nodes=104):
        """C[j,k] is the coefficient of h_j(V) h_k(Z) in E[g|V,Z]."""
        assert abs(np.r_[self.rho, self.innovation] @ b) < 2e-10
        assert b @ b <= 1 + 2e-10
        degree = self.degree
        coeff = np.zeros((degree + 1, degree + 1))
        x, weights = roots_jacobi(resolvent_nodes, 0, self.resolvent - 1)
        r = (x + 1) / 2
        weights = weights / 2 ** self.resolvent
        ar = self.R2 + self.innovation ** 2 * r
        br = self.innovation * (r - 1) * b[-1]
        ap = ar[None, :] ** np.arange(degree + 1)[:, None]
        bp = br[None, :] ** np.arange(degree + 1)[:, None]
        for d in range(0, degree + 1, 2):
            for k in range(d + 1):
                moment = weights @ (ap[d - k] * bp[k])
                coeff[d - k, k] += self.innovation / self.norm * self.beta[d] * np.sqrt(float(comb(d, k))) * moment
        for j, children in enumerate(CHILDREN):
            mult = Counter(children)
            d = len(children)
            pol = np.array([1.0])
            denominator = 1
            for child, count in mult.items():
                denominator *= factorial(count)
                for _ in range(count):
                    pol = np.polynomial.polynomial.polymul(pol, [self.rho[child], b[child]])
            for k in range(d + 1):
                normalized = np.sqrt(factorial(d - k) * factorial(k) / denominator)
                coeff[d - k, k] += self.anchor_scale[j] * pol[k] * normalized
        error = np.max(np.abs(coeff[:, 0] - self.mcoef))
        assert error < 5e-12, error
        assert np.sum(coeff * coeff) <= 1 + 1e-10
        return coeff, error


def run(source, nodes=256, cutoff=12, surrogate_scale=1.0, simple_surrogate=False, polynomial_surrogate=False, rectangle_bins=0, return_context=False):
    core = RichCore(source)
    degree = core.degree
    alpha = .728
    x1, x2 = .01145, -.01294
    x0 = np.sqrt(1 - x1*x1 - x2*x2)
    z, wg = np.polynomial.legendre.leggauss(nodes)
    # Positive-half symmetry; split at the sole old mask discontinuity.
    vc = alpha * (z + 1) / 2
    vt = alpha + (cutoff - alpha) * (z + 1) / 2
    v = np.r_[vc, vt]
    weights = np.r_[alpha * wg * phi(vc), (cutoff - alpha) * wg * phi(vt)]
    H = np.r_[np.ones(nodes), np.zeros(nodes)]
    hv = hermites(v, degree)
    m = hv @ core.mcoef
    p = 2 * ndtr(alpha) - 1
    center_mean = float(weights @ (H * m))
    center_second = float(weights @ (H * m * m))
    norm_m2 = float(core.mcoef @ core.mcoef)
    var1 = center_second - center_mean**2 / p
    var2 = norm_m2 - center_second
    # u = tail_scale*m + H*(central_linear*m + central_constant).
    tail_scale = x2 / np.sqrt(var2)
    central_linear = x1 / np.sqrt(var1) - tail_scale
    central_constant = x0 / np.sqrt(p) - x1 * center_mean / (p * np.sqrt(var1))
    u = tail_scale * m + H * (central_linear * m + central_constant)
    c = x0 * center_mean / np.sqrt(p) + x1 * np.sqrt(var1) + x2 * np.sqrt(var2)
    d = np.sqrt(1 - c*c)
    lam = 2 * phi(alpha) * (2 * ndtr(c * alpha / d) - 1)
    gam = 4 * PHI0 * ndtr(-alpha / d)
    t = 1 - p - lam*lam - 2*c*lam*gam - gam*gam
    Kbar = lam * m + gam * u
    B = float(weights @ (H * phi(Kbar / np.sqrt(t)))) / np.sqrt(t)
    # The parameter modifies only a legal old-measurable surrogate.
    f_tail_scale = -2 * B * (lam + gam * tail_scale)
    f_correction = H * (surrogate_scale * (2 * ndtr(Kbar / np.sqrt(t)) - 1)
                         -2 * B * gam * (central_linear * m + central_constant))
    surrogate_polynomial = None
    if polynomial_surrogate:
        target_values = f_correction[:nodes] + (f_tail_scale + .036) * m[:nodes]
        surrogate_polynomial = np.array([.9969843, -.0728676, .5715628, -1.4557984])
        f_tail_scale = -.036
        f_correction = H * np.polynomial.polynomial.polyval(v*v, surrogate_polynomial)
    elif simple_surrogate:
        f_tail_scale = -.036
        f_correction = H
    f = f_tail_scale * m + f_correction
    mf = f_tail_scale * norm_m2 + float(weights @ (m * f_correction))
    normf2 = (f_tail_scale*f_tail_scale * norm_m2
              +2*f_tail_scale*float(weights @ (m*f_correction))
              +float(weights @ (f_correction*f_correction)))
    nu = np.sqrt(normf2 - mf*mf)
    # Covariances of L=Uf with every independent coordinate of the rich core.
    core_l_cov = []
    for j, children in enumerate(CHILDREN):
        degree_j = len(children)
        integral = f_tail_scale * core.mcoef[degree_j] + float(weights @ (f_correction * hv[:, degree_j]))
        core_l_cov.append(core.anchor_projection[j] * integral)
    core_l_cov = np.array(core_l_cov)
    star_cov = (mf - core.rho @ core_l_cov) / core.innovation
    cov = np.r_[core_l_cov, star_cov]
    b = (cov - mf * np.r_[core.rho, core.innovation]) / nu
    coeff, projection_error = core.conditional_polynomial(b)
    conditional_coefficients = hv @ coeff
    # New tested signing is an actual function of V and L=mf*V+nu*Z.
    threshold = np.maximum(psi(Kbar, t) - B, 0)
    mean = mf * v
    lower = (-threshold - mean) / nu
    upper = (threshold - mean) / nu
    qnew = ndtr(-upper) + ndtr(lower)
    fnew_mean = ndtr(-upper) - ndtr(lower)
    new_p = float(weights @ qnew)
    new_aV = float(weights @ (v * fnew_mean))
    new_aZ = float(weights @ (phi(upper) + phi(lower)))
    new_a = new_aV - new_aZ * mf / nu
    new_fscale = new_aZ / nu
    mass_new_H = 1 - qnew
    # Exact truncated-Hermite endpoint identity, evaluated in floating point.
    # Clip only where Gaussian tails are far below this diagnostic's precision.
    lower_clip, upper_clip = np.clip(lower, -40, 40), np.clip(upper, -40, 40)
    hl, hu = hermites(lower_clip, degree - 1), hermites(upper_clip, degree - 1)
    moments = np.zeros_like(conditional_coefficients)
    moments[:, 0] = mass_new_H
    for k in range(1, degree + 1):
        moments[:, k] = (phi(lower_clip) * hl[:, k-1] - phi(upper_clip) * hu[:, k-1]) / np.sqrt(k)
    weighted_g = np.sum(conditional_coefficients * moments, axis=1)
    normold = lam*lam + 2*c*lam*gam + gam*gam
    normnew = new_aV*new_aV + new_aZ*new_aZ
    # <u,f> known exactly from m norm and compact-center correction.
    uf = (tail_scale*f_tail_scale*norm_m2
          + float(weights @ (H*(central_linear*m+central_constant)*f_tail_scale*m))
          + float(weights @ (u*f_correction)))
    cross = lam*new_a + lam*new_fscale*mf + gam*new_a*c + gam*new_fscale*uf

    rectangle_records = None
    if rectangle_bins:
        assert polynomial_surrogate
        rational_alpha = Fraction(91, 125)
        endpoints = ([rational_alpha * j / rectangle_bins for j in range(rectangle_bins + 1)]
                     + [rational_alpha + (2-rational_alpha) * j / rectangle_bins for j in range(1, rectangle_bins + 1)])
        left = np.array([float(x) for x in endpoints[:-1]])
        right = np.array([float(x) for x in endpoints[1:]])
        midpoint = (left+right)/2
        mid_m = hermites(midpoint, degree) @ core.mcoef
        old_h = np.r_[np.ones(rectangle_bins), np.zeros(rectangle_bins)]
        mid_u = tail_scale*mid_m + old_h*(central_linear*mid_m+central_constant)
        mid_k = lam*mid_m + gam*mid_u
        thresholds = np.maximum(psi(mid_k, t)-B, 0)
        lower_int = np.rint(np.clip((-thresholds-mf*midpoint)/nu, -8, 8)*10000).astype(int)
        upper_int = np.rint(np.clip((thresholds-mf*midpoint)/nu, -8, 8)*10000).astype(int)
        lo, hi = lower_int/10000, upper_int/10000
        mass_v = 2*(ndtr(right)-ndtr(left))
        moment_v = 2*(phi(left)-phi(right))
        mass_u = ndtr(hi)-ndtr(lo)
        qnew_rect = 1-mass_u
        new_p = float(mass_v @ qnew_rect + 2*ndtr(-2))
        new_aV = float(moment_v @ (ndtr(-hi)-ndtr(lo)) + 2*phi(2))
        new_aZ = float(mass_v @ (phi(lo)+phi(hi)))
        new_a = new_aV-new_aZ*mf/nu
        new_fscale = new_aZ/nu
        normnew = new_aV*new_aV+new_aZ*new_aZ
        cross = lam*new_a+lam*new_fscale*mf+gam*new_a*c+gam*new_fscale*uf
        integral_v = np.zeros((len(left), degree+1))
        integral_u = np.zeros_like(integral_v)
        integral_v[:, 0] = mass_v
        integral_u[:, 0] = mass_u
        hv_left, hv_right = hermites(left, degree-1), hermites(right, degree-1)
        hu_left, hu_right = hermites(lo, degree-1), hermites(hi, degree-1)
        for k in range(1, degree+1):
            integral_v[:, k] = 2*(phi(left)*hv_left[:, k-1]-phi(right)*hv_right[:, k-1])/np.sqrt(k)
            integral_u[:, k] = (phi(lo)*hu_left[:, k-1]-phi(hi)*hu_right[:, k-1])/np.sqrt(k)
        bin_m = integral_v @ core.mcoef
        bin_g = np.sum((integral_v @ coeff)*integral_u, axis=1)
        # Low power polynomial converted by exact Hermite identities.
        poly_hermite = np.array([surrogate_polynomial[0]+surrogate_polynomial[1]+3*surrogate_polynomial[2]+15*surrogate_polynomial[3],
                                np.sqrt(2)*(surrogate_polynomial[1]+6*surrogate_polynomial[2]+45*surrogate_polynomial[3]),
                                np.sqrt(24)*(surrogate_polynomial[2]+15*surrogate_polynomial[3]),
                                np.sqrt(720)*surrogate_polynomial[3]])
        bin_poly = integral_v[:, :7:2] @ poly_hermite
        bin_f = -.036*bin_m+old_h*bin_poly
        bin_u = tail_scale*bin_m+old_h*(central_linear*bin_m+central_constant*mass_v)
        rectangle_records = [{"left":str(endpoints[j]), "right":str(endpoints[j+1]),
                              "lower_Z":str(Fraction(int(lower_int[j]),10000)),
                              "upper_Z":str(Fraction(int(upper_int[j]),10000))}
                             for j in range(len(left))]

    def line(theta, wrong_independence=False):
        support = (1-theta)*(1-p) + theta*new_p
        norm = (1-theta)**2*normold + 2*theta*(1-theta)*cross + theta*theta*normnew
        variance = support - norm
        if rectangle_bins:
            factor = (1-theta)*old_h+theta*mass_u
            rectangle_mass = factor*mass_v
            rectangle_g = (1-theta)*old_h*bin_m+theta*bin_g
            rectangle_scalar = factor*((1-theta)*gam*bin_u+theta*new_fscale*bin_f)
            rectangle_K = ((1-theta)*lam+theta*new_a)*rectangle_g+rectangle_scalar
            select = rectangle_mass > 1e-15
            return float(np.sum(rectangle_mass[select]*psi(rectangle_K[select]/rectangle_mass[select], variance)))
        mass = (1-theta)*H + theta*mass_new_H
        integral_g = (1-theta)*H*m + theta*weighted_g
        if wrong_independence:
            integral_g = mass*m
        scale_g = (1-theta)*lam + theta*new_a
        scalar = (1-theta)*gam*u + theta*new_fscale*f
        integral_K = scale_g*integral_g + mass*scalar
        select = mass > 1e-15
        return float(weights[select] @ (mass[select] * psi(integral_K[select]/mass[select], variance)))

    opt = minimize_scalar(lambda theta: -line(theta), bounds=(0, 1), method="bounded",
                          options={"xatol": 2e-12})
    candidates = [(0., line(0.)), (1., line(1.)), (float(opt.x), -float(opt.fun))]
    theta_best, value_best = max(candidates, key=lambda pair: pair[1])
    if return_context:
        return {"core":core,"conditional_coefficients":coeff,"alpha":alpha,"mass":p,
                "lambda":lam,"gamma":gam,"c":c,"old_variance":t,"A":mf,"nu":nu,
                "uf":uf,"tail_scale":tail_scale,"center_linear":central_linear,
                "center_constant":central_constant,"polynomial":surrogate_polynomial,
                "b_norm_squared":float(b@b)}
    return {
        "status": "floating_diagnostic_not_certificate", "nodes_each_interval": nodes,
        "cutoff": cutoff, "surrogate_scale": surrogate_scale,
        "simple_surrogate": simple_surrogate,
        "surrogate_even_power_coefficients": None if surrogate_polynomial is None else surrogate_polynomial.tolist(),
        "rectangle_policy": rectangle_records,
        "fixed_rational_theta_value": line(.437),
        "old_conditional_value": line(0), "best_theta": theta_best,
        "best_value": value_best, "gain": value_best-line(0),
        "full_step": line(1), "wrong_independence_at_best": line(theta_best, True),
        "line_derivative": (line(1e-5)-line(0))/1e-5,
        "line_values": {str(theta): line(theta) for theta in np.linspace(0, 1, 11)},
        "lambda": lam, "gamma": gam, "old_variance": t, "B_surrogate": B,
        "birth_covariance_V": mf, "birth_innovation_sd": nu,
        "birth_core_covariances": cov.tolist(), "innovation_core_covariances": b.tolist(),
        "innovation_projection_norm_squared": float(b @ b),
        "conditional_polynomial_norm_squared": float(np.sum(coeff*coeff)),
        "conditional_V_coefficient_replay_error": float(projection_error),
        "conditional_U_degree_norms_squared": np.sum(coeff*coeff, axis=0).tolist(),
        "new_support": new_p, "new_aV": new_aV, "new_aZ": new_aZ,
        "new_g_coefficient": new_a, "new_f_coefficient": new_fscale,
        "old_new_first_chaos_inner_product": cross,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--projection", type=Path, default=Path("computations/results/resumed_response_conditional_v_certificate_2026_09_06.json"))
    parser.add_argument("--output", type=Path)
    parser.add_argument("--nodes", type=int, default=256)
    parser.add_argument("--cutoff", type=float, default=12.)
    parser.add_argument("--surrogate-scale", type=float, default=1.)
    parser.add_argument("--simple-surrogate", action="store_true")
    parser.add_argument("--polynomial-surrogate", action="store_true")
    parser.add_argument("--rectangle-bins", type=int, default=0)
    parser.add_argument("--single", action="store_true")
    args = parser.parse_args()
    source = json.loads(args.projection.read_text())
    first = run(source, args.nodes, args.cutoff, args.surrogate_scale, args.simple_surrogate, args.polynomial_surrogate, args.rectangle_bins)
    result = {"first": first}
    if not args.single:
        replay = run(source, 2*args.nodes, args.cutoff, args.surrogate_scale, args.simple_surrogate, args.polynomial_surrogate, args.rectangle_bins)
        result["double_node_replay"] = replay
        result["best_value_difference"] = replay["best_value"] - first["best_value"]
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+"\n")
    for record in result.values():
        if isinstance(record, dict):
            print(json.dumps({key: record[key] for key in ["nodes_each_interval", "old_conditional_value", "best_theta", "best_value", "gain", "full_step", "line_derivative", "wrong_independence_at_best", "conditional_polynomial_norm_squared", "innovation_projection_norm_squared"]}), flush=True)


if __name__ == "__main__":
    main()
