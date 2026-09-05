"""Deterministic numerical exploration of a two-field Gaussian mask ascent.

This is NOT an interval certificate.  It only computes Gaussian integrals
for masks covered by the proved hierarchical Gaussianization theorem.
"""

import json
import math
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
from scipy.special import ndtr
from numpy.polynomial.legendre import leggauss


ALPHA = 37 / 50
RESOLVENT = 97 / 10
DEGREE = 200
PHI0 = 1 / math.sqrt(2 * math.pi)


def phi(x):
    return PHI0 * np.exp(-np.asarray(x) ** 2 / 2)


def folded(mean, sd):
    if sd < 1e-15:
        return np.abs(mean)
    z = mean / sd
    return 2 * sd * phi(z) + mean * (2 * ndtr(z) - 1)


def tail_abs(alpha, variance, covariance):
    rem = math.sqrt(max(0, variance - covariance * covariance))
    if rem < 1e-14:
        return 2 * abs(covariance) * float(phi(alpha))
    return (2 * covariance * float(phi(alpha))
            * (2 * ndtr(covariance * alpha / rem) - 1)
            + 4 * math.sqrt(variance) * PHI0
            * ndtr(-alpha * math.sqrt(variance) / rem))


class Ascent:
    def __init__(self, inner_order=96, threshold_scale=1, threshold_shift=0):
        self.inner_order = inner_order
        self.threshold_scale = threshold_scale
        self.threshold_shift = threshold_shift
        he = np.zeros(DEGREE + 1)
        he[0], he[1] = 1, ALPHA
        for r in range(1, DEGREE):
            he[r + 1] = (ALPHA * he[r] - math.sqrt(r) * he[r - 1]) / math.sqrt(r + 1)
        self.mass = 2 * ndtr(ALPHA) - 1
        c = np.zeros(DEGREE + 1)
        c[0] = self.mass
        for r in range(2, DEGREE + 1, 2):
            c[r] = -2 * float(phi(ALPHA)) * he[r - 1] / math.sqrt(r)
        u = c / (RESOLVENT + np.arange(DEGREE + 1))
        self.u = u / np.linalg.norm(u)
        self.w = float(np.dot(self.u, c))
        self.s = math.sqrt(self.mass - self.w ** 2)
        self.a = 2 * float(phi(ALPHA)) * (2 * ndtr(self.w * ALPHA / self.s) - 1)
        self.b = 4 * PHI0 / math.sqrt(self.mass) * ndtr(-ALPHA * math.sqrt(self.mass) / self.s)
        self.nodes, self.weights = leggauss(inner_order)
        self.p1 = self.integrate(lambda v: self.conditional_mass(v))
        self.r1 = self.integrate(lambda v: self.g(v) * self.conditional_mass(v))
        self.r0 = 2 * quad(lambda v: float(phi(v)) * self.conditional_mass(v),
                           0, ALPHA, epsabs=1e-11, epsrel=1e-11, limit=200)[0]
        inv_cov = np.linalg.inv([[1, self.w], [self.w, self.mass]])
        self.regression = inv_cov @ np.array([self.r1, self.r0])
        self.residual = self.p1 - np.dot([self.r1, self.r0], self.regression)
        assert self.residual > 0

    def g(self, v):
        prev, current = 1.0, float(v)
        result = self.u[0]
        for r in range(1, DEGREE):
            nxt = (v * current - math.sqrt(r) * prev) / math.sqrt(r + 1)
            result += self.u[r + 1] * nxt
            prev, current = current, nxt
        return result

    def threshold(self, v):
        return max(0, self.threshold_scale * self.a * self.g(v)
                   + self.b * (v <= ALPHA) + self.threshold_shift)

    def conditional_mass(self, v):
        threshold = self.threshold(v)
        return ndtr((threshold - self.w * v) / self.s) - ndtr((-threshold - self.w * v) / self.s)

    def conditional_signed_mass(self, v):
        threshold = self.threshold(v)
        mean = self.w * v
        return (ndtr((threshold - mean) / self.s)
                + ndtr((-threshold - mean) / self.s) - 2 * ndtr(-mean / self.s))

    def conditional_absolute_moment(self, v):
        threshold = self.threshold(v)
        mean = self.w * v
        return (mean * self.conditional_signed_mass(v)
                + self.s * (2 * phi(mean / self.s)
                            - phi((threshold - mean) / self.s)
                            - phi((threshold + mean) / self.s)))

    def minorant(self):
        baseline = tail_abs(ALPHA, self.mass, self.w)
        old_absolute = 2 * PHI0 * math.sqrt(self.mass) - baseline
        new_absolute = self.integrate(self.conditional_absolute_moment)
        absolute_difference = new_absolute - old_absolute
        old_signed = 2 * PHI0 * self.w / math.sqrt(self.mass) - self.a - self.b * self.w
        new_signed = self.integrate(lambda v: v * self.conditional_signed_mass(v))
        signed_difference = new_signed - old_signed
        gamma = (self.a * (self.r1 - self.w) + self.b * (self.r0 - self.mass)
                 - absolute_difference)
        mass_difference = self.mass + self.p1 - 2 * self.r0
        c, d = self.regression
        quadratic = c * signed_difference + (d - 1) * absolute_difference
        epsilon = min(1, max(0, gamma / (2 * quadratic))) if quadratic > 0 else 1
        return {
            "gradient": gamma, "mask_L2_difference_squared": mass_difference,
            "exact_quadratic_loss": quadratic,
            "old_absolute_moment": old_absolute,
            "new_absolute_moment": new_absolute,
            "old_V_signed_moment": old_signed,
            "new_V_signed_moment": new_signed,
            "minorant_epsilon": epsilon,
            "quadratic_minorant_value": baseline + epsilon * gamma - epsilon ** 2 * quadratic,
            "coarse_minorant_value": baseline + gamma ** 2 / (4 * mass_difference),
        }

    @staticmethod
    def integrate(fn):
        # Numerical truncation only; this script does not claim a tail certificate.
        return 2 * sum(quad(lambda v: float(phi(v)) * fn(v), lo, hi,
                            epsabs=1e-10, epsrel=1e-10, limit=400)[0]
                       for lo, hi in [(0, ALPHA), (ALPHA, 2), (2, 4), (4, 8)])

    def masked_folded(self, v, epsilon):
        threshold = self.threshold(v)
        mean = self.w * v
        lo = max(-10, (-threshold - mean) / self.s)
        hi = min(10, (threshold - mean) / self.s)
        if hi <= lo:
            return 0.0
        z = lo + (self.nodes + 1) * (hi - lo) / 2
        w_values = mean + self.s * z
        c, d = self.regression
        new_mean = epsilon * c * v + (1 - epsilon + epsilon * d) * w_values
        new_sd = epsilon * math.sqrt(self.residual)
        return float((hi - lo) / 2 * np.dot(self.weights, phi(z) * folded(new_mean, new_sd)))

    def objective(self, epsilon):
        variance = ((1 - epsilon) ** 2 * self.mass
                    + 2 * epsilon * (1 - epsilon) * self.r0
                    + epsilon ** 2 * self.p1)
        covariance = (1 - epsilon) * self.w + epsilon * self.r1
        outside_old = tail_abs(ALPHA, variance, covariance)
        inside_new = self.integrate(lambda v: self.masked_folded(v, epsilon))
        return ((1 - epsilon) * outside_old
                + epsilon * (2 * PHI0 * math.sqrt(variance) - inside_new))

    def report(self):
        solution = minimize_scalar(lambda e: -self.objective(e), bounds=(0.001, 1),
                                   method="bounded", options={"xatol": 1e-8})
        return {
            "status": "numerical_only_not_a_certificate",
            "alpha": ALPHA, "resolvent": RESOLVENT, "degree": DEGREE,
            "inner_legendre_order": self.inner_order,
            "threshold_scale": self.threshold_scale,
            "threshold_shift": self.threshold_shift,
            "original_mask_mass": self.mass, "original_covariance": self.w,
            "gradient_a": self.a, "gradient_b": self.b,
            "new_mask_mass": self.p1, "covariance_V_new": self.r1,
            "covariance_W_new": self.r0,
            "conditional_regression": self.regression.tolist(),
            "conditional_residual_variance": self.residual,
            "epsilon": float(solution.x), "new_objective": float(-solution.fun),
            "hard_update_objective": self.objective(1),
            "baseline": self.objective(0),
            "quadratic_minorant": self.minorant(),
            "grid": [{"epsilon": e, "objective": self.objective(e)}
                     for e in [0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1]],
        }


if __name__ == "__main__":
    print(json.dumps(Ascent().report(), indent=2))
