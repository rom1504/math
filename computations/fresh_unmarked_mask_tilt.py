"""Numerical diagnostic, not an interval certificate: unmarked mask tilt.

Uses the proved conference/all-power AMP coupling. The base mask is exactly
the compact-support two-field mask in fresh_limit_gaussian_mask_ascent.md.
"""

import json
import math
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
from scipy.special import ndtr
from numpy.polynomial.legendre import leggauss

from fresh_limit_fixed_point_mask_ascent import Ascent, ALPHA, PHI0, phi


class Tilt(Ascent):
    @staticmethod
    def integrate(fn):
        return 2 * sum(quad(lambda v: float(phi(v)) * fn(v), lo, hi,
                            epsabs=2e-10, epsrel=2e-10, limit=200)[0]
                       for lo, hi in [(0, ALPHA), (ALPHA, 1), (1, 2),
                                      (2, 4), (4, 8)])

    def threshold(self, v):
        if abs(v) > 1:
            return 0.0
        return max(0.0, 2.5 * self.a * self.g(v)
                   + self.b * (abs(v) <= ALPHA) - 0.8)

    def __init__(self, inner_order=128):
        super().__init__(inner_order=inner_order)
        self.m_v = (2 * PHI0 * self.w / math.sqrt(self.mass)
                    - self.integrate(lambda v: v * self.conditional_signed_mass(v)))
        self.m_w = (2 * PHI0 * math.sqrt(self.mass)
                    - self.integrate(self.conditional_absolute_moment))
        inv = np.linalg.inv([[1, self.w], [self.w, self.mass]])
        self.linear = inv @ np.array([self.m_v, self.m_w])
        self.residual_f = 1 - self.p1 - np.dot(self.linear, [self.m_v, self.m_w])
        assert self.residual_f > 0
        self.baseline_minorant = float(np.dot(self.regression, [self.m_v, self.m_w]))

    def conditional(self, v, tau):
        """Return E[M erf], E[W M erf], E[M sign(W) erf], E[M phi]."""
        mean, sd, threshold = self.w * v, self.s, self.threshold(v)
        # Integrate the two unmasked W tails, splitting at zero if needed.
        cut = [mean - 10 * sd, -threshold, 0.0, threshold, mean + 10 * sd]
        cut = sorted(set(x for x in cut if mean - 10*sd <= x <= mean + 10*sd))
        total = np.zeros(4)
        for lo, hi in zip(cut, cut[1:]):
            if hi <= lo or abs((lo + hi) / 2) < threshold:
                continue
            values = lo + (self.nodes + 1) * (hi - lo) / 2
            density = phi((values - mean) / sd) / sd
            smooth = 2 * ndtr(values / tau) - 1
            components = np.array([smooth, values*smooth,
                                   np.sign(values)*smooth, phi(values/tau)])
            total += (hi-lo)/2 * (components @ (self.weights*density))
        return total

    def objective(self, tau):
        first = self.integrate(lambda v: v * self.conditional(v, tau)[0])
        second = self.integrate(lambda v: self.conditional(v, tau)[1])
        absolute = self.integrate(lambda v: self.conditional(v, tau)[2])
        density = self.integrate(lambda v: self.conditional(v, tau)[3])
        feature_correlation = (absolute - self.linear[0]*first
                               - self.linear[1]*second) / math.sqrt(self.residual_f)
        unmarked_derivative = 2 * density
        paired = float(np.dot(self.regression, [first, second]))
        extra = unmarked_derivative * feature_correlation
        return dict(value=paired+extra, paired=paired, extra=extra,
                    derivative=unmarked_derivative,
                    feature_correlation=feature_correlation)

    def report_tilt(self):
        sol = minimize_scalar(lambda t: -self.objective(t)["value"],
                              bounds=(0.002, 0.5), method="bounded",
                              options={"xatol": 1e-6})
        return dict(status="numerical_only_not_an_interval_certificate",
                    inner_order=self.inner_order,
                    p=self.mass, w=self.w, s=self.s,
                    new_mask_mass=self.p1,
                    regression=self.regression.tolist(),
                    response_linear_projection=self.linear.tolist(),
                    response_residual_variance=self.residual_f,
                    baseline=self.baseline_minorant,
                    optimized_tau=float(sol.x), result=self.objective(sol.x),
                    grid=[dict(tau=t, **self.objective(t))
                          for t in [0.01, 0.03, 0.05, 0.1, 0.15, 0.2]])


if __name__ == "__main__":
    import sys
    print(json.dumps(Tilt(int(sys.argv[1]) if len(sys.argv)>1 else 128).report_tilt(),
                     indent=2))
