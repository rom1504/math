"""Finite diagnostics for the proved signed two-spin cavity falsifier."""

import itertools
import json
import math

from scipy.integrate import quad
from scipy.special import ndtr


def phi(t):
    return math.exp(-t * t / 2) / math.sqrt(2 * math.pi)


def shifted_absolute(s, variance):
    scale = math.sqrt(variance)
    return 2 * scale * phi(s / scale) + s * (2 * ndtr(s / scale) - 1)


def clipped_expectation(function):
    interior = quad(lambda s: function(s) * phi(s / 10) / 10, -1, 1)[0]
    endpoint_mass = 2 * (1 - ndtr(0.1))
    return interior + endpoint_mass * (function(-1) + function(1)) / 2


def response(variance):
    return clipped_expectation(lambda s: shifted_absolute(s, variance))


grid = [-20, -3, -1, -0.5, 0, 0.5, 1, 3, 20]
for g1, g2 in itertools.product(grid, repeat=2):
    values = [y1 * y2 + g1 * y1 + g2 * y2
              for y1, y2 in itertools.product([-1, 1], repeat=2)]
    signed_formula = max(1, abs(g2)) + abs(g1 + max(-1, min(1, g2)))
    absolute_formula = 1 + abs(g1) + abs(g2)
    assert max(values) == signed_formula
    assert max(map(abs, values)) == absolute_formula

variance = 0.25
curvature = clipped_expectation(
    lambda s: phi(s / math.sqrt(variance)) * (s * s - variance)
) / (2 * variance ** 2.5)
step = 0.01
mixture_gap = (response(variance - step) + response(variance + step)) / 2 - response(variance)
assert curvature > 0.424
assert mixture_gap > 0
print(json.dumps({
    "identity_grid_pairs": len(grid) ** 2,
    "both_identities_pass": True,
    "interior_probability": 2 * ndtr(0.1) - 1,
    "signed_curvature_at_one_quarter": curvature,
    "rigorous_lower_bound": 0.424,
    "variance_mixture_step": step,
    "signed_variance_mixture_gap": mixture_gap,
    "scope": "signed branch; full absolute cap is concave"
}, indent=2))
