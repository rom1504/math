"""Independent finite checks; run with the repository-local .venv Python.

No saved report is overwritten. Exact matrix counts precede floating tests.
The stdout JSON is a reproducible diagnostic, not a numerical proof.
"""

from fractions import Fraction
import json
import math

import numpy as np
from scipy.linalg import hadamard
from scipy.special import ndtri


def twin_check(m):
    h = hadamard(m).astype(np.int64)
    c = np.kron(np.ones((2, 2), dtype=np.int64), h)
    a = c - np.diag(np.diag(c))
    n = 2 * m
    a2 = a @ a
    target = np.kron(np.ones((2, 2), dtype=np.int64), np.eye(m, dtype=np.int64))
    delta = a2 - (n - 1) * target
    identity_delta = a2 - (n - 1) * np.eye(n, dtype=np.int64)
    trace_a3 = int(np.sum(a * a2.T))
    assert trace_a3 == 0
    assert np.all(np.diag(a) == 0)
    assert np.all(a[np.triu_indices(n, 1)] ** 2 == 1)
    q = a2.astype(float) / (n - 1)
    t = math.pi / 2
    cos_factors = np.cos(t * q)
    np.fill_diagonal(cos_factors, 1.0)
    nonlinear_root_test = float(np.mean(np.sin(t * np.diag(q)) * np.prod(cos_factors, axis=1)))
    row_var = np.sum(q * q, axis=1)
    fourth = 3 * row_var**2 - 2 * np.sum(q**4, axis=1)
    return {
        "m": m,
        "n": n,
        "trace_A3_exact": trace_a3,
        "normalized_return_error_to_twin_exact": str(Fraction(int(np.sum(delta * delta)), n * (n - 1)**2)),
        "normalized_return_error_to_identity_exact": str(Fraction(int(np.sum(identity_delta * identity_delta)), n * (n - 1)**2)),
        "mean_return_variance": float(np.mean(row_var)),
        "mean_return_fourth_moment": float(np.mean(fourth)),
        "mean_S_sin_pi_return_over_2": nonlinear_root_test,
        "involution_value_for_same_test": 1.0,
        "operator_norm_upper_bound": (2 * math.sqrt(m) + 1) / math.sqrt(n - 1),
    }


def profile_check():
    u = np.linspace(-1, 1, 100001)
    phi = np.exp(-ndtri((1 + u) / 2)**2 / 2) / math.sqrt(2 * math.pi)
    lower = (1 - u * u) / math.sqrt(2 * math.pi)
    assert np.min(phi - lower) >= -1e-14
    return {"grid_size": len(u), "minimum_profile_slack": float(np.min(phi - lower)), "center_phi": float(phi[len(u) // 2])}


def one_root_mask_model():
    # (R1,R2) are independent standard normals; (X1,X2) have correlation b.
    # Own-root R/X correlations vanish; crossed ones equal a.
    # H(x)=(1+cos(tx))/2 is a bounded even [0,1] mask.
    a, b, t = 0.5, 0.5, 1.0
    sigma = np.array([[1, 0, 0, a], [0, 1, a, 0], [0, a, 1, b], [a, 0, b, 1]], dtype=float)
    assert np.min(np.linalg.eigvalsh(sigma)) > 0
    covariance = a*a*t*t*math.exp(-t*t)*math.sinh(b*t*t) / 4
    return {
        "coordinate_order": ["R1", "R2", "X1", "X2"],
        "covariance_matrix": sigma.tolist(),
        "crossed_mask_covariance_exact_formula": "a^2 t^2 exp(-t^2) sinh(b t^2)/4",
        "crossed_mask_covariance": covariance,
        "independent_block_model_mask_covariance": 0.0,
        "warning": "A Gaussian covariance-information counterexample, not an actual-signing family or an energy counterexample.",
    }


if __name__ == "__main__":
    print(json.dumps({
        "twin_signings": [twin_check(m) for m in [2, 4, 8, 16, 32, 64, 128, 256]],
        "rounding_profile_grid": profile_check(),
        "one_root_mask_model": one_root_mask_model(),
    }, indent=2))
