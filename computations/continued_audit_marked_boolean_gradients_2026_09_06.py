"""Complete Boolean check of marked gradients and matrix Poincare."""

import json

import numpy as np

from continued_director_steiner_signing_2026_09_06 import make_signing


def check(a):
    n = len(a)
    count = 2**n
    labels = np.arange(count)
    seeds = 1 - 2 * ((labels[:, None] >> np.arange(n)) & 1)
    b = a / np.sqrt(n - 1)
    q = b @ b
    g = seeds @ b
    marked = seeds * (g**2 - 1) / np.sqrt(2)
    y, v = marked @ b, marked @ q
    w = g + y + v
    cube = w**3
    gradient_gram = np.zeros((n, n))
    worst_marked, worst_cube = 0.0, 0.0
    for s_index, s in enumerate(seeds):
        jd = (np.diag((g[s_index]**2 - 1) / np.sqrt(2))
              + np.sqrt(2) * (s * g[s_index])[:, None] * b
              - np.sqrt(2) * s[:, None] * b**2 * s[None, :])
        flip = s_index ^ (1 << np.arange(n))
        actual_jd = (marked[s_index, :, None] - marked[flip].T) / (2 * s[None, :])
        jw = b + (b + q) @ jd
        jc = (3 * w[s_index, :, None]**2 * jw
              - 6 * w[s_index, :, None] * jw**2 * s[None, :]
              + 4 * jw**3)
        actual_jc = (cube[s_index, :, None] - cube[flip].T) / (2 * s[None, :])
        worst_marked = max(worst_marked, float(np.max(np.abs(jd - actual_jd))))
        worst_cube = max(worst_cube, float(np.max(np.abs(jc - actual_jc))))
        gradient_gram += jc @ jc.T / count
    covariance = (cube - np.mean(cube, axis=0)).T @ (cube - np.mean(cube, axis=0)) / count
    min_poincare = float(np.linalg.eigvalsh(gradient_gram - covariance)[0])
    exact_marked_cov = (1 - 3 / (n - 1)) * np.eye(n) + 2 / (n - 1) * q
    cov_error = float(np.max(np.abs(marked.T @ marked / count - exact_marked_cov)))
    assert worst_marked < 1e-11
    assert worst_cube < 1e-9
    assert min_poincare > -1e-8
    assert cov_error < 1e-11
    return {"order": n, "all_seed_count": count,
            "marked_gradient_error": worst_marked,
            "composite_cube_gradient_error": worst_cube,
            "marked_covariance_error": cov_error,
            "minimum_matrix_Poincare_difference_eigenvalue": min_poincare}


if __name__ == "__main__":
    steiner, _ = make_signing(2)
    rng = np.random.default_rng(9060920)
    upper = np.triu(rng.choice((-1, 1), size=(8, 8)), 1)
    print(json.dumps({"status": "complete finite-seed floating identity checks",
                      "records": [check(steiner), check(upper + upper.T)]}, indent=2))
