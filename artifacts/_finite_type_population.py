#!/usr/bin/env python3
"""Population state evolution for the explicit three-fibre obstruction.

The base conference dimension tends to infinity while the fibre size
stays three.  This script evaluates the resulting finite-dimensional
paired cavity functional by Sobol quadrature.
"""

from __future__ import annotations

import math
import numpy as np
from scipy.special import ndtri
from scipy.stats import qmc


T = (np.ones((3, 3)) - 2 * np.eye(3)) / math.sqrt(3)
CG = T @ T
INV_CG = np.linalg.inv(CG)
S_STATES = np.array(
    [[1 if (mask >> j) & 1 else -1 for j in range(3)] for mask in range(8)],
    dtype=float,
)


def sobol_normals(power: int, dim: int, seed: int) -> np.ndarray:
    u = qmc.Sobol(dim, scramble=True, seed=seed).random_base2(power)
    eps = np.finfo(float).eps
    return ndtri(np.clip(u, eps, 1 - eps))


def first_stage(t: float, z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return stacked S,G,R and residual covariance K."""
    chol_g = np.linalg.cholesky(CG)
    g0 = z @ chol_g.T
    ss = np.repeat(S_STATES, len(z), axis=0)
    gg = np.tile(g0, (len(S_STATES), 1))
    ff = np.sign(gg + t * ss)
    ff[ff == 0] = 1
    count = len(ss)
    ls = (ff.T @ ss) / count
    eg = (ff.T @ gg) / count
    lg = eg @ INV_CG
    rr = ff - ss @ ls.T - gg @ lg.T
    kk = (rr.T @ rr) / count
    return ss, gg, rr, kk


def population_value(
    t: float,
    coeffs: tuple[float, float, float, float],
    power: int = 16,
    seed: int = 7,
) -> tuple[float, dict[str, np.ndarray]]:
    z = sobol_normals(power, 6, seed)
    ss, gg, rr, kk = first_stage(t, z[:, :3])
    cw = T @ kk @ T
    chol_w = np.linalg.cholesky(cw)
    w0 = z[:, 3:] @ chol_w.T
    ww = np.tile(w0, (len(S_STATES), 1))
    p, q, r, d = coeffs
    yy = np.sign(p * ss + q * gg + r * rr + d * ww)
    yy[yy == 0] = 1
    count = len(ss)

    aa = (yy.T @ ss) / count
    hh = ((yy.T @ gg) / count) @ INV_CG
    cc = ((yy.T @ rr) / count) @ np.linalg.inv(kk)
    jj = ((yy.T @ ww) / count) @ np.linalg.inv(cw)
    channel1 = np.trace(T @ aa @ T @ hh.T)
    channel2 = np.trace(T @ cc @ kk @ T @ jj.T)
    value = (2 / 3) * (channel1 + channel2)
    return value, {
        "K": kk,
        "A": aa,
        "H": hh,
        "C": cc,
        "J": jj,
        "channel1": np.array(channel1),
        "channel2": np.array(channel2),
    }


if __name__ == "__main__":
    conference_coeffs = (
        0.5859761744,
        0.6179560304,
        0.2396817825,
        0.4661704739,
    )
    value, data = population_value(0.8414699114, conference_coeffs, power=18)
    print(f"value={value:.12f}")
    for key, item in data.items():
        print(key)
        print(item)
