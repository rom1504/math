#!/usr/bin/env python3
"""Numerical audit of the compatible homotopy B_u=sgn(B)(I+uQ)^(1/2).

This is a falsification tool, not a proof.  It optimizes the finite-channel
paired response by repeated conditional-score updates.
"""

from __future__ import annotations

import itertools
import math
import numpy as np
from scipy.special import ndtri
from scipy.stats import qmc


def sobol_normals(power: int, dim: int, seed: int) -> np.ndarray:
    uu = qmc.Sobol(dim, scramble=True, seed=seed).random_base2(power)
    eps = np.finfo(float).eps
    return ndtri(np.clip(uu, eps, 1 - eps))


def matrix_sqrt_psd(mat: np.ndarray) -> np.ndarray:
    ee, vv = np.linalg.eigh(mat)
    return (vv * np.sqrt(np.maximum(ee, 0))) @ vv.T


def matrix_sign(mat: np.ndarray) -> np.ndarray:
    ee, vv = np.linalg.eigh(mat)
    signs = np.where(ee >= 0, 1.0, -1.0)
    return (vv * signs) @ vv.T


def response_moments(
    yy: np.ndarray,
    ss: np.ndarray,
    gg: np.ndarray,
    rr: np.ndarray,
    ww: np.ndarray,
    inv_c: np.ndarray,
    inv_d: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    nn = len(yy)
    aa = yy.T @ ss / nn
    hh = (yy.T @ gg / nn) @ inv_c
    mm = yy.T @ rr / nn
    jj = (yy.T @ ww / nn) @ inv_d
    return aa, hh, mm, jj


def value_and_score(
    tmat: np.ndarray,
    yy: np.ndarray,
    ss: np.ndarray,
    gg: np.ndarray,
    rr: np.ndarray,
    ww: np.ndarray,
    inv_c: np.ndarray,
    inv_d: np.ndarray,
) -> tuple[float, np.ndarray]:
    aa, hh, mm, jj = response_moments(
        yy, ss, gg, rr, ww, inv_c, inv_d
    )
    size = len(tmat)
    val = (2.0 / size) * (
        np.trace(tmat @ aa @ tmat @ hh.T)
        + np.trace(tmat @ mm @ tmat @ jj.T)
    )
    # Row-sample form of
    # T H T S + T A T C^{-1}G + T J T R + T M T D^{-1}W.
    score = (
        ss @ tmat @ hh.T @ tmat
        + gg @ inv_c @ tmat @ aa.T @ tmat
        + rr @ tmat @ jj.T @ tmat
        + ww @ inv_d @ tmat @ mm.T @ tmat
    )
    return float(val), score


def optimize_at_tmat(
    tmat: np.ndarray,
    threshold: float = 0.8414699114,
    power: int = 13,
    seed: int = 17,
    steps: int = 30,
    restarts: int = 3,
) -> float:
    size = len(tmat)
    cc = tmat @ tmat
    inv_c = np.linalg.inv(cc)
    states = np.array(list(itertools.product((-1.0, 1.0), repeat=size)))
    zz = sobol_normals(power, 2 * size, seed)
    base_count = len(zz)
    sqrt_c = matrix_sqrt_psd(cc)
    gg0 = zz[:, :size] @ sqrt_c
    ss = np.repeat(states, base_count, axis=0)
    gg = np.tile(gg0, (len(states), 1))
    ff = np.sign(gg + threshold * ss)
    ff[ff == 0] = 1
    nn = len(ss)
    ls = ff.T @ ss / nn
    lg = (ff.T @ gg / nn) @ inv_c
    rr = ff - ss @ ls.T - gg @ lg.T
    kk = rr.T @ rr / nn
    dd = tmat @ kk @ tmat
    inv_d = np.linalg.inv(dd)
    sqrt_d = matrix_sqrt_psd(dd)
    ww0 = zz[:, size:] @ sqrt_d
    ww = np.tile(ww0, (len(states), 1))

    scalar = (0.5859761744, 0.6179560304, 0.2396817825, 0.4661704739)
    best = -1e100
    rng = np.random.default_rng(seed + 1000)
    for restart in range(restarts):
        if restart == 0:
            p, q, r, d = scalar
            yy = np.sign(p * ss + q * gg + r * rr + d * ww)
        else:
            coeff = rng.normal(size=4)
            yy = np.sign(
                coeff[0] * ss
                + coeff[1] * gg
                + coeff[2] * rr
                + coeff[3] * ww
            )
        yy[yy == 0] = 1
        old = -1e100
        for _ in range(steps):
            val, score = value_and_score(
                tmat, yy, ss, gg, rr, ww, inv_c, inv_d
            )
            best = max(best, val)
            new_y = np.sign(score)
            new_y[new_y == 0] = 1
            # A simultaneous best response can two-cycle.  Retain the better
            # endpoint and use a deterministic half update if needed.
            new_val, _ = value_and_score(
                tmat, new_y, ss, gg, rr, ww, inv_c, inv_d
            )
            if new_val + 1e-10 < val:
                mask = rng.random(len(yy)) < 0.5
                new_y[~mask] = yy[~mask]
            yy = new_y
            if abs(val - old) < 1e-10:
                break
            old = val
    return best


def homotopy_from_fibre(
    fibre: np.ndarray, grid: tuple[float, ...], **kwargs
) -> list[float]:
    size = len(fibre)
    t1 = fibre / math.sqrt(size)
    c1 = t1 @ t1
    oo = matrix_sign(t1)
    out = []
    for uu in grid:
        cu = np.eye(size) + uu * (c1 - np.eye(size))
        tu = oo @ matrix_sqrt_psd(cu)
        out.append(optimize_at_tmat(tu, seed=17, **kwargs))
    return out


def symmetric_fibres(size: int):
    coords = [(i, j) for i in range(size) for j in range(i, size)]
    for bits in itertools.product((-1.0, 1.0), repeat=len(coords)):
        aa = np.empty((size, size))
        for (i, j), val in zip(coords, bits):
            aa[i, j] = aa[j, i] = val
        if abs(np.linalg.det(aa)) > 1e-8:
            yield aa


if __name__ == "__main__":
    grid = (0.0, 0.25, 0.5, 0.75, 1.0)
    for size in (2, 3):
        worst_drop = (0.0, None, None)
        for idx, fibre in enumerate(symmetric_fibres(size)):
            vals = homotopy_from_fibre(
                fibre, grid, power=11, steps=12, restarts=1
            )
            drops = np.diff(vals)
            md = float(np.min(drops))
            if md < worst_drop[0]:
                worst_drop = (md, fibre.copy(), vals)
        print("size", size, "worst_drop", worst_drop[0])
        print(worst_drop[1])
        print(worst_drop[2])
