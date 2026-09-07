#!/usr/bin/env python3
"""Exhaustive row-sign shore spectra for saved conference matrices."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np


ROOT = Path("/home/math/quadra")
FILES = [
    "conference_double_p5.json",
    "conference_order10_gf9.json",
    "conference_double_p13.json",
    "conference_double_p17.json",
]
KAPPA = math.pi / 2 - 1


def projector(block: np.ndarray, sign: int) -> float:
    vals = np.linalg.eigvalsh(sign * block)
    vals = vals[vals > 1e-9]
    if not len(vals):
        return 0.0
    theta = min(1.0, vals.sum() / (2 * KAPPA * len(vals)))
    return (2 / math.pi) * (theta * vals.sum() - KAPPA * theta**2 * len(vals))


def nuclear(block: np.ndarray) -> float:
    size = len(block)
    if size < 2:
        return 0.0
    return max(0.0, np.abs(np.linalg.eigvalsh(block)).sum() / math.pi - (1 - 2 / math.pi) * size)


def one_sided(block: np.ndarray, sign: int) -> int:
    size = len(block)
    if size < 2:
        return 0
    codes = np.arange(1 << (size - 1), dtype=np.uint32)[:, None]
    bits = ((codes >> np.arange(size - 1, dtype=np.uint32)) & 1).astype(np.int16)
    spins = np.concatenate((np.ones((len(codes), 1), np.int16), 1 - 2 * bits), axis=1)
    values = np.einsum("bi,ij,bj->b", spins, block, spins, optimize=True)
    return int(values.max() if sign > 0 else -values.min())


def main() -> None:
    for filename in FILES:
        matrix = np.asarray(json.loads((ROOT / "computations/results" / filename).read_text())["conference_matrix"], dtype=np.int16)
        n = len(matrix)
        codes = np.arange(1 << (n - 1), dtype=np.uint32)[:, None]
        bits = ((codes >> np.arange(n - 1, dtype=np.uint32)) & 1).astype(np.int16)
        spins = np.concatenate((np.ones((len(codes), 1), np.int16), 1 - 2 * bits), axis=1)
        sums = {key: 0.0 for key in ["bilinear", "nuclear", "projector", "defect_projector", "exact", "defect_exact", "size_i", "size_i2", "collapsed_nuclear", "collapsed_projector", "collapsed_exact", "h_l2_sq", "h_max"]}
        masks: dict[tuple[int, ...], tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = {}
        records = []
        for x in spins:
            field = matrix @ x
            y = np.sign(field).astype(np.int16)
            g = x * y
            ii = np.flatnonzero(g > 0)
            jj = np.flatnonzero(g < 0)
            switched = matrix * x[:, None] * x[None, :]
            bi = switched[np.ix_(ii, ii)].astype(float)
            bj = switched[np.ix_(jj, jj)].astype(float)
            p = int(np.ones(len(ii), int) @ bi @ np.ones(len(ii), int)) if len(ii) else 0
            r = int(np.ones(len(jj), int) @ bj @ np.ones(len(jj), int)) if len(jj) else 0
            response = int(np.abs(field).sum())
            if p * r < 0:
                pi = projector(bi, 1 if r > 0 else -1)
                pj = projector(bj, 1 if p > 0 else -1)
                ni = nuclear(bi); nj = nuclear(bj)
                ei = one_sided(bi.astype(np.int16), 1 if r > 0 else -1)
                ej = one_sided(bj.astype(np.int16), 1 if p > 0 else -1)
                defect_p = min(max(0, abs(p) - pi), max(0, abs(r) - pj))
                defect_e = min(max(0, abs(p) - ei), max(0, abs(r) - ej))
            else:
                pi = pj = ni = nj = defect_p = defect_e = 0.0
                ei = ej = 0
            # J collapsed over I. The exact cap is feasible through n=18 because J is typically small;
            # skip enormous degenerate shores and use zero where recoupling already has zero defect.
            if len(jj) and len(ii):
                h = switched[np.ix_(jj, ii)] @ np.ones(len(ii))
                e = np.zeros((len(jj) + 1, len(jj) + 1), float)
                e[:-1, :-1] = bj; e[:-1, -1] = h; e[-1, :-1] = h
                cn = max(0.0, np.abs(np.linalg.eigvalsh(e)).sum() / math.pi - (1 - 2 / math.pi) * len(e) * max(1.0, np.abs(e).max()))
                cp = projector(e, 1 if p > 0 else -1)
                ce = one_sided(e.astype(np.int16), 1 if p > 0 else -1)
                hl2 = float(h @ h); hmax = float(np.max(np.abs(h)))
            else:
                cn = cp = ce = hl2 = hmax = 0.0
            data = dict(bilinear=response,nuclear=ni+nj,projector=pi+pj,defect_projector=defect_p,exact=ei+ej,defect_exact=defect_e,size_i=len(ii),size_i2=len(ii)**2,collapsed_nuclear=cn,collapsed_projector=cp,collapsed_exact=ce,h_l2_sq=hl2,h_max=hmax)
            for key, value in data.items(): sums[key] += value
        count = len(spins)
        means = {key: value / count for key, value in sums.items()}
        m=(n-2)//2
        eta=math.comb(2*m,m)/2**(2*m)
        rho=(math.comb(m,m//2)/2**m)**2 if m%2==0 else 0
        output={"file":filename,"n":n,"samples":count,"means":means,"theory":{"bilinear":n*(n-1)*eta,"size_i":n/2,"size_i2":n*n/4+(n+n*(n-1)*rho)/4}}
        print(json.dumps(output,sort_keys=True),flush=True)


if __name__ == "__main__":
    main()
