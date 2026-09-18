#!/usr/bin/env python3
"""Exact two-near-state polarization and shore-resource audit (Wave 17).

All arithmetic is integral.  Oriented states are projective Boolean spins.
Scratch-only: this file is intentionally under the ignored tmp/ directory.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from itertools import product

import numpy as np

from check_finite_bridge_r16 import (
    A8,
    A9,
    EXPECTED_CODES,
    V,
    all_minimizers,
    edge_code,
    energies,
    hybrid,
    qnorm,
    spins,
    wall_matrices,
)


A8_BLOCKS = ((0, 3, 4, 5), (1, 2, 6, 7))


def oriented_states(M: np.ndarray) -> list[dict]:
    q = qnorm(M)
    e = energies(M)
    return [
        {"sigma": sigma, "z": z.copy(), "e": sigma * int(value), "delta": q - sigma * int(value)}
        for z, value in zip(spins(len(M)), e)
        for sigma in (-1, 1)
    ]


def near_states(C: np.ndarray, t: int) -> dict[int, list[dict]]:
    by_orientation: dict[int, list[dict]] = {-1: [], 1: []}
    for row in oriented_states(C):
        if row["delta"] <= t:
            by_orientation[row["sigma"]].append(row)
    return by_orientation


def shore_data(A: np.ndarray, sigma: int, z: np.ndarray, w: np.ndarray) -> dict:
    """Compute the exact identities for the opposite-oriented pair.

    The first state is (sigma,z), the second (-sigma,w).  S is their agreement
    shore and T their disagreement shore.  The raw tolled residual sums both
    orientations on both shores; it is not an allocated endpoint resource.
    """
    n = len(A)
    q = qnorm(A)
    S = tuple(i for i in range(n) if int(z[i] * w[i]) == 1)
    T = tuple(i for i in range(n) if int(z[i] * w[i]) == -1)
    ez = sigma * int(z @ A @ z)
    ew = -sigma * int(w @ A @ w)
    alpha = q - ez
    beta = q - ew
    b = sigma * int(z @ A @ w)
    H = (beta - alpha) // 2
    assert 2 * H == beta - alpha

    def h(X: tuple[int, ...]) -> int:
        if not X:
            return 0
        zx = z[list(X)]
        return sigma * int(zx @ A[np.ix_(X, X)] @ zx)

    hS, hT = h(S), h(T)
    assert hS + hT == H
    assert hS - hT == b

    def cross_exposure(X: tuple[int, ...], Y: tuple[int, ...]) -> tuple[int, int]:
        if not X or not Y:
            return 0, 0
        zx, zy = z[list(X)], z[list(Y)]
        B = A[np.ix_(Y, X)]
        gamma = sigma * int(zy @ B @ zx)
        L = int(np.abs(B @ zx).sum())
        return gamma, L

    gammaS, LS = cross_exposure(S, T)
    gammaT, LT = cross_exposure(T, S)
    assert gammaS == gammaT
    gamma = gammaS
    assert 4 * gamma == 2 * q - alpha - beta
    assert LS >= abs(gamma) and LT >= abs(gamma)
    ellS, ellT = LS - gamma, LT - gamma
    assert ellS >= 0 and ellT >= 0

    def pos(x: int) -> int:
        return max(x, 0)

    rSplus = pos(2 * LS + hS - q)
    rSminus = pos(2 * LS - hS - q)
    rTplus = pos(2 * LT + hT - q)
    rTminus = pos(2 * LT - hT - q)
    residual = rSplus + rSminus + rTplus + rTminus

    # Exact deficit/polarization rewrite of all four residuals.
    assert rSplus == pos(2 * ellS - hT - alpha)
    assert rSminus == pos(2 * ellS + hT - beta)
    assert rTplus == pos(2 * ellT - hS - alpha)
    assert rTminus == pos(2 * ellT + hS - beta)

    d = alpha + beta
    lower_exact = (
        pos(abs(b - H) - d) + pos(abs(b + H) - d)
    ) // 2
    # Each bracket is even in the signing normalization.
    assert 2 * lower_exact == pos(abs(b - H) - d) + pos(abs(b + H) - d)
    lower_simple = pos(abs(b) - d)
    assert residual >= lower_exact >= lower_simple
    return {
        "S": S,
        "T": T,
        "alpha": alpha,
        "beta": beta,
        "defsum": d,
        "b": b,
        "H": H,
        "hS": hS,
        "hT": hT,
        "gamma": gamma,
        "LS": LS,
        "LT": LT,
        "ellS": ellS,
        "ellT": ellT,
        "residual": residual,
        "lower_exact": lower_exact,
        "lower_simple": lower_simple,
        "bucket_residuals": (rSplus, rSminus, rTplus, rTminus),
    }


def pair_audit(A: np.ndarray, C: np.ndarray, t: int) -> dict:
    q, qc = qnorm(A), qnorm(C)
    near = near_states(C, t)
    pair_rows = []
    for sigma in (-1, 1):
        for zr in near[sigma]:
            z = zr["z"]
            for wr in near[-sigma]:
                w = wr["z"]
                row = shore_data(A, sigma, z, w)
                row.update(
                    {
                        "sigma": sigma,
                        "z": tuple(map(int, z)),
                        "w": tuple(map(int, w)),
                        "u": zr["delta"],
                        "v": wr["delta"],
                    }
                )
                # A-deficits expressed through C-deficits and D=A-C.
                D = A - C
                dz = sigma * int(z @ D @ z)
                dw = -sigma * int(w @ D @ w)
                assert row["alpha"] == q - qc + row["u"] - dz
                assert row["beta"] == q - qc + row["v"] - dw
                row["dz"] = dz
                row["dw"] = dw
                row["response_sum"] = (dz - row["u"]) + (dw - row["v"])
                assert row["defsum"] == 2 * (q - qc) - row["response_sum"]
                pair_rows.append(row)

    single_response = {}
    for sigma in (-1, 1):
        vals = []
        D = A - C
        for row in near[sigma]:
            z = row["z"]
            vals.append(sigma * int(z @ D @ z) - row["delta"])
        single_response[sigma] = max(vals) if vals else None
    finite = [x for x in single_response.values() if x is not None]
    R = max(finite) if finite else None
    Gamma = max(-R, 0) if R is not None else None

    def maximum(key: str):
        return max((r[key] for r in pair_rows), default=None)

    def minimum(key: str):
        return min((r[key] for r in pair_rows), default=None)

    best_lb = max(pair_rows, key=lambda r: (r["lower_simple"], r["residual"])) if pair_rows else None
    best_raw = max(pair_rows, key=lambda r: r["residual"]) if pair_rows else None
    best_def = min(pair_rows, key=lambda r: (r["defsum"], -abs(r["b"]))) if pair_rows else None
    return {
        "q": q,
        "qc": qc,
        "t": t,
        "near_counts": {sigma: len(rows) for sigma, rows in near.items()},
        "pair_count": len(pair_rows),
        "single_response": single_response,
        "R": R,
        "Gamma": Gamma,
        "min_defsum": minimum("defsum"),
        "max_response_sum": maximum("response_sum"),
        "max_abs_b": max((abs(r["b"]) for r in pair_rows), default=None),
        "max_lower_simple": maximum("lower_simple"),
        "max_lower_exact": maximum("lower_exact"),
        "max_raw_residual": maximum("residual"),
        "best_lb": best_lb,
        "best_raw": best_raw,
        "best_def": best_def,
        "pair_rows": pair_rows,
    }


def compact(audit: dict) -> dict:
    keys = (
        "q",
        "qc",
        "t",
        "near_counts",
        "pair_count",
        "single_response",
        "R",
        "Gamma",
        "min_defsum",
        "max_response_sum",
        "max_abs_b",
        "max_lower_simple",
        "max_lower_exact",
        "max_raw_residual",
    )
    return {k: audit[k] for k in keys}


def a8_cross() -> np.ndarray:
    C = A8.copy()
    for block in A8_BLOCKS:
        C[np.ix_(block, block)] = 0
    return C


def a8_best_hybrids(C: np.ndarray) -> list[np.ndarray]:
    q4, minimizers = all_minimizers(4)
    assert q4 == 8
    best = 10**9
    out: list[np.ndarray] = []
    for G0 in minimizers:
        for G1 in minimizers:
            H = C.copy()
            H[np.ix_(A8_BLOCKS[0], A8_BLOCKS[0])] = G0
            H[np.ix_(A8_BLOCKS[1], A8_BLOCKS[1])] = G1
            qh = qnorm(H)
            if qh < best:
                best, out = qh, [H]
            elif qh == best:
                out.append(H)
    assert best == 20 and len(out) == 16
    return out


def main() -> None:
    C8 = a8_cross()
    e8 = energies(C8)
    assert (max(map(int, e8)), -min(map(int, e8)), qnorm(C8)) == (16, 16, 16)
    a8 = pair_audit(A8, C8, 0)
    assert compact(a8) == {
        "q": 20,
        "qc": 16,
        "t": 0,
        "near_counts": {-1: 12, 1: 12},
        "pair_count": 288,  # directed; 144 after fixing sigma=+1
        "single_response": {-1: 4, 1: 4},
        "R": 4,
        "Gamma": 0,
        "min_defsum": 0,
        "max_response_sum": 8,
        "max_abs_b": 16,
        "max_lower_simple": 8,
        "max_lower_exact": 8,
        "max_raw_residual": 8,
    }
    a8_response_pairs = [
        r
        for r in a8["pair_rows"]
        if r["sigma"] == 1
        and r["dz"] - r["u"] == a8["single_response"][1]
        and r["dw"] - r["v"] == a8["single_response"][-1]
    ]
    assert len(a8_response_pairs) == 16
    assert Counter((abs(r["b"]), r["residual"]) for r in a8_response_pairs) == Counter({(0, 0): 16})
    a8_best_pairs = [r for r in a8["pair_rows"] if r["sigma"] == 1 and r["residual"] == 8]
    assert len(a8_best_pairs) == 4
    assert all(
        (r["u"], r["v"], r["alpha"], r["beta"], abs(r["b"]), r["response_sum"])
        == (0, 0, 4, 4, 16, 0)
        for r in a8_best_pairs
    )
    # These profitable cuts have two exact q4 shores and raw buckets (4,0,0,4),
    # but they are not endpoint cuts and carry no pre-existing allocation.
    for r in a8_best_pairs:
        assert tuple(qnorm(A8[np.ix_(X, X)]) for X in (r["S"], r["T"])) == (8, 8)
        assert sorted(r["bucket_residuals"]) == [0, 0, 4, 4]
    print("A8 t=0:", compact(a8))
    print("A8 response-max pair histogram (|b|,raw):", Counter((abs(r["b"]), r["residual"]) for r in a8_response_pairs))
    print("A8 max-raw pair count:", len(a8_best_pairs))

    C9, _ = wall_matrices()
    q6, minimizers = all_minimizers(6)
    assert q6 == 10
    completions = [G for G in minimizers if qnorm(hybrid(C9, G)) == 24]
    assert tuple(sorted(edge_code(G) for G in completions)) == EXPECTED_CODES
    e9 = energies(C9)
    assert (max(map(int, e9)), -min(map(int, e9)), qnorm(C9)) == (18, 22, 22)
    print("A9 40 completions")
    for t in (0, 4, 8):
        summaries = Counter()
        response_pair_counts = Counter()
        response_pair_hist = Counter()
        for G in completions:
            H = hybrid(C9, G)
            audit = pair_audit(H, C9, t)
            summaries[repr(compact(audit))] += 1
            if t == 0:
                assert audit["near_counts"] == {-1: 3, 1: 0}
                assert audit["R"] == -6 and audit["Gamma"] == 6 and not audit["pair_rows"]
                continue
            response_pairs = [
                r
                for r in audit["pair_rows"]
                if r["sigma"] == 1
                and r["dz"] - r["u"] == audit["single_response"][1]
                and r["dw"] - r["v"] == audit["single_response"][-1]
            ]
            assert all((r["alpha"], r["beta"]) == (0, 0) for r in response_pairs)
            response_pair_counts[len(response_pairs)] += 1
            response_pair_hist.update((abs(r["b"]), r["residual"]) for r in response_pairs)
            assert audit["single_response"] == {-1: 2, 1: 2}
            assert audit["min_defsum"] == 0 and audit["max_response_sum"] == 4
            assert audit["max_raw_residual"] == (8 if t == 4 else 16)
        if t == 0:
            assert len(summaries) == 1 and next(iter(summaries.values())) == 40
        elif t == 4:
            assert response_pair_counts == Counter({12: 24, 24: 16})
            assert response_pair_hist == Counter({(0, 0): 336, (4, 4): 192, (8, 8): 144})
        else:
            assert response_pair_counts == Counter({108: 24, 156: 16})
            assert response_pair_hist == Counter(
                {(0, 0): 1944, (4, 4): 1152, (8, 8): 1896, (16, 16): 96}
            )
        print("t", t, "summary profiles", summaries)
        if t:
            print("  response-max pair counts/completion:", response_pair_counts)
            print("  response-max histogram (|b|,raw):", response_pair_hist)

    print("PASS: every polarization, bucket, response, and census assertion")


if __name__ == "__main__":
    main()
