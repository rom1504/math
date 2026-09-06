#!/usr/bin/env python3
"""Wave 25 endpoint-transport identities and finite audits.

The script enumerates the full oriented cut space.  Its Hellinger convention is

    H^2(mu,nu) = 1/2 sum_x (sqrt(mu_x)-sqrt(nu_x))^2

also for finite measures of unequal mass.  Thus the raw Johnson ordinary-LS
functional has coefficient m(n-m)/(2n), exactly as printed below.
"""

from __future__ import annotations

from itertools import product
import math

import numpy as np
from scipy.special import logsumexp


A4 = np.array([
    [0, -1, 1, 1],
    [-1, 0, 1, -1],
    [1, 1, 0, -1],
    [1, -1, -1, 0],
], dtype=np.int64)

A6 = np.array([
    [0, 1, 1, 1, 1, 1],
    [1, 0, -1, -1, 1, 1],
    [1, -1, 0, 1, -1, 1],
    [1, -1, 1, 0, 1, -1],
    [1, 1, -1, 1, 0, -1],
    [1, 1, 1, -1, -1, 0],
], dtype=np.int64)

A8 = np.array([
    [0, 1, 1, 1, 1, 1, 1, 1],
    [1, 0, 1, -1, 1, 1, -1, -1],
    [1, 1, 0, 1, -1, 1, -1, -1],
    [1, -1, 1, 0, -1, -1, -1, 1],
    [1, 1, -1, -1, 0, -1, 1, -1],
    [1, 1, 1, -1, -1, 0, 1, 1],
    [1, -1, -1, -1, 1, 1, 0, 1],
    [1, -1, -1, 1, -1, 1, 1, 0],
], dtype=np.int64)

A9 = np.array([
    [0, 1, -1, 1, 1, 1, 1, -1, -1],
    [1, 0, 1, -1, -1, -1, 1, -1, -1],
    [-1, 1, 0, 1, 1, 1, 1, 1, -1],
    [1, -1, 1, 0, 1, 1, 1, -1, 1],
    [1, -1, 1, 1, 0, 1, -1, 1, -1],
    [1, -1, 1, 1, 1, 0, -1, -1, -1],
    [1, 1, 1, 1, -1, -1, 0, 1, 1],
    [-1, -1, 1, -1, 1, -1, 1, 0, -1],
    [-1, -1, -1, 1, -1, -1, 1, -1, 0],
], dtype=np.int64)


def oriented_cuts(n: int) -> np.ndarray:
    """All 2^n oriented rank-one cuts as symmetric zero-diagonal matrices."""
    ans = []
    for sigma in (-1, 1):
        for tail in product((-1, 1), repeat=n - 1):
            x = np.asarray((1,) + tail, dtype=np.int64)
            d = sigma * np.outer(x, x)
            np.fill_diagonal(d, 0)
            ans.append(d)
    return np.asarray(ans, dtype=np.int64)


def restriction_key(d: np.ndarray, ids: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(int(d[ids[i], ids[j]]) for i in range(len(ids)) for j in range(i + 1, len(ids)))


def entropy_term(x: np.ndarray) -> np.ndarray:
    return np.where(x > 0, x * np.log(x), 0.0)


def tau_ls_upper(n: int, m: int) -> float:
    """The ordinary Johnson log-Sobolev constant used in the memo."""
    return (2.0 / math.log(2.0)) * math.log(n * n / (m * (n - m)))


def endpoint_metrics(a_mat: np.ndarray, beta: float, gamma: float) -> dict[str, float]:
    n = len(a_mat)
    m = n - 1
    selectors = [tuple(j for j in range(n) if j != omitted) for omitted in range(n)]
    cuts = oriented_cuts(n)
    energies = np.einsum("bij,ij->b", cuts, a_mat, optimize=True).astype(float)
    log_za = float(logsumexp(beta * energies))
    nu = np.exp(beta * energies - log_za)

    num_s, num_d = len(selectors), len(cuts)
    c = np.empty((num_s, num_d), dtype=float)
    ell = np.empty_like(c)
    log_zs = np.empty(num_s, dtype=float)
    q = np.empty_like(c)

    # Build each conditional partition K_{beta,S}(y) by grouping full cuts by
    # their restriction.  Every child cut occurs with equal multiplicity.
    for si, ids in enumerate(selectors):
        c[si] = np.asarray([
            np.sum(a_mat[np.ix_(ids, ids)] * d[np.ix_(ids, ids)]) for d in cuts
        ])
        groups: dict[tuple[int, ...], list[int]] = {}
        for di, d in enumerate(cuts):
            groups.setdefault(restriction_key(d, ids), []).append(di)
        assert len(groups) == 2**m

        child_cs = np.asarray([c[si, ds[0]] for ds in groups.values()])
        log_zs[si] = float(logsumexp(gamma * child_cs))
        for ds in groups.values():
            log_k = float(logsumexp(beta * (energies[ds] - c[si, ds])))
            ell[si, ds] = log_k + (beta - gamma) * c[si, ds]
            q[si, ds] = np.exp(
                gamma * c[si, ds] - log_zs[si]
                + beta * (energies[ds] - c[si, ds]) - log_k
            )
        assert abs(float(np.sum(q[si])) - 1.0) < 2e-12

    log_z0 = float(logsumexp(log_zs) - math.log(num_s))
    cap_c = math.exp(log_za - log_z0)
    h = np.exp(-ell)
    av_h = np.mean(h, axis=0)
    f_cut = cap_c * av_h
    assert abs(float(np.dot(nu, f_cut)) - 1.0) < 3e-11

    ent_cut = float(np.dot(nu, entropy_term(f_cut)))
    rho = h / (num_s * av_h[None, :])
    selector_kl_by_d = np.sum(np.where(rho > 0, rho * np.log(rho * num_s), 0.0), axis=0)
    selector_entropy = float(np.dot(nu * f_cut, selector_kl_by_d))
    total = ent_cut + selector_entropy
    total_direct = float(np.sum((nu[None, :] / num_s) * (cap_c * h) * np.log(cap_c * h)))
    assert abs(total - total_direct) < 2e-11

    # The lifted endpoint identity C nu(d) h_S(d) = w_S q_S(d).
    w = np.exp(log_zs - log_z0)
    lhs = cap_c * nu[None, :] * h
    rhs = w[:, None] * q
    assert np.max(np.abs(lhs - rhs)) < 3e-12

    # Directed Johnson edges.  For m=n-1 this is every ordered S != T.
    directed = [(s, t) for s in range(num_s) for t in range(num_s) if s != t]
    for s, t in directed:
        jeff_pointwise = (h[t] - h[s]) * (ell[s] - ell[t])
        hell_pointwise = (np.sqrt(h[t]) - np.sqrt(h[s])) ** 2
        assert np.min(jeff_pointwise - 4.0 * hell_pointwise) > -2e-12
    edge_av_jeff = np.mean([
        np.dot(nu, (h[t] - h[s]) * (ell[s] - ell[t])) for s, t in directed
    ])
    edge_av_sqrt = np.mean([
        np.dot(nu, (np.sqrt(h[t]) - np.sqrt(h[s])) ** 2) for s, t in directed
    ])
    coeff = m * (n - m) / (2.0 * n)
    d_mls_raw = coeff * cap_c * edge_av_jeff
    h_ls_raw = coeff * cap_c * edge_av_sqrt
    assert d_mls_raw + 3e-12 >= 4.0 * h_ls_raw

    tau = tau_ls_upper(n, m)
    ls_rhs_upper = tau * h_ls_raw
    # This is an actual proved upper RHS, so it must dominate the exact
    # conditional selector entropy (up to floating error).
    assert selector_entropy <= ls_rhs_upper + 2e-10, (
        n, beta, selector_entropy, ls_rhs_upper
    )

    # Audit finite-measure Hellinger and its mass/channel decomposition on
    # every edge.  H^2(mu,nu)=1/2 sum(sqrt(mu)-sqrt(nu))^2.
    finite_h_av = 0.0
    for s, t in directed:
        bc = float(np.sum(np.sqrt((w[s] * q[s]) * (w[t] * q[t]))))
        finite_h = (w[s] + w[t]) / 2.0 - bc
        finite_h_square = 0.5 * float(np.sum((np.sqrt(w[s] * q[s]) - np.sqrt(w[t] * q[t])) ** 2))
        assert abs(finite_h - finite_h_square) < 2e-12
        bc_q = float(np.sum(np.sqrt(q[s] * q[t])))
        h_q = 1.0 - bc_q
        decomposed = 0.5 * (math.sqrt(w[s]) - math.sqrt(w[t])) ** 2 + math.sqrt(w[s] * w[t]) * h_q
        assert abs(finite_h - decomposed) < 2e-12
        finite_h_av += finite_h / len(directed)
    assert abs(h_ls_raw - (m * (n - m) / n) * finite_h_av) < 3e-12

    # Check the adjacent conditional-probability formula
    # ell_S-ell_T = log P(D_{S\C}|D_C)-log P(D_{T\C}|D_C)
    #               - gamma(c_S-c_T).
    for s, t in directed:
        core = tuple(sorted(set(selectors[s]).intersection(selectors[t])))
        p_core: dict[tuple[int, ...], float] = {}
        for di, d in enumerate(cuts):
            key = restriction_key(d, core)
            p_core[key] = p_core.get(key, 0.0) + nu[di]
        p_s: dict[tuple[int, ...], float] = {}
        p_t: dict[tuple[int, ...], float] = {}
        for di, d in enumerate(cuts):
            ks = restriction_key(d, selectors[s])
            kt = restriction_key(d, selectors[t])
            p_s[ks] = p_s.get(ks, 0.0) + nu[di]
            p_t[kt] = p_t.get(kt, 0.0) + nu[di]
        for di, d in enumerate(cuts):
            kc = restriction_key(d, core)
            cond_s = p_s[restriction_key(d, selectors[s])] / p_core[kc]
            cond_t = p_t[restriction_key(d, selectors[t])] / p_core[kc]
            predicted = math.log(cond_s / cond_t) - gamma * (c[s, di] - c[t, di])
            assert abs((ell[s, di] - ell[t, di]) - predicted) < 2e-11

    return {
        "ent_cut": ent_cut,
        "selector_entropy": selector_entropy,
        "d_mls_raw": d_mls_raw,
        "h_ls_raw": h_ls_raw,
        "tau_ls_upper": tau,
        "ls_rhs_upper": ls_rhs_upper,
        "total": total,
    }


def main() -> None:
    assert int(np.max(np.einsum("bij,ij->b", oriented_cuts(4), A4))) == 8

    print("Convention: H^2=1/2 sum(sqrt(mu)-sqrt(nu))^2")
    print("Columns: beta Ent_cut selector_KL Jeff_raw H_raw tau_LS*tH_raw total")
    all_rows: dict[str, list[tuple[float, dict[str, float]]]] = {}
    for name, matrix in (("A6", A6), ("A8", A8), ("A9", A9)):
        rows = []
        print(name, "matched beta=gamma, m=n-1")
        for beta in (0.1, 0.5, 1.0, 2.0, 4.0, 8.0):
            z = endpoint_metrics(matrix, beta, beta)
            rows.append((beta, z))
            print(
                f"{beta:4.1f} {z['ent_cut']:.9f} {z['selector_entropy']:.9f} "
                f"{z['d_mls_raw']:.9f} {z['h_ls_raw']:.9f} "
                f"{z['ls_rhs_upper']:.9f} {z['total']:.9f}"
            )
        all_rows[name] = rows

    # At low temperature the Jeffreys form continues to grow while both the
    # actual entropy and the Hellinger form have essentially saturated.
    for name, rows in all_rows.items():
        at4 = rows[-2][1]
        at8 = rows[-1][1]
        assert at8["d_mls_raw"] > 1.9 * at4["d_mls_raw"]
        assert abs(at8["h_ls_raw"] - at4["h_ls_raw"]) < 4e-4
        assert abs(at8["selector_entropy"] - at4["selector_entropy"]) < 1e-4
        print(
            f"{name} beta 4->8 ratios: Jeff={at8['d_mls_raw']/at4['d_mls_raw']:.6f}, "
            f"H={at8['h_ls_raw']/at4['h_ls_raw']:.6f}, "
            f"selector={at8['selector_entropy']/at4['selector_entropy']:.6f}"
        )

    z4 = endpoint_metrics(A4, 0.5, 0.1)
    print("A4 exact minimizer, beta=.5 gamma=.1:", z4)
    assert abs(z4["ent_cut"] - 0.6276970227309644) < 2e-11
    assert abs(z4["selector_entropy"] - 0.14716108819833879) < 2e-11
    assert abs(z4["d_mls_raw"] - 0.3636529347339671) < 2e-11
    assert abs(z4["h_ls_raw"] - 0.0804378663432645) < 2e-11
    assert abs(z4["total"] - 0.7748581109293031) < 2e-11
    assert z4["ent_cut"] > z4["d_mls_raw"]
    assert z4["ent_cut"] > z4["ls_rhs_upper"]

    print("all endpoint identities and inequalities verified")


if __name__ == "__main__":
    main()
