#!/usr/bin/env python3
"""Wave 35 checks for endpoint information and orientation marginalization.

All calculations stay in the repository-local scratch directory.  Floating
enumeration verifies the KL chain rules for the canonical harmonic mixture.
An independent Fraction calculation gives an exact orientation counterexample
at beta=(log 2)/2.
"""

from __future__ import annotations

import itertools
import math
import sys
from fractions import Fraction

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import logsumexp
from scipy.special import expit

sys.path.insert(0, "/home/math/quadra/tmp")
from endpoint_transport_r25 import A4, A6, A8, A9, oriented_cuts, restriction_key
from harmonic_resonance_r34_check import Endpoint


def kl(p: np.ndarray, q: np.ndarray) -> float:
    mask = p > 0
    return float(np.sum(p[mask] * np.log(p[mask] / q[mask])))


def information_audit(A: np.ndarray, beta: float, m: int) -> dict[str, object]:
    H = Endpoint(A, beta, m)
    mu = H.nu * H.f
    component = H.nu[None, :] * H.h
    q = H.q

    dfull = kl(mu, H.nu)
    component_d = np.asarray([kl(component[s], H.nu) for s in range(len(q))])
    joint_d = float(q @ component_d)
    ifull = joint_d - dfull

    C = []
    J = []
    Icond = []
    Dminus = []
    Iminus = []
    omitted_component_cost = []
    for bit, edges in enumerate(H.edges):
        nu_y = np.asarray([H.nu[x] + H.nu[y] for x, y in edges])
        mu_y = np.asarray([mu[x] + mu[y] for x, y in edges])
        dminus = kl(mu_y, nu_y)
        Dminus.append(dminus)
        C.append(dfull - dminus)

        comp_y = np.asarray(
            [[component[s, x] + component[s, y] for x, y in edges]
             for s in range(len(q))]
        )
        comp_dminus = np.asarray([kl(comp_y[s], nu_y) for s in range(len(q))])
        J.append(float(q @ (component_d - comp_dminus)))

        # I(S;D_-j) under P(S,D)=q(S)mu_S(D).
        iy = 0.0
        for s in range(len(q)):
            positive = comp_y[s] > 0
            iy += float(q[s] * np.sum(
                comp_y[s, positive]
                * np.log(comp_y[s, positive] / mu_y[positive])
            ))
        Iminus.append(iy)
        Icond.append(ifull - iy)

        # For vertex bit i>=1, h_S is bit-invariant whenever i is omitted.
        max_omit = 0.0
        if bit >= 1:
            for s, S in enumerate(H.selectors):
                if bit not in S:
                    max_omit = max(max_omit, abs(component_d[s] - comp_dminus[s]))
        omitted_component_cost.append(max_omit)

    C = np.asarray(C)
    J = np.asarray(J)
    Icond = np.asarray(Icond)
    Dminus = np.asarray(Dminus)
    Iminus = np.asarray(Iminus)
    assert np.max(np.abs(C - (J - Icond))) < 2e-9
    assert abs(float(np.sum(C)) - (H.n * dfull - float(np.sum(Dminus)))) < 2e-9
    assert np.max(np.asarray(omitted_component_cost)) < 2e-9

    # Marginalizing the orientation bit has exact entropy loss C_0.
    quotient_ent = float(Dminus[0])
    assert abs(dfull - quotient_ent - C[0]) < 2e-9

    return {
        "n": H.n,
        "m": m,
        "beta": beta,
        "D": dfull,
        "Dquot": quotient_ent,
        "C": C,
        "Csum": float(np.sum(C)),
        "C0": float(C[0]),
        "Jsum": float(np.sum(J)),
        "Ifull": ifull,
        "Icondsum": float(np.sum(Icond)),
        "identity_error": float(np.max(np.abs(C - (J - Icond)))),
    }


def pow2_fraction(k: int) -> Fraction:
    return Fraction(2**k, 1) if k >= 0 else Fraction(1, 2 ** (-k))


def exact_U_log2_half(A: np.ndarray, m: int) -> tuple[list[np.ndarray], list[Fraction]]:
    """Return exact U_beta for beta=(log 2)/2 on all oriented cuts."""
    A = np.asarray(A, dtype=int)
    n = len(A)
    cuts = list(oriented_cuts(n))
    energies = [int(np.sum(A * d)) for d in cuts]
    selectors = list(itertools.combinations(range(n), m))
    accum = [Fraction(0, 1) for _ in cuts]
    for S in selectors:
        groups: dict[tuple[int, ...], list[int]] = {}
        for di, d in enumerate(cuts):
            groups.setdefault(restriction_key(d, S), []).append(di)
        for inds in groups.values():
            child = int(np.sum(A[np.ix_(S, S)] * cuts[inds[0]][np.ix_(S, S)]))
            external = [energies[di] - child for di in inds]
            assert all(x % 2 == 0 for x in external)
            denom = sum((pow2_fraction(x // 2) for x in external), Fraction(0, 1))
            e_minus_F = Fraction(len(inds), 1) / denom
            for di in inds:
                accum[di] += e_minus_F
    U = [x / len(selectors) for x in accum]
    return cuts, U


def exact_orientation_counterexample() -> dict[str, object]:
    cuts, U = exact_U_log2_half(A6, 3)
    by_key = {restriction_key(d, tuple(range(6))): i for i, d in enumerate(cuts)}
    for i, d in enumerate(cuts):
        j = by_key[restriction_key(-d, tuple(range(6)))]
        if U[i] != U[j]:
            ratio = U[i] / U[j]
            # Endpoint/base conditional odds ratio is exactly f(d)/f(-d)=U(d)/U(-d).
            assert ratio != 1
            answer = {
                "i": i,
                "j": j,
                "energy_i": int(np.sum(A6 * d)),
                "energy_j": int(np.sum(A6 * (-d))),
                "U_i": U[i],
                "U_j": U[j],
                "ratio": ratio,
                "key_i": restriction_key(d, tuple(range(6))),
            }
            assert answer["i"] == 0 and answer["j"] == 32
            assert answer["U_i"] == Fraction(847888, 3590575)
            assert answer["U_j"] == Fraction(111448, 963325)
            assert answer["ratio"] == Fraction(1165846, 571171)
            return answer
    raise AssertionError("unexpected orientation invariance")


def migration_state(H: Endpoint, t: float) -> dict[str, np.ndarray | float]:
    logmu = H.lognu + t * H.g
    logmu -= float(logsumexp(logmu))
    mu = np.exp(logmu)
    d_full = kl(mu, H.nu)
    mean_g = float(mu @ H.g)
    var_g = float(mu @ (H.g - mean_g) ** 2)

    qdim = len(H.edges)
    d_marg = np.zeros(qdim)
    c_cond = np.zeros(qdim)
    cov = np.zeros(qdim)
    var_context_score = np.zeros(qdim)
    energy = np.zeros(qdim)
    d_marg_derivative = np.zeros(qdim)
    c_derivative = np.zeros(qdim)
    for bit, edges in enumerate(H.edges):
        x = edges[:, 0]
        y = edges[:, 1]
        M = mu[x] + mu[y]
        M0 = H.nu[x] + H.nu[y]
        d_marg[bit] = kl(M, M0)
        c_cond[bit] = d_full - d_marg[bit]
        u0 = H.lognu[y] - H.lognu[x]
        ut = u0 + t * (H.g[y] - H.g[x])
        p = expit(ut)
        Aprime = (1.0 - p) * H.g[x] + p * H.g[y]
        log_context_likelihood = np.log(M / M0)
        kval = p * (ut - u0) - (np.logaddexp(0.0, ut) - np.logaddexp(0.0, u0))
        mean_A = float(M @ Aprime)
        mean_k = float(M @ kval)
        cov[bit] = float(M @ ((Aprime - mean_A) * (kval - mean_k)))
        var_context_score[bit] = float(M @ (Aprime - mean_A) ** 2)
        energy[bit] = float(M @ (expit(ut) * expit(-ut) * (H.g[y] - H.g[x]) ** 2))
        d_marg_derivative[bit] = float(
            M @ ((Aprime - mean_A) * (log_context_likelihood - d_marg[bit]))
        )
        c_derivative[bit] = t * var_g - d_marg_derivative[bit]

    assert np.max(np.abs(cov - (t * var_context_score - d_marg_derivative))) < 3e-9
    assert np.max(np.abs(c_derivative - (t * energy + cov))) < 3e-9
    return {
        "d_full": d_full,
        "var_g": var_g,
        "d_marg": d_marg,
        "c_cond": c_cond,
        "cov": cov,
        "var_context_score": var_context_score,
        "energy": energy,
        "d_marg_derivative": d_marg_derivative,
        "c_derivative": c_derivative,
    }


def migration_audit(A: np.ndarray, beta: float, m: int, order: int = 64) -> dict[str, object]:
    H = Endpoint(A, beta, m)
    nodes, weights = leggauss(order)
    ts = (nodes + 1.0) / 2.0
    ws = weights / 2.0
    signed_cov = np.zeros(H.n)
    adverse_cov = np.zeros(H.n)
    positive_cov = np.zeros(H.n)
    weighted_energy = np.zeros(H.n)
    positive_d_marg_variation = np.zeros(H.n)
    negative_d_marg_variation = np.zeros(H.n)
    min_d_marg_derivative = np.full(H.n, math.inf)
    min_c_derivative = np.full(H.n, math.inf)
    for t, w in zip(ts, ws):
        z = migration_state(H, float(t))
        cov = np.asarray(z["cov"])
        dmd = np.asarray(z["d_marg_derivative"])
        cd = np.asarray(z["c_derivative"])
        signed_cov += w * cov
        adverse_cov += w * np.maximum(-cov, 0.0)
        positive_cov += w * np.maximum(cov, 0.0)
        weighted_energy += w * t * np.asarray(z["energy"])
        positive_d_marg_variation += w * np.maximum(dmd, 0.0)
        negative_d_marg_variation += w * np.maximum(-dmd, 0.0)
        min_d_marg_derivative = np.minimum(min_d_marg_derivative, dmd)
        min_c_derivative = np.minimum(min_c_derivative, cd)
    end = migration_state(H, 1.0)
    C = np.asarray(end["c_cond"])
    Dmarg = np.asarray(end["d_marg"])
    assert np.max(np.abs(C - (weighted_energy + signed_cov))) < 3e-8
    # Pointwise [-Cov]_+ <= [D'_marg]_+ follows from
    # -Cov=D'_marg-t Var(A') and t Var(A')>=0.
    assert np.max(adverse_cov - positive_d_marg_variation) < 3e-8
    return {
        "C": C,
        "Dmarg": Dmarg,
        "weighted_energy": weighted_energy,
        "signed_cov": signed_cov,
        "adverse_cov": adverse_cov,
        "positive_cov": positive_cov,
        "positive_Dmarg_variation": positive_d_marg_variation,
        "negative_Dmarg_variation": negative_d_marg_variation,
        "min_Dmarg_derivative": min_d_marg_derivative,
        "min_C_derivative": min_c_derivative,
    }


def main() -> None:
    for name, A, beta, m in [
        ("A4", A4, 0.5, 3),
        ("A6", A6, 0.5, 3),
        ("A8", A8, 0.5, 4),
        ("A9", A9, 0.5, 3),
        ("A9", A9, 2.0, 7),
    ]:
        z = information_audit(A, beta, m)
        print(
            name,
            "n", z["n"], "m", z["m"], "beta", z["beta"],
            "D", f'{z["D"]:.10g}',
            "Dquot", f'{z["Dquot"]:.10g}',
            "C0", f'{z["C0"]:.10g}',
            "Csum", f'{z["Csum"]:.10g}',
            "Jsum", f'{z["Jsum"]:.10g}',
            "I", f'{z["Ifull"]:.10g}',
            "sumIcond", f'{z["Icondsum"]:.10g}',
        )
    exact = exact_orientation_counterexample()
    print("exact A6 beta=log(2)/2 m=3 orientation counterexample", exact)
    for name, A, beta, m in [
        ("A6", A6, 0.5, 3),
        ("A8", A8, 0.5, 4),
        ("A9", A9, 0.5, 3),
        ("A9", A9, 2.0, 3),
        ("A9", A9, 2.0, 7),
    ]:
        z = migration_audit(A, beta, m, order=48)
        C = np.asarray(z["C"])
        adverse = np.asarray(z["adverse_cov"])
        energy = np.asarray(z["weighted_energy"])
        nonzero = C > 1e-12
        print(
            name, "migration", "beta", beta, "m", m,
            "sumC", f"{float(np.sum(C)):.10g}",
            "sumE", f"{float(np.sum(energy)):.10g}",
            "sumSignedCov", f'{float(np.sum(z["signed_cov"])):.10g}',
            "sumAdverse", f"{float(np.sum(adverse)):.10g}",
            "maxEdgeE/C", f"{float(np.max(energy[nonzero] / C[nonzero])):.6g}",
            "maxEdgeAdverse/C", f"{float(np.max(adverse[nonzero] / C[nonzero])):.6g}",
            "minDmarginPrime", f'{float(np.min(z["min_Dmarg_derivative"])):.6g}',
            "minCPrime", f'{float(np.min(z["min_C_derivative"])):.6g}',
        )
    print("PASS harmonic_information_r35_check")


if __name__ == "__main__":
    main()
