#!/usr/bin/env python3
"""Wave 44 checks: intrinsic-load flips, deletion Hessian, and A8 wall."""

from __future__ import annotations

import itertools
import math
import sys
from collections import Counter, defaultdict
from fractions import Fraction

import numpy as np
from numpy.polynomial.legendre import leggauss

sys.path.insert(0, "/home/math/quadra/tmp")
from endpoint_transport_r25 import A6, A8, A9, restriction_key
from harmonic_resonance_r34_check import Endpoint
from harmonic_integrated_context_r43_check import (
    audit,
    edge_costs,
    endpoint_system,
    law,
    node_data,
    wave36_system,
)


def directed_kv(H: Endpoint, s: float) -> tuple[np.ndarray, np.ndarray]:
    """Statewise conditional KL and curvature, one row per chart bit."""
    k = np.zeros((H.n, H.N))
    v = np.zeros((H.n, H.N))
    for bit in range(H.n):
        nei = H.neighbors[bit]
        omega = H.lognu[nei] - H.lognu
        chi = H.g[nei] - H.g
        u = omega + s * chi
        p = 1.0 / (1.0 + np.exp(-np.clip(u, -700.0, 700.0)))
        kval = p * (u - omega) - (
            np.logaddexp(0.0, u) - np.logaddexp(0.0, omega)
        )
        k[bit] = np.maximum(kval, 0.0)
        v[bit] = p * (1.0 - p) * chi**2
    return k, v


def flip_audit(H: Endpoint, s: float = 0.63, order: int = 64) -> dict[str, float]:
    k, _ = directed_kv(H, s)
    load = np.sum(k[1:], axis=0)
    flip_error = 0.0
    base_square_error = 0.0
    score_square_error = 0.0
    curvature_integral_error = 0.0
    nodes, weights = leggauss(order)
    ts = s * (nodes + 1.0) / 2.0
    ws = s * weights / 2.0
    integral_k = np.zeros_like(k)
    for t, w in zip(ts, ws):
        _, v = directed_kv(H, float(t))
        integral_k += float(w * t) * v
    curvature_integral_error = float(np.max(np.abs(k - integral_k)))

    for j in range(1, H.n):
        dj = H.neighbors[j]
        rhs = np.sum(k[1:, dj] - k[1:], axis=0) - (k[j, dj] - k[j])
        flip_error = max(flip_error, float(np.max(np.abs(load[dj] - load - rhs))))
        for i in range(1, H.n):
            if i == j:
                continue
            di = H.neighbors[i]
            omega_i = H.lognu[di] - H.lognu
            chi_i = H.g[di] - H.g
            base_increment = omega_i[dj] - omega_i
            predicted = np.asarray(
                [8.0 * H.beta * H.A[i, j] * H.cuts[d, i, j] for d in range(H.N)]
            )
            base_square_error = max(
                base_square_error, float(np.max(np.abs(base_increment - predicted)))
            )
            mixed_ij = chi_i[dj] - chi_i
            chi_j = H.g[H.neighbors[j]] - H.g
            mixed_ji = chi_j[di] - chi_j
            score_square_error = max(
                score_square_error, float(np.max(np.abs(mixed_ij - mixed_ji)))
            )
    assert curvature_integral_error < 2e-10
    assert flip_error < 2e-12
    assert base_square_error < 2e-11
    assert score_square_error < 2e-12
    return {
        "curvature_integral_error": curvature_integral_error,
        "flip_increment_error": flip_error,
        "quadratic_base_increment_error": base_square_error,
        "mixed_score_symmetry_error": score_square_error,
    }


def marginal(p: np.ndarray, n: int, keep: tuple[int, ...]) -> np.ndarray:
    ans = np.zeros(2 ** len(keep))
    for mass, state in zip(p, itertools.product((0, 1), repeat=n)):
        key = sum(state[j] << q for q, j in enumerate(keep))
        ans[key] += mass
    return ans


def kl(p: np.ndarray, q: np.ndarray) -> float:
    mask = p > 0.0
    return float(np.sum(p[mask] * np.log(p[mask] / q[mask])))


def deletion_audit(H: Endpoint, s: float = 1.0) -> dict[str, float | int]:
    system = endpoint_system(H)
    mu, nu = law(system, s), np.exp(system.lognu)
    full = tuple(range(H.n))
    F: dict[tuple[int, ...], float] = {}
    for r in range(3):
        for erased in itertools.combinations(full, r):
            keep = tuple(i for i in full if i not in erased)
            F[erased] = kl(marginal(mu, H.n, keep), marginal(nu, H.n, keep))
    C = {i: F[()] - F[(i,)] for i in range(H.n)}
    mixed = {
        (i, j): F[()] - F[(i,)] - F[(j,)] + F[(i, j)]
        for i, j in itertools.combinations(range(H.n), 2)
    }
    vertex = range(1, H.n)
    deletion_increments = []
    identity_error = 0.0
    for j in vertex:
        predicted = C[j] + sum(mixed[tuple(sorted((i, j)))] for i in vertex if i != j)
        erased_j = F[(j,)]
        coarse_load = sum(
            erased_j - F[tuple(sorted((i, j)))] for i in vertex if i != j
        )
        direct = sum(C[i] for i in vertex) - coarse_load
        identity_error = max(identity_error, abs(predicted - direct))
        deletion_increments.append(direct)
    vals = [mixed[(i, j)] for i, j in itertools.combinations(vertex, 2)]
    assert identity_error < 2e-12
    return {
        "mixed_min": min(vals),
        "mixed_max": max(vals),
        "mixed_negative": sum(x < -1e-12 for x in vals),
        "mixed_positive": sum(x > 1e-12 for x in vals),
        "deletion_min": min(deletion_increments),
        "deletion_max": max(deletion_increments),
        "deletion_identity_error": identity_error,
    }


def load_dirichlet(H: Endpoint, s: float) -> tuple[float, float]:
    system = endpoint_system(H)
    costs, scores, ctotal, scoretotal, _ = edge_costs(system)
    z = node_data(system, costs, scores, ctotal, scoretotal, s)
    mu, load = np.asarray(z["mu"]), np.asarray(z["load"])
    energy = 0.0
    for edges in H.edges:  # include orientation; it contributes zero in the wall
        x, y = edges[:, 0], edges[:, 1]
        M = mu[x] + mu[y]
        energy += float(np.sum(mu[x] * mu[y] / M * (load[y] - load[x]) ** 2))
    return float(z["load_var"]), energy


def integrated_dirichlet(H: Endpoint, order: int = 160) -> tuple[float, float]:
    nodes, weights = leggauss(order)
    J2 = D2 = 0.0
    for node, weight0 in zip(nodes, weights):
        s = float((node + 1.0) / 2.0)
        var, energy = load_dirichlet(H, s)
        weight = float(weight0 / 2.0)
        J2 += weight * var / s
        D2 += weight * energy / s
    return J2, D2


def tropical_data(H: Endpoint) -> tuple[np.ndarray, list[Fraction]]:
    exponent = np.full(H.N, 10**9, dtype=int)
    coefficient = [Fraction() for _ in range(H.N)]
    for S in H.selectors:
        groups: dict[tuple[int, ...], list[int]] = defaultdict(list)
        for j, d in enumerate(H.cuts):
            groups[restriction_key(d, S)].append(j)
        for inds in groups.values():
            first = inds[0]
            child = int(np.sum(H.A[np.ix_(S, S)] * H.cuts[first][np.ix_(S, S)]))
            external = [int(H.energies[j] - child) for j in inds]
            maximum, multiplicity = max(external), external.count(max(external))
            for j in inds:
                if maximum < exponent[j]:
                    exponent[j] = maximum
                    coefficient[j] = Fraction(1, multiplicity)
                elif maximum == exponent[j]:
                    coefficient[j] += Fraction(1, multiplicity)
    return exponent, coefficient


def a8_tropical_wall() -> dict[str, object]:
    H = Endpoint(A8, 1.0, 4)
    exponent, coeff = tropical_data(H)
    energy = H.energies.astype(int)
    line = energy - exponent
    top = int(np.max(line))
    active = np.flatnonzero(line == top)
    load_slope = [Fraction() for _ in range(H.N)]
    for bit in range(1, H.n):
        for x, y in H.edges[bit]:
            if energy[x] == energy[y]:
                value = Fraction()
            else:
                low = x if energy[x] < energy[y] else y
                if line[x] > line[y]:
                    px, py = Fraction(1), Fraction()
                elif line[y] > line[x]:
                    px, py = Fraction(), Fraction(1)
                else:
                    px = coeff[x] / (coeff[x] + coeff[y])
                    py = coeff[y] / (coeff[x] + coeff[y])
                value = abs(int(energy[x] - energy[y])) * (px if low == x else py)
            load_slope[x] += value
            load_slope[y] += value
    histogram = Counter((int(energy[d]), str(coeff[d]), str(load_slope[d])) for d in active)
    differing_edges = []
    for bit in range(H.n):
        for x, y in H.edges[bit]:
            if load_slope[x] != load_slope[y]:
                differing_edges.append(top - min(int(line[x]), int(line[y])))
    assert top == 12 and len(active) == 20
    assert histogram == Counter({
        (16, "1/6", "4/7"): 7,
        (20, "1", "4/7"): 7,
        (16, "1/6", "0"): 5,
        (20, "1", "0"): 1,
    })
    assert min(differing_edges) == 4
    variance_coefficient = Fraction(11, 225)
    J_coefficient = 7.0 * math.log(7.0 / 6.0) + (217.0 / 36.0) * math.log(4.0 / 5.0) + 49.0 / 180.0
    return {
        "endpoint_line": top,
        "active_states": len(active),
        "active_histogram": dict(histogram),
        "minimum_different-load_edge_gap": min(differing_edges),
        "Var_L_over_beta2_limit": str(variance_coefficient),
        "J_L2_over_beta_limit": J_coefficient,
        "endpoint_cost_over_beta_limit": "7/15",
        "J2_over_endpoint_limit": J_coefficient * 15.0 / 7.0,
        "Dirichlet_ratio": "Theta(exp(4 beta))",
    }


def main() -> None:
    cases = (
        ("A6", A6, 0.5, 3),
        ("A8", A8, 2.0, 4),
        ("A9", A9, 2.0, 4),
        ("A9", A9, 2.0, 7),
    )
    for name, A, beta, m in cases:
        H = Endpoint(A, beta, m)
        print(name, beta, m, "flip", flip_audit(H), "deletion", deletion_audit(H))
    print("A8 tropical wall", a8_tropical_wall())
    for beta in (1.0, 2.0, 4.0, 8.0):
        H = Endpoint(A8, beta, 4)
        var, energy = load_dirichlet(H, 1.0)
        J2, D2 = integrated_dirichlet(H, order=128)
        endpoint = audit(endpoint_system(H), order=128)["endpoint_vertex_cost"]
        print("A8 beta", beta, {
            "pointwise_Var_over_Dirichlet": var / energy,
            "J2": J2,
            "integrated_Dirichlet": D2,
            "integrated_ratio": J2 / D2,
            "J2_over_endpoint": J2 / endpoint,
        })
    for M in (100.0, 1000.0, 10000.0):
        z = audit(wave36_system(M), order=192)
        print("Wave36 M", int(M), {
            "J2_over_endpoint": z["sharp_J2"] / z["endpoint_vertex_cost"],
            "adverse_over_endpoint": z["adverse_signed"] / z["endpoint_vertex_cost"],
        })
    print("PASS harmonic_load_increments_r44_check")


if __name__ == "__main__":
    main()
