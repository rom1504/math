#!/usr/bin/env python3
"""Wave 43 audit of integrated context-mass harmonic migration bounds.

All files and numerical work stay below /home/math/quadra/tmp.  The checker
tests the exact covariance normalization, the statewise total-cost envelope,
the global (cross-coordinate) Cauchy bound, the separate orientation identity,
the A9 rare edge, and the abstract Wave-36 transient-mode-transfer wall.
"""

from __future__ import annotations

import itertools
import math
import sys
from dataclasses import dataclass

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import expit, logsumexp

sys.path.insert(0, "/home/math/quadra/tmp")
from endpoint_transport_r25 import A6, A8, A9
from harmonic_resonance_r34_check import Endpoint


@dataclass
class CubeSystem:
    lognu: np.ndarray
    g: np.ndarray
    edges: list[np.ndarray]
    exclusion: np.ndarray
    vertex_bits: tuple[int, ...]
    orientation_bit: int | None = None

    @property
    def N(self) -> int:
        return len(self.g)


def endpoint_system(H: Endpoint) -> CubeSystem:
    return CubeSystem(
        lognu=np.asarray(H.lognu, dtype=float),
        g=np.asarray(H.g, dtype=float),
        edges=H.edges,
        exclusion=np.asarray(H.excl, dtype=float),
        vertex_bits=tuple(range(1, H.n)),
        orientation_bit=0,
    )


def wave36_system(M: float) -> CubeSystem:
    """Reconstruct the strictly positive canonical Wave-36 migration wall."""
    noise = M ** -4
    orient = M ** -2
    states = list(itertools.product((0, 1), repeat=4))
    index = {x: k for k, x in enumerate(states)}
    nu = np.zeros(16)
    f = np.zeros(16)
    q = np.asarray([M / (M + 2.0), 0.5 / (M + 2.0), 1.5 / (M + 2.0)])
    hs = np.zeros((3, 16))
    for k, (o, u, v, y) in enumerate(states):
        nu[k] = 0.125 * ((1.0 - noise) if v == u else noise)
        ro = (1.0 - orient) if o == 0 else (1.0 + orient)
        hu = (1.0 / M) if u == 0 else (2.0 - 1.0 / M)
        hv = 1.0
        hy = (1.0 / 3.0) if y == 0 else (5.0 / 3.0)
        hs[:, k] = ro * np.asarray([hu, hv, hy])
        raw = ((1.0, 2.0), (M, M + 1.0))[u][y]
        f[k] = ro * 2.0 * raw / (M + 2.0)
    posterior = q[:, None] * hs / f[None, :]
    # Vertex coordinates are u,v,y = 1,2,3.  A singleton selector omits the
    # other two coordinates, hence these three omission probabilities.
    exclusion = np.ones((4, 16))
    exclusion[1] = posterior[1] + posterior[2]
    exclusion[2] = posterior[0] + posterior[2]
    exclusion[3] = posterior[0] + posterior[1]
    assert np.max(np.abs(np.sum(exclusion[1:], axis=0) - 2.0)) < 1e-11
    assert abs(float(np.sum(nu)) - 1.0) < 1e-12
    assert abs(float(nu @ f) - 1.0) < 1e-11
    edges: list[np.ndarray] = []
    for bit in range(4):
        pairs = []
        for k, state in enumerate(states):
            other = list(state)
            other[bit] ^= 1
            j = index[tuple(other)]
            if k < j:
                pairs.append((k, j))
        edges.append(np.asarray(pairs, dtype=int))
    return CubeSystem(
        lognu=np.log(nu),
        g=np.log(f),
        edges=edges,
        exclusion=exclusion,
        vertex_bits=(1, 2, 3),
        orientation_bit=0,
    )


def flat_fixed_density_system(nbits: int = 4, omission: float = 0.5) -> CubeSystem:
    """Uniform likelihood: raw C has a baseline although all migration is zero."""
    states = list(itertools.product((0, 1), repeat=nbits))
    index = {x: k for k, x in enumerate(states)}
    edges: list[np.ndarray] = []
    for bit in range(nbits):
        pairs = []
        for k, state in enumerate(states):
            other = list(state)
            other[bit] ^= 1
            j = index[tuple(other)]
            if k < j:
                pairs.append((k, j))
        edges.append(np.asarray(pairs, dtype=int))
    return CubeSystem(
        lognu=np.full(2**nbits, -nbits * math.log(2.0)),
        g=np.zeros(2**nbits),
        edges=edges,
        exclusion=np.full((nbits, 2**nbits), omission),
        vertex_bits=tuple(range(nbits)),
        orientation_bit=None,
    )


def law(system: CubeSystem, s: float) -> np.ndarray:
    z = system.lognu + s * system.g
    z -= float(logsumexp(z))
    return np.exp(z)


def binary_kl_logits(u: np.ndarray, u0: np.ndarray) -> np.ndarray:
    p = expit(u)
    ans = p * (u - u0) - (np.logaddexp(0.0, u) - np.logaddexp(0.0, u0))
    # Roundoff can create negative values of size about 1e-16 near s=0.
    return np.maximum(ans, 0.0)


def edge_costs(
    system: CubeSystem,
) -> tuple[list[np.ndarray], list[np.ndarray], np.ndarray, np.ndarray, float]:
    costs: list[np.ndarray] = [np.zeros(len(edges)) for edges in system.edges]
    scores: list[np.ndarray] = [np.zeros(len(edges)) for edges in system.edges]
    ctotal = np.zeros(system.N)
    scoretotal = np.zeros(system.N)
    exclusion_error = 0.0
    for bit in system.vertex_bits:
        edges = system.edges[bit]
        x, y = edges[:, 0], edges[:, 1]
        rx, ry = system.exclusion[bit, x], system.exclusion[bit, y]
        assert float(np.min(rx)) > 0.0 and float(np.max(rx)) <= 1.0 + 1e-11
        assert float(np.min(ry)) > 0.0 and float(np.max(ry)) <= 1.0 + 1e-11
        c = -0.5 * (np.log(rx) + np.log(ry))
        chi = system.g[y] - system.g[x]
        exclusion_error = max(
            exclusion_error,
            float(np.max(np.abs(chi - np.log(rx / ry)))),
        )
        assert float(np.min(c)) >= -2e-12
        assert float(np.max(np.abs(chi) - 2.0 * c)) <= 2e-10
        costs[bit] = c
        scores[bit] = np.abs(chi)
        ctotal[x] += c
        ctotal[y] += c
        scoretotal[x] += np.abs(chi)
        scoretotal[y] += np.abs(chi)
    return costs, scores, ctotal, scoretotal, exclusion_error


def node_data(
    system: CubeSystem,
    costs: list[np.ndarray],
    scores: list[np.ndarray],
    ctotal: np.ndarray,
    scoretotal: np.ndarray,
    s: float,
) -> dict[str, float | np.ndarray]:
    mu = law(system, s)
    mean = float(mu @ system.g)
    centered = system.g - mean
    variance = float(mu @ centered**2)
    load = np.zeros(system.N)
    energy_vertex = 0.0
    context_cov_sum = 0.0
    flux = 0.0
    score_flux = 0.0
    clipped_flux = 0.0
    exact_flux = 0.0
    moment_q = 0.0
    clippedtotal = np.zeros(system.N)
    max_mass_derivative_error = 0.0
    max_k_bound_error = 0.0

    for bit in system.vertex_bits:
        edges = system.edges[bit]
        x, y = edges[:, 0], edges[:, 1]
        M = mu[x] + mu[y]
        chi = system.g[y] - system.g[x]
        u0 = system.lognu[y] - system.lognu[x]
        us = u0 + s * chi
        p = expit(us)
        kval = binary_kl_logits(us, u0)
        load[x] += kval
        load[y] += kval
        energy_vertex += float(np.sum(M * p * (1.0 - p) * chi**2))

        X = (mu[x] * system.g[x] + mu[y] * system.g[y]) / M
        Mprime = M * (X - mean)
        direct_Mprime = mu[x] * centered[x] + mu[y] * centered[y]
        max_mass_derivative_error = max(
            max_mass_derivative_error,
            float(np.max(np.abs(Mprime - direct_Mprime))),
        )
        # Each M is a probability law on contexts: sum_e M_e=1.  Therefore
        # Cov_M(X,k)=sum_e M'_e k_e, with no extra normalization.
        context_cov_sum += float(np.sum(Mprime * kval))
        flux += float(np.sum(costs[bit] * np.maximum(-Mprime, 0.0)))
        score_flux += float(np.sum(scores[bit] * np.maximum(-Mprime, 0.0)))
        clipped = np.minimum(scores[bit], s * scores[bit] ** 2 / 8.0)
        clipped_flux += float(np.sum(clipped * np.maximum(-Mprime, 0.0)))
        exact_flux += float(np.sum(kval * np.maximum(-Mprime, 0.0)))
        clippedtotal[x] += clipped
        clippedtotal[y] += clipped
        bminus = float(np.sum(M * (mean - X) ** 2 * (X < mean)))
        acost = float(np.sum(M * costs[bit] ** 2))
        moment_q += math.sqrt(max(acost * bminus, 0.0))
        max_k_bound_error = max(
            max_k_bound_error,
            float(np.max(kval - s * np.abs(chi))),
            float(np.max(kval - 2.0 * s * costs[bit])),
        )

    load_mean = float(mu @ load)
    load_var = float(mu @ (load - load_mean) ** 2)
    global_cov = float(mu @ (centered * (load - load_mean)))
    assert abs(global_cov - context_cov_sum) < 2e-9

    orientation = {
        "energy": 0.0,
        "cost": 0.0,
        "cov": 0.0,
    }
    if system.orientation_bit is not None:
        bit = system.orientation_bit
        edges = system.edges[bit]
        x, y = edges[:, 0], edges[:, 1]
        M = mu[x] + mu[y]
        chi = system.g[y] - system.g[x]
        u0 = system.lognu[y] - system.lognu[x]
        us = u0 + s * chi
        p = expit(us)
        kval = binary_kl_logits(us, u0)
        X = (mu[x] * system.g[x] + mu[y] * system.g[y]) / M
        Mprime = M * (X - mean)
        orientation = {
            "energy": float(np.sum(M * p * (1.0 - p) * chi**2)),
            "cost": float(np.sum(M * kval)),
            "cov": float(np.sum(Mprime * kval)),
        }

    return {
        "mu": mu,
        "variance": variance,
        "load": load,
        "load_var": load_var,
        "global_cov": global_cov,
        "energy_vertex": energy_vertex,
        "flux": flux,
        "score_flux": score_flux,
        "clipped_flux": clipped_flux,
        "exact_flux": exact_flux,
        "moment_q": moment_q,
        "cost_square": float(mu @ ctotal**2),
        "score_square": float(mu @ scoretotal**2),
        "clipped_score_square": float(mu @ clippedtotal**2),
        "orientation_energy": orientation["energy"],
        "orientation_cost": orientation["cost"],
        "orientation_cov": orientation["cov"],
        "max_mass_derivative_error": max_mass_derivative_error,
        "max_k_bound_error": max_k_bound_error,
    }


def audit(system: CubeSystem, order: int = 128) -> dict[str, float]:
    costs, scores, ctotal, scoretotal, exclusion_error = edge_costs(system)
    nodes, weights = leggauss(order)
    entropy = 0.0
    covariance_integral = 0.0
    vertex_energy_integral = 0.0
    sharp_J2 = 0.0
    cost_pressure = 0.0
    score_pressure = 0.0
    clipped_score_pressure = 0.0
    exact_flux_bound = 0.0
    score_flux_bound = 0.0
    clipped_flux_bound = 0.0
    clipped_flux_J2 = 0.0
    direct_flux_bound = 0.0
    flux_J2 = 0.0
    moment_J2 = 0.0
    orientation_energy_integral = 0.0
    orientation_cov_integral = 0.0
    maximum_identity_error = 0.0
    maximum_k_bound_error = 0.0
    for node, weight0 in zip(nodes, weights):
        s = float((node + 1.0) / 2.0)
        weight = float(weight0 / 2.0)
        z = node_data(system, costs, scores, ctotal, scoretotal, s)
        var = float(z["variance"])
        entropy += weight * s * var
        covariance_integral += weight * float(z["global_cov"])
        vertex_energy_integral += weight * s * float(z["energy_vertex"])
        sharp_J2 += weight * float(z["load_var"]) / s
        cost_pressure += weight * s * float(z["cost_square"])
        score_pressure += weight * s * float(z["score_square"])
        clipped_score_pressure += weight * s * float(z["clipped_score_square"])
        exact_flux_bound += weight * float(z["exact_flux"])
        score_flux_bound += weight * s * float(z["score_flux"])
        clipped_flux_bound += weight * s * float(z["clipped_flux"])
        direct_flux_bound += weight * 2.0 * s * float(z["flux"])
        if var > 1e-30:
            flux_J2 += weight * s * float(z["flux"]) ** 2 / var
            clipped_flux_J2 += weight * s * float(z["clipped_flux"]) ** 2 / var
            moment_J2 += weight * s * float(z["moment_q"]) ** 2 / var
        orientation_energy_integral += weight * s * float(z["orientation_energy"])
        orientation_cov_integral += weight * float(z["orientation_cov"])
        maximum_identity_error = max(
            maximum_identity_error, float(z["max_mass_derivative_error"])
        )
        maximum_k_bound_error = max(
            maximum_k_bound_error, float(z["max_k_bound_error"])
        )

    endpoint = node_data(system, costs, scores, ctotal, scoretotal, 1.0)
    endpoint_vertex_cost = float(endpoint["mu"] @ endpoint["load"])
    endpoint_orientation_cost = float(endpoint["orientation_cost"])
    adverse_signed = max(-covariance_integral, 0.0)
    sharp_bound = math.sqrt(max(entropy * sharp_J2, 0.0))
    pressure_bound = 2.0 * math.sqrt(max(entropy * cost_pressure, 0.0))
    flux_cauchy_bound = 2.0 * math.sqrt(max(entropy * flux_J2, 0.0))
    clipped_flux_cauchy_bound = math.sqrt(max(entropy * clipped_flux_J2, 0.0))
    moment_cauchy_bound = 2.0 * math.sqrt(max(entropy * moment_J2, 0.0))

    # Exact endpoint integration-by-parts identities, vertex and orientation
    # separately.  The latter has no exclusion-cost formula.
    vertex_ibp_error = abs(
        vertex_energy_integral - (endpoint_vertex_cost - covariance_integral)
    )
    orientation_ibp_error = abs(
        orientation_energy_integral
        - (endpoint_orientation_cost - orientation_cov_integral)
    )
    assert maximum_identity_error < 3e-10
    assert maximum_k_bound_error < 3e-10
    assert vertex_ibp_error < 3e-9
    assert orientation_ibp_error < 3e-9
    assert sharp_J2 <= clipped_score_pressure + 2e-9
    assert clipped_score_pressure <= score_pressure + 2e-9
    assert sharp_J2 <= 4.0 * cost_pressure + 2e-9
    assert adverse_signed <= sharp_bound + 2e-9
    assert adverse_signed <= exact_flux_bound + 2e-9
    assert exact_flux_bound <= clipped_flux_bound + 2e-9
    assert clipped_flux_bound <= score_flux_bound + 2e-9
    assert clipped_flux_bound <= clipped_flux_cauchy_bound + 2e-9
    assert score_flux_bound <= direct_flux_bound + 2e-9
    assert adverse_signed <= direct_flux_bound + 2e-9
    assert direct_flux_bound <= flux_cauchy_bound + 2e-9
    assert flux_cauchy_bound <= moment_cauchy_bound + 2e-9
    assert sharp_bound <= pressure_bound + 2e-9
    return {
        "entropy": entropy,
        "adverse_signed": adverse_signed,
        "sharp_J2": sharp_J2,
        "sharp_bound": sharp_bound,
        "cost_pressure": cost_pressure,
        "score_pressure": score_pressure,
        "clipped_score_pressure": clipped_score_pressure,
        "pressure_bound": pressure_bound,
        "direct_flux_bound": direct_flux_bound,
        "exact_flux_bound": exact_flux_bound,
        "score_flux_bound": score_flux_bound,
        "clipped_flux_bound": clipped_flux_bound,
        "clipped_flux_J2": clipped_flux_J2,
        "flux_J2": flux_J2,
        "moment_J2": moment_J2,
        "endpoint_vertex_cost": endpoint_vertex_cost,
        "vertex_energy_integral": vertex_energy_integral,
        "vertex_ibp_error": vertex_ibp_error,
        "orientation_energy_integral": orientation_energy_integral,
        "endpoint_orientation_cost": endpoint_orientation_cost,
        "orientation_cov_integral": orientation_cov_integral,
        "orientation_ibp_error": orientation_ibp_error,
        "exclusion_error": exclusion_error,
        "max_context_derivative_error": maximum_identity_error,
        "max_k_bound_error": maximum_k_bound_error,
    }


def a9_factor_costs(beta: float = 2.0) -> dict[str, float]:
    """Three-factor costs on the exact A9,m=4 bad edge."""
    Hm = Endpoint(A9, beta, 4)
    Hp = Endpoint(A9, beta, 5)
    x, y, bit = 292, 308, 4
    bm = np.exp(-Hm.Fall)
    bp = np.exp(-Hp.Fall)
    map_m = {S: j for j, S in enumerate(Hm.selectors)}
    map_p = {S: j for j, S in enumerate(Hp.selectors)}
    Zm = np.sum(bm, axis=0)
    geom = []
    reveal = []
    zpx = 0.0
    zpy = 0.0
    numerator = 0.0
    for R in Hp.selectors:
        if bit not in R:
            continue
        S = tuple(j for j in R if j != bit)
        bs = float(bm[map_m[S], x])
        brx = float(bp[map_p[R], x])
        bry = float(bp[map_p[R], y])
        gg = math.sqrt(brx * bry)
        geom.append(gg)
        reveal.append(bs / gg)
        numerator += bs
        zpx += brx
        zpy += bry
    geom_a = np.asarray(geom)
    affinity = float(np.sum(geom_a) / math.sqrt(zpx * zpy))
    reveal_mean = float(np.sum(geom_a * np.asarray(reveal)) / np.sum(geom_a))
    gamma = math.sqrt(zpx * zpy / (float(Zm[x]) * float(Zm[y])))
    costs = {
        "C_reveal": -math.log(reveal_mean),
        "C_affinity": -math.log(affinity),
        "C_level": -math.log(gamma),
    }
    rx = numerator / float(Zm[x])
    ry = numerator / float(Zm[y])
    costs["C_total"] = -0.5 * math.log(rx * ry)
    costs["chi"] = float(Hm.g[y] - Hm.g[x])
    costs["rx"] = rx
    costs["ry"] = ry
    assert costs["C_level"] < 0.0
    assert costs["C_total"] >= abs(costs["chi"]) / 2.0 - 2e-10
    assert abs(
        costs["C_reveal"]
        + costs["C_affinity"]
        + costs["C_level"]
        - costs["C_total"]
    ) < 3e-9
    return costs


def a9_bad_edge_mass(beta: float, order: int = 128) -> dict[str, float]:
    H = Endpoint(A9, beta, 4)
    system = endpoint_system(H)
    costs, _, ctotal, _, _ = edge_costs(system)
    x, y, bit = 292, 308, 4
    where = np.flatnonzero(
        (system.edges[bit][:, 0] == x) & (system.edges[bit][:, 1] == y)
    )
    assert len(where) == 1
    C = float(costs[bit][int(where[0])])
    nodes, weights = leggauss(order)
    diagonal = 0.0
    state_pressure = 0.0
    max_mass = max(
        float(np.sum(law(system, endpoint)[[x, y]])) for endpoint in (0.0, 1.0)
    )
    for node, weight0 in zip(nodes, weights):
        s = float((node + 1.0) / 2.0)
        weight = float(weight0 / 2.0)
        mu = law(system, s)
        mass = float(mu[x] + mu[y])
        max_mass = max(max_mass, mass)
        diagonal += weight * s * mass * C * C
        state_pressure += weight * s * (
            mu[x] * ctotal[x] ** 2 + mu[y] * ctotal[y] ** 2
        )
    return {
        "beta": beta,
        "C_bad": C,
        "max_context_mass": max_mass,
        "diagonal_pressure": diagonal,
        "two_state_total_pressure": state_pressure,
    }


def compact(z: dict[str, float]) -> dict[str, float]:
    keys = (
        "entropy",
        "adverse_signed",
        "sharp_J2",
        "sharp_bound",
        "cost_pressure",
        "score_pressure",
        "clipped_score_pressure",
        "exact_flux_bound",
        "clipped_flux_bound",
        "pressure_bound",
        "endpoint_vertex_cost",
        "vertex_energy_integral",
        "orientation_energy_integral",
        "vertex_ibp_error",
        "orientation_ibp_error",
    )
    return {key: z[key] for key in keys}


def main() -> None:
    cases = (
        ("A6", A6, 0.5, 3),
        ("A8", A8, 2.0, 4),
        ("A9", A9, 2.0, 4),
        ("A9", A9, 2.0, 7),
    )
    for name, matrix, beta, m in cases:
        z = audit(endpoint_system(Endpoint(matrix, beta, m)), order=96)
        print(name, "beta", beta, "m", m, compact(z))
    print("A9 bad-edge factorization", a9_factor_costs())
    for beta in (1.0, 2.0, 4.0):
        print("A9 bad-edge mass", a9_bad_edge_mass(beta, order=96))
    flat = audit(flat_fixed_density_system(), order=48)
    print("flat fixed-density baseline wall", compact(flat))
    assert flat["sharp_J2"] == 0.0
    assert flat["score_pressure"] == 0.0
    assert flat["cost_pressure"] > 1.0
    wall_rows = []
    for M in (10.0, 100.0, 1000.0, 10000.0):
        z = audit(wave36_system(M), order=192)
        wall_rows.append((M, z))
        print("Wave36 wall M", int(M), compact(z))
    # Endpoint cost collapses much faster than the transient norm.  This is
    # the intended generic wall: rarity only at s=1 cannot prove the theorem.
    assert math.sqrt(wall_rows[-1][1]["sharp_J2"]) / wall_rows[-1][1]["endpoint_vertex_cost"] > 20.0
    assert wall_rows[-1][1]["adverse_signed"] / wall_rows[-1][1]["endpoint_vertex_cost"] > 10.0
    print("PASS harmonic_integrated_context_r43_check")


if __name__ == "__main__":
    main()
