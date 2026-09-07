#!/usr/bin/env python3
"""Wave 36 checks for a scalable abstract backtracking obstruction.

The example is a strictly positive canonical fixed-size selector mixture.
It has a separate (nonzero) global-orientation coordinate, three vertex
coordinates, and a correlated base law.  It is not claimed to arise from a
quadratic-signing minimizer.  Its purpose is to test what follows from the
endpoint-erasure/selector identities and fixed-size omission alone.
"""

from __future__ import annotations

import itertools
import math

import numpy as np
from scipy.integrate import quad


def kl(p: np.ndarray, q: np.ndarray) -> float:
    mask = p > 0.0
    return float(np.sum(p[mask] * np.log(p[mask] / q[mask])))


def marginal(p: np.ndarray, states: list[tuple[int, ...]], keep: tuple[int, ...]) -> np.ndarray:
    ans = np.zeros(2 ** len(keep))
    for mass, x in zip(p, states):
        key = sum(x[j] << k for k, j in enumerate(keep))
        ans[key] += mass
    return ans


def binary_kl(p: float, q: float) -> float:
    return p * math.log(p / q) + (1.0 - p) * math.log((1.0 - p) / (1.0 - q))


def build_example(M: float, noise: float | None = None, orient: float | None = None) -> dict:
    # Coordinates: O is global orientation; u,v,y are vertex coordinates.
    # The base makes v a noisy redundant copy of u.  The endpoint likelihood
    # depends on (O,u,y), not v.
    if noise is None:
        noise = M ** -4
    if orient is None:
        orient = M ** -2
    states = list(itertools.product((0, 1), repeat=4))
    nu = np.zeros(len(states))
    f = np.zeros(len(states))

    # Fixed-size singleton selector components on {u,v,y}.  Each component
    # is normalized and is invariant under every omitted vertex coordinate.
    q = np.asarray([M / (M + 2.0), 0.5 / (M + 2.0), 1.5 / (M + 2.0)])
    hs = np.zeros((3, len(states)))
    for k, (o, u, v, y) in enumerate(states):
        nu[k] = 0.125 * ((1.0 - noise) if v == u else noise)
        r_o = (1.0 - orient) if o == 0 else (1.0 + orient)
        h_u = (1.0 / M) if u == 0 else (2.0 - 1.0 / M)
        h_v = 1.0
        h_y = (1.0 / 3.0) if y == 0 else (5.0 / 3.0)
        hs[:, k] = r_o * np.asarray([h_u, h_v, h_y])
        raw = ((1.0, 2.0), (M, M + 1.0))[u][y]
        f[k] = r_o * 2.0 * raw / (M + 2.0)

    assert abs(float(np.sum(nu)) - 1.0) < 1e-13
    assert np.max(np.abs(hs @ nu - 1.0)) < 1e-11
    assert abs(float(np.sum(q)) - 1.0) < 1e-13
    assert np.max(np.abs(q @ hs - f)) < 1e-12
    assert abs(float(nu @ f) - 1.0) < 1e-12
    assert np.min(nu) > 0 and np.min(hs) > 0 and np.min(f) > 0

    mu = nu * f
    Dfull = kl(mu, nu)
    coords = range(4)
    C = np.zeros(4)
    for j in coords:
        keep = tuple(i for i in coords if i != j)
        C[j] = Dfull - kl(marginal(mu, states, keep), marginal(nu, states, keep))
    assert np.min(C) > -2e-11

    # Lift P(S,d)=q(S)nu(d)h_S(d), verify J=C+I coordinatewise.
    component = nu[None, :] * hs
    joint_D = float(sum(q[s] * kl(component[s], nu) for s in range(3)))
    Ifull = joint_D - Dfull
    J = np.zeros(4)
    Icond = np.zeros(4)
    for j in coords:
        keep = tuple(i for i in coords if i != j)
        nu_minus = marginal(nu, states, keep)
        mu_minus = marginal(mu, states, keep)
        comp_minus = np.asarray([marginal(component[s], states, keep) for s in range(3)])
        J[j] = sum(q[s] * (kl(component[s], nu) - kl(comp_minus[s], nu_minus)) for s in range(3))
        Iminus = 0.0
        for s in range(3):
            mask = comp_minus[s] > 0
            Iminus += q[s] * float(np.sum(comp_minus[s, mask] * np.log(comp_minus[s, mask] / mu_minus[mask])))
        Icond[j] = Ifull - Iminus
    assert np.max(np.abs(C - (J - Icond))) < 2e-9

    # Exact singleton posterior: the three vertex omission probabilities sum
    # to n-m=2 at every state.  Selector order is u,v,y.
    posterior = q[:, None] * hs / f[None, :]
    assert np.max(np.abs(np.sum(posterior, axis=0) - 1.0)) < 1e-12
    omissions = np.vstack((posterior[1] + posterior[2],
                           posterior[0] + posterior[2],
                           posterior[0] + posterior[1]))
    assert np.max(np.abs(np.sum(omissions, axis=0) - 2.0)) < 1e-12

    # J_(S,j)=0 whenever the singleton S omits vertex j (vertex coordinates
    # u,v,y are chart coordinates 1,2,3).
    selector_vertex = (1, 2, 3)
    for s, active_j in enumerate(selector_vertex):
        for j in selector_vertex:
            if j == active_j:
                continue
            keep = tuple(i for i in coords if i != j)
            nu_minus = marginal(nu, states, keep)
            comp_minus = marginal(component[s], states, keep)
            assert abs(kl(component[s], nu) - kl(comp_minus, nu_minus)) < 2e-10

    g = np.log(f)
    state_index = {x: k for k, x in enumerate(states)}

    def mut(x: tuple[int, ...], j: int) -> tuple[int, ...]:
        z = list(x)
        z[j] ^= 1
        return tuple(z)

    edges_by_j: list[list[tuple[int, int]]] = []
    for j in coords:
        edges = []
        for k, x in enumerate(states):
            z = state_index[mut(x, j)]
            if k < z:
                edges.append((k, z))
        edges_by_j.append(edges)

    # Endpoint costs from conditional binary KL agree with erasure KL.
    edge_K: list[np.ndarray] = []
    for j, edges in enumerate(edges_by_j):
        vals = []
        csum = 0.0
        for x, y in edges:
            p0 = nu[y] / (nu[x] + nu[y])
            p1 = mu[y] / (mu[x] + mu[y])
            kval = binary_kl(float(p1), float(p0))
            vals.append(kval)
            csum += (mu[x] + mu[y]) * kval
        edge_K.append(np.asarray(vals))
        assert abs(csum - C[j]) < 2e-9

    def path(t: float) -> tuple[np.ndarray, float]:
        w = nu * np.exp(t * g)
        mt = w / np.sum(w)
        return mt, float(mt @ g)

    backward_by_j = np.zeros(4)
    adverse_by_j = np.zeros(4)
    for j, edges in enumerate(edges_by_j):
        for ei, (x, y) in enumerate(edges):
            def neg_derivative(t: float) -> float:
                mt, mean = path(t)
                deriv = mt[x] * (g[x] - mean) + mt[y] * (g[y] - mean)
                return max(-float(deriv), 0.0)

            variation, _ = quad(neg_derivative, 0.0, 1.0, epsabs=2e-12, epsrel=2e-10, limit=200)
            backward_by_j[j] += edge_K[j][ei] * variation

        def adverse_covariance(t: float) -> float:
            mt, mean = path(t)
            covariance = 0.0
            for x, y in edges:
                mass_derivative = mt[x] * (g[x] - mean) + mt[y] * (g[y] - mean)
                p0 = nu[y] / (nu[x] + nu[y])
                wt_x = nu[x] * math.exp(t * g[x])
                wt_y = nu[y] * math.exp(t * g[y])
                pt = wt_y / (wt_x + wt_y)
                kt = binary_kl(float(pt), float(p0))
                covariance += mass_derivative * kt
            return max(-float(covariance), 0.0)

        adverse_by_j[j], _ = quad(adverse_covariance, 0.0, 1.0,
                                   epsabs=2e-12, epsrel=2e-10, limit=200)

    k0 = binary_kl(2.0 / 3.0, 0.5)
    kM = binary_kl((M + 1.0) / (2.0 * M + 1.0), 0.5)
    low_endpoint = 3.0 / (2.0 * M + 4.0)
    cy_formula = low_endpoint * k0 + (1.0 - low_endpoint) * kM
    by_lower = k0 * (0.5 - low_endpoint)
    assert abs(C[3] - cy_formula) < 2e-10
    assert backward_by_j[3] + 2e-9 >= by_lower
    assert abs(C[2]) < 2e-9  # v is an omitted redundant copy channel.
    assert C[0] > 0.0  # the separate global-orientation charge is retained.
    assert C[1] <= M ** -2

    return {
        "M": M,
        "noise": noise,
        "orient": orient,
        "q": q,
        "C": C,
        "Csum": float(np.sum(C)),
        "J": J,
        "Icond": Icond,
        "backward": backward_by_j,
        "Bsum": float(np.sum(backward_by_j)),
        "ratio": float(np.sum(backward_by_j) / np.sum(C)),
        "adverse": adverse_by_j,
        "Asum": float(np.sum(adverse_by_j)),
        "Aratio": float(np.sum(adverse_by_j) / np.sum(C)),
        "k0": k0,
        "kM": kM,
        "Cy_formula": cy_formula,
        "By_lower": by_lower,
    }


def main() -> None:
    last = 0.0
    for M in (10.0, 30.0, 100.0, 300.0, 1000.0):
        z = build_example(M)
        print(
            "M", int(M),
            "C", np.array2string(z["C"], precision=5),
            "J_y", f'{z["J"][3]:.8g}',
            "I_y", f'{z["Icond"][3]:.8g}',
            "B", np.array2string(z["backward"], precision=5),
            "B/C", f'{z["ratio"]:.8g}',
            "(B/C)/M", f'{z["ratio"] / M:.8g}',
            "A/C", f'{z["Aratio"]:.8g}',
        )
        assert z["ratio"] > last
        assert z["ratio"] / M > 0.25
        assert z["adverse"][3] * math.log(M) ** 2 > 0.02
        assert z["Aratio"] * math.log(M) ** 2 / M > 0.25
        last = z["ratio"]
    assert last > 100.0
    print("PASS harmonic_erasure_r36_check")


if __name__ == "__main__":
    main()
