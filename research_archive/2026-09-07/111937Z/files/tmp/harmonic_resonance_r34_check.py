#!/usr/bin/env python3
"""Wave 34 exact/numerical harmonic-resonance audits.

All scratch output stays below /home/math/quadra/tmp.  The script checks:
  * the exact vertex posterior-exclusion response;
  * the edgewise binary-KL integral and the moving-edge-mass correction;
  * orientation crossings and their energy contribution;
  * optimal and fixed-kappa restoring defects.
"""

from __future__ import annotations

import itertools
import math
import sys
from dataclasses import dataclass

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import logsumexp

sys.path.insert(0, "/home/math/quadra/tmp")
from endpoint_transport_r25 import (
    A4,
    A6,
    A8,
    A9,
    oriented_cuts,
    restriction_key,
)
from harmonic_composition_r31 import flip_vertex


def bernoulli_kl(p: float, q: float) -> float:
    ans = 0.0
    if p > 0.0:
        ans += p * math.log(p / q)
    if p < 1.0:
        ans += (1.0 - p) * math.log((1.0 - p) / (1.0 - q))
    return ans


def bernoulli_kl_logits(u1: float, u0: float) -> float:
    """D(Ber(logistic(u1)) || Ber(logistic(u0))), stably."""
    p1 = float(1.0 / (1.0 + math.exp(-u1))) if u1 >= 0.0 else float(math.exp(u1) / (1.0 + math.exp(u1)))
    return p1 * (u1 - u0) - float(np.logaddexp(0.0, u1) - np.logaddexp(0.0, u0))


def log_cosh(x: float) -> float:
    return float(np.logaddexp(x, -x) - math.log(2.0))


@dataclass
class Endpoint:
    A: np.ndarray
    beta: float
    m: int

    def __post_init__(self) -> None:
        self.A = np.asarray(self.A, dtype=np.int64)
        self.n = len(self.A)
        self.cuts = oriented_cuts(self.n)
        self.N = len(self.cuts)
        self.full = tuple(range(self.n))
        self.energies = np.einsum(
            "bij,ij->b", self.cuts, self.A, optimize=True
        ).astype(float)
        self.logZA = float(logsumexp(self.beta * self.energies))
        self.lognu = self.beta * self.energies - self.logZA
        self.nu = np.exp(self.lognu)
        self.selectors = list(itertools.combinations(range(self.n), self.m))
        self.selmask = np.zeros((len(self.selectors), self.n), dtype=bool)
        for si, S in enumerate(self.selectors):
            self.selmask[si, list(S)] = True

        Fall = np.empty((len(self.selectors), self.N), dtype=float)
        for si, S in enumerate(self.selectors):
            groups: dict[tuple[int, ...], list[int]] = {}
            for j, d in enumerate(self.cuts):
                groups.setdefault(restriction_key(d, S), []).append(j)
            assert len(groups) == 2**self.m
            for inds in groups.values():
                child = float(
                    np.sum(
                        self.A[np.ix_(S, S)]
                        * self.cuts[inds[0]][np.ix_(S, S)]
                    )
                )
                ext = self.energies[inds] - child
                val = float(logsumexp(self.beta * ext) - (self.n - self.m) * math.log(2.0))
                Fall[si, inds] = val
        self.Fall = Fall
        eF = np.exp(-Fall)
        z = eF @ self.nu
        self.q = z / np.sum(z)
        self.h = eF / z[:, None]
        self.f = self.q @ self.h
        assert abs(float(self.nu @ self.f) - 1.0) < 3e-10
        self.posterior = self.q[:, None] * self.h / self.f[None, :]
        self.g = np.log(self.f)

        idx = {
            restriction_key(d, self.full): j for j, d in enumerate(self.cuts)
        }
        neighbors: list[np.ndarray] = []
        neighbors.append(
            np.asarray([idx[restriction_key(-d, self.full)] for d in self.cuts])
        )
        for i in range(1, self.n):
            neighbors.append(
                np.asarray(
                    [
                        idx[restriction_key(flip_vertex(d, i), self.full)]
                        for d in self.cuts
                    ]
                )
            )
        self.neighbors = neighbors
        self.edges = [
            np.asarray([(j, int(jp)) for j, jp in enumerate(nei) if j < jp], dtype=int)
            for nei in neighbors
        ]

        # Posterior selector-exclusion probabilities, one row per vertex.
        self.excl = np.asarray(
            [np.sum(self.posterior[~self.selmask[:, i], :], axis=0) for i in range(self.n)]
        )

    def state(self, t: float) -> dict[str, object]:
        logw = self.lognu + t * self.g
        logw -= float(logsumexp(logw))
        mu = np.exp(logw)
        mean = float(mu @ self.g)
        X = self.g - mean
        var = float(mu @ X**2)
        Lg = np.zeros(self.N, dtype=float)
        energy_by_bit = np.zeros(self.n, dtype=float)
        cross_energy_by_bit = np.zeros(self.n, dtype=float)
        cov_correction_by_bit = np.zeros(self.n, dtype=float)
        conditional_kl_mu1_by_bit = np.zeros(self.n, dtype=float)
        conditional_kl_nuweight_by_bit = np.zeros(self.n, dtype=float)
        max_excl_err = 0.0
        max_cross_excl_violation = 0.0
        max_edge_formula_err = 0.0
        max_binary_kl_err = 0.0
        max_base_odds_err = 0.0
        crossing_mass_by_bit = np.zeros(self.n, dtype=float)

        for bit, (nei, edges) in enumerate(zip(self.neighbors, self.edges)):
            chi_dir = self.g[nei] - self.g
            omega_dir = self.lognu[nei] - self.lognu
            u_dir = omega_dir + t * chi_dir
            pjump = 1.0 / (1.0 + np.exp(-u_dir))
            Lg += pjump * chi_dir

            for x, y in edges:
                chi = float(self.g[y] - self.g[x])
                omega = float(self.lognu[y] - self.lognu[x])
                u = omega + t * chi
                M = float(mu[x] + mu[y])
                c = float(mu[x] * mu[y] / M)
                psi = chi * chi / (4.0 * math.cosh(u / 2.0) ** 2)
                max_edge_formula_err = max(max_edge_formula_err, abs(c * chi * chi - M * psi))
                energy_by_bit[bit] += c * chi * chi
                crossing = omega * (omega + chi) <= 0.0 and abs(chi) > 1e-14
                if crossing:
                    cross_energy_by_bit[bit] += c * chi * chi
                    crossing_mass_by_bit[bit] += M

                # Conditional binary path on this coordinate edge.
                d10 = bernoulli_kl_logits(omega + chi, omega)
                d10_explicit = (
                    0.5 * chi * math.tanh((omega + chi) / 2.0)
                    - log_cosh((omega + chi) / 2.0)
                    + log_cosh(omega / 2.0)
                )
                max_binary_kl_err = max(max_binary_kl_err, abs(d10 - d10_explicit))
                M1 = float(self.nu[x] * self.f[x] + self.nu[y] * self.f[y])
                M0 = float(self.nu[x] + self.nu[y])
                conditional_kl_mu1_by_bit[bit] += M1 * d10
                conditional_kl_nuweight_by_bit[bit] += M0 * d10

                # Exact derivative correction for moving edge marginals:
                # d/dt E_M[k_e(t)] = t M A'' + Cov_M(A',k).
                A0 = float(np.logaddexp(self.lognu[x], self.lognu[y]))
                At = float(np.logaddexp(self.lognu[x] + t * self.g[x], self.lognu[y] + t * self.g[y]) - A0)
                mt = float((mu[x] * self.g[x] + mu[y] * self.g[y]) / M)
                kt = t * mt - At
                # Accumulate ingredients; covariance is computed after edge loop.
                # Store through temporary lists local to this bit.
            # Recompute covariance over the edge marginal distribution cleanly.
            ms = []
            ks = []
            masses = []
            for x, y in edges:
                M = float(mu[x] + mu[y])
                mt = float((mu[x] * self.g[x] + mu[y] * self.g[y]) / M)
                A0 = float(np.logaddexp(self.lognu[x], self.lognu[y]))
                At = float(np.logaddexp(self.lognu[x] + t * self.g[x], self.lognu[y] + t * self.g[y]) - A0)
                masses.append(M)
                ms.append(mt)
                ks.append(t * mt - At)
            masses_a = np.asarray(masses)
            ms_a = np.asarray(ms)
            ks_a = np.asarray(ks)
            assert abs(float(np.sum(masses_a)) - 1.0) < 2e-10
            cov_correction_by_bit[bit] = float(
                np.sum(masses_a * (ms_a - np.sum(masses_a * ms_a)) * (ks_a - np.sum(masses_a * ks_a)))
            )

            if bit >= 1 and self.m < self.n:
                for x, y in edges:
                    chi = float(self.g[y] - self.g[x])
                    rhs = math.log(float(self.excl[bit, x] / self.excl[bit, y]))
                    max_excl_err = max(max_excl_err, abs(chi - rhs))
                    omega = float(self.lognu[y] - self.lognu[x])
                    if omega > 0.0 and omega + chi <= 0.0:
                        max_cross_excl_violation = max(
                            max_cross_excl_violation,
                            math.log(float(self.excl[bit, x])) + omega,
                        )
                    if omega < 0.0 and omega + chi >= 0.0:
                        max_cross_excl_violation = max(
                            max_cross_excl_violation,
                            math.log(float(self.excl[bit, y])) - omega,
                        )

            # Quadratic-signing base odds, including the distinct orientation bit.
            for x, y in edges:
                omega = float(self.lognu[y] - self.lognu[x])
                if bit == 0:
                    max_base_odds_err = max(
                        max_base_odds_err, abs(omega + 2.0 * self.beta * self.energies[x])
                    )
                else:
                    h_i = float(
                        sum(
                            self.A[bit, k] * self.cuts[x, bit, k]
                            for k in range(self.n)
                            if k != bit
                        )
                    )
                    max_base_odds_err = max(
                        max_base_odds_err, abs(omega + 4.0 * self.beta * h_i)
                    )

        B = -Lg
        E = float(np.sum(energy_by_bit))
        ibp = float(mu @ (X * B))
        W = float(mu @ B**2)
        assert abs(E - ibp) < 2e-9
        k_reg = E / var if var > 1e-30 else math.nan
        k_geom = math.sqrt(W / var) if var > 1e-30 else math.nan
        def ratio(kappa: float) -> float:
            return float(mu @ (B - kappa * X) ** 2) / E if E > 1e-30 else math.nan

        return {
            "mu": mu,
            "var": var,
            "energy": E,
            "energy_by_bit": energy_by_bit,
            "cross_energy_by_bit": cross_energy_by_bit,
            "crossing_mass_by_bit": crossing_mass_by_bit,
            "cov_correction_by_bit": cov_correction_by_bit,
            "conditional_kl_mu1_by_bit": conditional_kl_mu1_by_bit,
            "conditional_kl_nuweight_by_bit": conditional_kl_nuweight_by_bit,
            "B2": W,
            "C_scr": var / E if E > 1e-30 else math.nan,
            "k_reg": k_reg,
            "k_geom": k_geom,
            "defect1_E": ratio(1.0),
            "defect_kreg_E": ratio(k_reg),
            "defect_kgeom_E": ratio(k_geom),
            "max_excl_err": max_excl_err,
            "max_cross_excl_violation": max_cross_excl_violation,
            "max_edge_formula_err": max_edge_formula_err,
            "max_binary_kl_err": max_binary_kl_err,
            "max_base_odds_err": max_base_odds_err,
        }

    def integrated(self, order: int = 64) -> dict[str, object]:
        nodes, weights = leggauss(order)
        ts = (nodes + 1.0) / 2.0
        ws = weights / 2.0
        I = 0.0
        Iv = 0.0
        Ibits = np.zeros(self.n)
        Icross = np.zeros(self.n)
        Icov = np.zeros(self.n)
        max_def1 = 0.0
        min_kreg = math.inf
        max_C = 0.0
        for t, w in zip(ts, ws):
            z = self.state(float(t))
            I += float(w * t * z["energy"])
            Iv += float(w * t * z["var"])
            Ibits += w * t * np.asarray(z["energy_by_bit"])
            Icross += w * t * np.asarray(z["cross_energy_by_bit"])
            Icov += w * np.asarray(z["cov_correction_by_bit"])
            max_def1 = max(max_def1, float(z["defect1_E"]))
            min_kreg = min(min_kreg, float(z["k_reg"]))
            max_C = max(max_C, float(z["C_scr"]))
        mu1 = self.nu * self.f
        entropy = float(mu1 @ self.g)
        # For each coordinate: endpoint conditional KL = integral(t E_j + Cov_j).
        zend = self.state(1.0)
        condkl = np.asarray(zend["conditional_kl_mu1_by_bit"])
        correction_err = float(np.max(np.abs(condkl - (Ibits + Icov))))
        return {
            "energy_integral": I,
            "variance_integral": Iv,
            "entropy": entropy,
            "energy_bits": Ibits,
            "cross_bits": Icross,
            "cov_integrals": Icov,
            "condkl": condkl,
            "correction_err": correction_err,
            "max_defect1_E": max_def1,
            "min_kreg": min_kreg,
            "max_C": max_C,
        }


def run_case(name: str, A: np.ndarray, beta: float, m: int, order: int = 48) -> None:
    H = Endpoint(A, beta, m)
    z = H.integrated(order)
    orient = float(z["energy_bits"][0])
    orient_cross = float(z["cross_bits"][0])
    vertex = float(np.sum(z["energy_bits"][1:]))
    vertex_cross = float(np.sum(z["cross_bits"][1:]))
    print(
        name,
        "n", len(A),
        "beta", beta,
        "m", m,
        "Ent", f"{z['entropy']:.10g}",
        "I_E", f"{z['energy_integral']:.10g}",
        "I/Ent", f"{z['energy_integral']/z['entropy']:.6g}",
        "orient", f"{orient:.6g}",
        "orient_cross", f"{orient_cross:.6g}",
        "vertex", f"{vertex:.6g}",
        "vertex_cross", f"{vertex_cross:.6g}",
        "min_kreg", f"{z['min_kreg']:.6g}",
        "max_C", f"{z['max_C']:.6g}",
        "max_def1/E", f"{z['max_defect1_E']:.6g}",
        "Kcond", f"{float(np.sum(z['condkl'])):.6g}",
        "migration", f"{float(np.sum(z['cov_integrals'])):.6g}",
        "corr_err", f"{z['correction_err']:.3g}",
    )
    assert abs(float(z["variance_integral"]) - float(z["entropy"])) < 2e-8
    assert float(z["correction_err"]) < 2e-8
    s05 = H.state(0.5)
    assert float(s05["max_excl_err"]) < 2e-9
    assert float(s05["max_cross_excl_violation"]) < 2e-9
    assert float(s05["max_edge_formula_err"]) < 2e-9
    assert float(s05["max_binary_kl_err"]) < 2e-9
    assert float(s05["max_base_odds_err"]) < 2e-9


def main() -> None:
    # Exact minimizers recorded in the ledger.  Moderate grids keep the audit
    # fast enough for repeated independent checks.
    cases = [
        ("A4", A4, 0.5, 3),
        ("A6", A6, 0.1, 3),
        ("A8", A8, 0.5, 4),
        ("A9", A9, 0.5, 3),
        ("A9", A9, 2.0, 3),
        ("A9", A9, 0.5, 4),
        ("A9", A9, 2.0, 7),
    ]
    for args in cases:
        run_case(*args)
    print("PASS harmonic_resonance_r34_check")


if __name__ == "__main__":
    main()
