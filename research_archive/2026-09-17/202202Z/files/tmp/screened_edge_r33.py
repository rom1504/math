#!/usr/bin/env python3
"""Exact edge-symmetrization and finite coercivity audits for Wave 33."""

from __future__ import annotations

import itertools
import math
import sys

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from endpoint_transport_r25 import A4, A6, restriction_key
from harmonic_composition_r31 import HarmonicData, flip_vertex


def kl(p: np.ndarray, q: np.ndarray) -> float:
    return float(np.sum(p * np.log(p / q)))


def edge_audit(A: np.ndarray, beta: float, m: int, t: float) -> dict[str, float]:
    H = HarmonicData(A, beta)
    Hneg = HarmonicData(A, -beta)
    n = H.n
    selectors = list(itertools.combinations(range(n), m))
    Fall = np.asarray([[H.f(S, d) for d in H.cuts] for S in selectors])
    eF = np.exp(-Fall)
    z = eF @ H.nu
    q = z / np.sum(z)
    h = eF / z[:, None]
    f = q @ h
    posterior = q[:, None] * h / f[None, :]
    g = np.log(f)

    full = tuple(range(n))
    idx = {restriction_key(d, full): j for j, d in enumerate(H.cuts)}
    neighbors: list[list[int]] = []
    # T_0 is the global orientation involution; T_i, i>=1, flips vertex i.
    neighbors.append([idx[restriction_key(-d, full)] for d in H.cuts])
    for i in range(1, n):
        neighbors.append([
            idx[restriction_key(flip_vertex(d, i), full)] for d in H.cuts
        ])

    mu = H.nu * f**t
    mu /= np.sum(mu)
    mean_g = float(mu @ g)
    var_g = float(mu @ (g - mean_g) ** 2)
    Lg = np.zeros_like(g)
    energy = 0.0
    orient_energy = 0.0
    max_err = 0.0
    max_orientation_beta_err = 0.0
    min_conductance = math.inf
    max_resonant_cost = 0.0

    for bit, nei in enumerate(neighbors):
        for j, jp in enumerate(nei):
            # Directed heat-bath drift.
            chi = float(g[jp] - g[j])
            omega = float(math.log(H.nu[jp] / H.nu[j]))
            u = omega + t * chi
            pjump = 1.0 / (1.0 + math.exp(-u))
            Lg[j] += pjump * chi

            if j >= jp:
                continue
            Z = np.log(h[:, jp] / h[:, j])
            chi_mgf = math.log(float(posterior[:, j] @ np.exp(Z)))
            k_minus = kl(posterior[:, j], posterior[:, jp])
            k_plus = kl(posterior[:, jp], posterior[:, j])
            mean_minus = float(posterior[:, j] @ Z)
            mean_plus = float(posterior[:, jp] @ Z)
            sym_chi = 0.5 * (mean_minus + mean_plus + k_minus - k_plus)
            jeff = mean_plus - mean_minus

            mass = float(mu[j] + mu[jp])
            conductance = float(mu[j] * mu[jp] / mass)
            conductance_cosh = mass / (4.0 * math.cosh(u / 2.0) ** 2)
            term = conductance * chi * chi
            energy += term
            if bit == 0:
                orient_energy += term
                # Reversing global cut orientation is beta -> -beta in F.
                for si, S in enumerate(selectors):
                    rhs = H.f(S, H.cuts[j]) - Hneg.f(S, H.cuts[j])
                    max_orientation_beta_err = max(
                        max_orientation_beta_err, abs(float(Z[si] - rhs))
                    )
            else:
                local_field = float(sum(A[bit, k] * H.cuts[j][bit, k]
                                        for k in range(n) if k != bit))
                max_err = max(max_err, abs(omega + 4.0 * beta * local_field))

            max_err = max(
                max_err,
                abs(chi - chi_mgf),
                abs(mean_minus - (chi - k_minus)),
                abs(mean_plus - (chi + k_plus)),
                abs(jeff - (k_minus + k_plus)),
                abs(chi - sym_chi),
                abs(math.log(mu[jp] / mu[j]) - u),
                abs(conductance - conductance_cosh),
            )
            if bit == 0:
                max_err = max(max_err, abs(omega + 2.0 * beta * H.energy[j]))
            min_conductance = min(min_conductance, conductance)
            if abs(u) <= 1.0:
                max_resonant_cost = max(max_resonant_cost, term / mass)

    # Reversibility/integration by parts for this one endpoint function.
    assert abs(float(mu @ Lg)) < 2e-11
    ibp = -float(mu @ ((g - mean_g) * Lg))
    max_err = max(max_err, abs(ibp - energy))
    restoring_drift = -Lg
    defect_one = float(mu @ (restoring_drift - (g - mean_g)) ** 2)
    drift_sq = float(mu @ restoring_drift**2)
    if energy > 0.0:
        delta_one = defect_one / energy
        comparison_bound = (
            (math.sqrt(delta_one) + math.sqrt(delta_one + 4.0)) / 2.0
        ) ** 2
        assert var_g / energy <= comparison_bound + 2e-11

    # A concrete endpoint-only restoring-drift certificate.  It need not hold.
    ratios = []
    sign_fail = 0
    for x, drift in zip(g - mean_g, -Lg):
        if abs(float(x)) <= 2e-12:
            continue
        ratios.append(float(drift / x))
        if float(x * drift) < -2e-12:
            sign_fail += 1

    return {
        "n": float(n),
        "beta": beta,
        "m": float(m),
        "t": t,
        "var": var_g,
        "energy": energy,
        "C_scr": var_g / energy if energy else math.nan,
        "orient_fraction": orient_energy / energy if energy else math.nan,
        "restore_min": min(ratios) if ratios else math.inf,
        "restore_max": max(ratios) if ratios else -math.inf,
        "restore_sign_fail": float(sign_fail),
        "defect_one_over_energy": defect_one / energy if energy else math.nan,
        "drift_sq_over_energy": drift_sq / energy if energy else math.nan,
        "min_conductance": min_conductance,
        "max_resonant_cost_per_mass": max_resonant_cost,
        "max_err": max_err,
        "orientation_beta_err": max_orientation_beta_err,
    }


def abstract_bottleneck(eps: float, a: float = 1.0) -> float:
    """One-function comparison can fail for an arbitrary endpoint likelihood."""
    # The base law is uniform on the square.  Multiplying f by a constant does
    # not change g differences or the t=1 law, so leave it unnormalized.
    f = np.asarray([1.0, eps, eps, math.exp(a)])
    mu = f / np.sum(f)
    g = np.log(f)
    mean = float(mu @ g)
    var = float(mu @ (g - mean) ** 2)
    edges = ((0, 1), (0, 2), (1, 3), (2, 3))
    energy = 0.0
    for x, y in edges:
        c = float(mu[x] * mu[y] / (mu[x] + mu[y]))
        energy += c * float(g[y] - g[x]) ** 2
    return var / energy


def main() -> None:
    bottlenecks = [abstract_bottleneck(eps) for eps in (1e-2, 1e-4, 1e-6)]
    assert bottlenecks[0] < bottlenecks[1] < bottlenecks[2]
    print("abstract_bottleneck_Cscr", *bottlenecks)
    rows = []
    for name, A in (("A4", A4), ("A6", A6)):
        ms = (3,) if len(A) == 4 else (3, 4, 5)
        for beta in (0.1, 0.5, 1.0, 2.0, 4.0):
            for m in ms:
                for t in (0.0, 0.5, 1.0):
                    out = edge_audit(A, beta, m, t)
                    assert out["max_err"] < 3e-9, (name, out)
                    assert out["orientation_beta_err"] < 3e-9, (name, out)
                    rows.append((name, out))
                    print(
                        name,
                        "beta", beta,
                        "m", m,
                        "t", t,
                        "C_scr", f"{out['C_scr']:.9g}",
                        "orient_frac", f"{out['orient_fraction']:.3g}",
                        "restore_min", f"{out['restore_min']:.9g}",
                        "sign_fail", int(out["restore_sign_fail"]),
                        "def1/E", f"{out['defect_one_over_energy']:.5g}",
                    )
    finite = [out for _, out in rows if math.isfinite(out["C_scr"])]
    print("C_scr_range", min(x["C_scr"] for x in finite),
          max(x["C_scr"] for x in finite))
    print("restore_min_global", min(x["restore_min"] for x in finite))
    print("max_identity_error", max(x["max_err"] for x in finite))
    print("max_orientation_beta_error",
          max(x["orientation_beta_err"] for x in finite))
    print("PASS screened_edge_r33")


if __name__ == "__main__":
    main()
