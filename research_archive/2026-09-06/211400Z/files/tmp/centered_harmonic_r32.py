#!/usr/bin/env python3
"""Exact finite audit of canonical harmonic centering and screening response."""

from __future__ import annotations

import itertools
import math
import sys
from collections import Counter
from fractions import Fraction

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from endpoint_transport_r25 import A4, A6, oriented_cuts, restriction_key
from harmonic_composition_r31 import HarmonicData, flip_vertex


def state_index(cuts: np.ndarray) -> dict[tuple[int, ...], int]:
    ids = tuple(range(len(cuts[0])))
    return {restriction_key(d, ids): j for j, d in enumerate(cuts)}


def audit(A: np.ndarray, beta: float, m: int) -> None:
    H = HarmonicData(A, beta)
    n = H.n
    selectors = list(itertools.combinations(range(n), m))
    Fall = np.asarray([[H.f(S, d) for d in H.cuts] for S in selectors])
    eF = np.exp(-Fall)

    # Canonical gauge: h_S has nu-mean one and q is the resulting barycenter
    # law.  The selector posterior q_S h_S/f is independent of this gauge.
    z = eF @ H.nu
    alpha = -np.log(z)
    q = z / np.sum(z)
    h = eF / z[:, None]
    f = q @ h
    assert abs(float(H.nu @ f) - 1.0) < 2e-12
    posterior = q[:, None] * h / f[None, :]
    direct_posterior = eF / np.sum(eF, axis=0)[None, :]
    assert np.allclose(posterior, direct_posterior, atol=2e-12, rtol=2e-12)
    Phi = -np.log(np.mean(eF, axis=0))
    assert np.ptp(Phi + np.log(f)) < 3e-12

    idx = state_index(H.cuts)
    full = tuple(range(n))
    neighbors: list[list[int]] = []
    # One orientation bit followed by n-1 independent projective vertex bits.
    neighbors.append([idx[restriction_key(-d, full)] for d in H.cuts])
    for i in range(1, n):
        neighbors.append([
            idx[restriction_key(flip_vertex(d, i), full)] for d in H.cuts
        ])

    print(f"n={n} beta={beta} m={m}")
    print("alpha_range", float(np.ptp(alpha)), "q_range", float(np.ptp(q)))
    for t in (0.0, 0.5, 1.0):
        mu = H.nu * f**t
        mu /= np.sum(mu)
        mean_phi = float(mu @ Phi)
        var_phi = float(mu @ (Phi - mean_phi) ** 2)
        energy = 0.0
        vertex_energy = 0.0
        orient_energy = 0.0
        common_sq = 0.0
        susceptibility_sq = 0.0
        post_var = 0.0
        max_identity_error = 0.0
        for bit, nei in enumerate(neighbors):
            for j, jp in enumerate(nei):
                if j >= jp:
                    continue
                conductance = float(mu[j] * mu[jp] / (mu[j] + mu[jp]))
                Z = np.log(h[:, jp]) - np.log(h[:, j])
                chi = math.log(float(posterior[:, j] @ np.exp(Z)))
                direct = math.log(float(f[jp] / f[j]))
                mean_z = float(posterior[:, j] @ Z)
                kl_back = chi - mean_z
                assert kl_back >= -2e-11
                post_prime = posterior[:, jp]
                kl_direct = float(np.sum(
                    posterior[:, j] * np.log(posterior[:, j] / post_prime)
                ))
                max_identity_error = max(max_identity_error,
                                         abs(chi - direct), abs(kl_back - kl_direct))
                term = conductance * chi * chi
                energy += term
                if bit == 0:
                    orient_energy += term
                else:
                    vertex_energy += term
                common_sq += conductance * mean_z * mean_z
                susceptibility_sq += conductance * kl_back * kl_back
                centered = Z - mean_z
                post_var += conductance * float(posterior[:, j] @ centered**2)
        ratio = var_phi / energy if energy > 0 else math.nan
        print("t", t, "varPhi", var_phi, "HBenergy", energy,
              "ratio", ratio, "vertex", vertex_energy, "orientation", orient_energy)
        print("  common_sq", common_sq, "KL_sq", susceptibility_sq,
              "posterior_var", post_var, "identity_error", max_identity_error)
    ent = float(H.nu @ (f * np.log(f)))
    print("entropy", ent, "oscPhi", float(np.ptp(Phi)))


def abstract_common_mode() -> None:
    # Any positive likelihood f on a two-point cut space can be duplicated
    # across selectors.  The selector posterior is then constant although the
    # parent likelihood and its entropy are nonconstant.
    nu = np.asarray([0.5, 0.5])
    f = np.asarray([1.5, 0.5])
    h = np.vstack([f, f])
    q = np.asarray([0.5, 0.5])
    post = q[:, None] * h / (q @ h)[None, :]
    assert np.allclose(post, 0.5)
    Z = np.log(h[:, 1]) - np.log(h[:, 0])
    chi = math.log(float(post[:, 0] @ np.exp(Z)))
    assert abs(chi - math.log(f[1] / f[0])) < 1e-12
    assert float(np.ptp(Z)) == 0.0
    ent = float(nu @ (f * np.log(f)))
    print("abstract_common_mode entropy", ent, "chi", chi,
          "posterior_score_variance", float(np.var(Z)))


def zero_temperature_a6_m3() -> None:
    """Exact max-plus data used for the centered-flatness wall."""
    H = HarmonicData(A6, 0.5)  # beta is irrelevant for the combinatorial data.
    selectors = list(itertools.combinations(range(6), 3))
    rows: list[tuple[int, int, Fraction]] = []
    profiles: dict[int, set[tuple[tuple[tuple[int, int], ...], ...]]] = {}
    for j, d in enumerate(H.cuts):
        heights: list[tuple[int, int]] = []
        profile: list[tuple[tuple[int, int], ...]] = []
        for S in selectors:
            inds = H.groups[S][restriction_key(d, S)]
            child = int(np.sum(A6[np.ix_(S, S)] * d[np.ix_(S, S)]))
            external = [int(H.energy[ii] - child) for ii in inds]
            profile.append(tuple(sorted(Counter(external).items())))
            height = max(external)
            heights.append((height, external.count(height)))
        lam = min(height for height, _ in heights)
        coeff = Fraction(2**3, len(selectors)) * sum(
            (Fraction(1, multiplicity)
             for height, multiplicity in heights if height == lam),
            Fraction(),
        )
        rows.append((int(H.energy[j]), lam, coeff))
        profiles.setdefault(int(H.energy[j]), set()).add(tuple(sorted(profile)))
    expected = Counter({
        (10, 4, Fraction(2, 3)): 12,
        (6, 4, Fraction(8, 15)): 20,
        (-6, 4, Fraction(2, 15)): 20,
        (-10, 8, Fraction(1, 1)): 12,
    })
    assert Counter(rows) == expected
    # U_beta is exactly constant inside each displayed energy class, not only
    # at leading exponential order.  The first omitted relative exponent in
    # the two leading classes is exp(-4 beta).
    assert set(profiles) == {10, 6, -6, -10}
    assert all(len(class_profiles) == 1 for class_profiles in profiles.values())
    entropy_coefficient = Fraction(5, 3) * (
        1 / 5 + (4 / 5) * math.log(4 / 5)
    )
    print("A6_m3_zero_temperature", expected)
    print("osc_Phi = 4 beta + log(2/3) + o(1)")
    print("entropy_coefficient_e^-4beta", float(entropy_coefficient))


def main() -> None:
    abstract_common_mode()
    zero_temperature_a6_m3()
    audit(A4, 0.5, 3)
    for m in (3, 4, 5):
        audit(A6, 0.5, m)
    print("PASS centered_harmonic_r32")


if __name__ == "__main__":
    main()
