#!/usr/bin/env python3
"""Wave 27 audit of exact signed deletion-deck cavity identities."""

from __future__ import annotations

import math

import numpy as np

from check_finite_bridge_r16 import A5
from check_response_dual_r16 import A8, A9, spins
from verify_compatible_replacement_r12 import A6


def oriented_states(A: np.ndarray):
    """Oriented projective states: 2^n states for an order-n signing."""
    X0 = spins(len(A))
    X = np.tile(X0, (2, 1))
    sigma = np.repeat(np.asarray((-1, 1), dtype=np.int64), len(X0))
    energy = sigma * np.einsum("bi,ij,bj->b", X, A, X)
    return X, sigma, energy.astype(float)


def qnorm(A: np.ndarray) -> int:
    return int(np.max(oriented_states(A)[2]))


def normalized_partition(A: np.ndarray, beta: float) -> float:
    energy = oriented_states(A)[2]
    return float(np.exp(beta * energy).mean())


def cycle_polynomial(A: np.ndarray, beta: float) -> float:
    n = len(A)
    return normalized_partition(A, beta) / math.cosh(2 * beta) ** (n * (n - 1) // 2)


def child_data(A: np.ndarray, i: int, beta: float):
    n = len(A)
    keep = [j for j in range(n) if j != i]
    C = A[np.ix_(keep, keep)]
    row = A[i, keep]
    Y, sigma, energy = oriented_states(C)
    weights = np.exp(beta * energy - np.max(beta * energy))
    weights /= weights.sum()
    fields = Y @ row
    v = float(weights @ fields**2)

    q_parent = qnorm(A)
    q_child = qnorm(C)
    d = q_parent - q_child
    deficit = q_child - energy
    mean_deficit = float(weights @ deficit)
    second_deficit = float(weights @ deficit**2)

    # The two parent-extension deficits from a child state.  The sigma in
    # the oriented field only swaps the two values, so the product is as below.
    extension_product = (d + deficit) ** 2 - 4 * fields**2
    assert np.min(extension_product) >= -1e-10
    extension_product_mean = float(weights @ extension_product)
    assert math.isclose(
        4 * v,
        d * d + 2 * d * mean_deficit + second_deficit - extension_product_mean,
        rel_tol=2e-11,
        abs_tol=2e-11,
    )

    # Partition derivative formula for the child deficit law.
    logD_prime = -mean_deficit
    logD_second = second_deficit - mean_deficit**2
    derivative_envelope = (d - logD_prime) ** 2 + logD_second
    assert math.isclose(
        derivative_envelope,
        d * d + 2 * d * mean_deficit + second_deficit,
        rel_tol=2e-11,
        abs_tol=2e-11,
    )
    assert 4 * v <= derivative_envelope + 1e-10

    # A child-only upper bound for the actual cavity reward.  It is the
    # pointwise extension cap 2|L| <= d+delta integrated without discarding
    # the signed child Gibbs measure.
    D0 = float(len(deficit))
    Dbeta = float(np.exp(-beta * deficit).sum())
    D2beta = float(np.exp(-2 * beta * deficit).sum())
    profile_upper = math.log(
        (math.exp(beta * d) * D0 + math.exp(-beta * d) * D2beta)
        / (2 * Dbeta)
    ) / beta
    kappa = math.log(float(weights @ np.cosh(2 * beta * fields))) / beta
    assert kappa <= profile_upper + 1e-10

    return {
        "C": C,
        "row": row,
        "Y": Y,
        "sigma": sigma,
        "weights": weights,
        "fields": fields,
        "v": v,
        "d": d,
        "mean_deficit": mean_deficit,
        "second_deficit": second_deficit,
        "extension_product_mean": extension_product_mean,
        "derivative_envelope": derivative_envelope,
        "kappa": kappa,
        "profile_upper": profile_upper,
    }


def star_twist(A: np.ndarray, i: int, mask: int) -> np.ndarray:
    out = A.copy()
    keep = [j for j in range(len(A)) if j != i]
    for bit, j in enumerate(keep):
        if (mask >> bit) & 1:
            out[i, j] *= -1
            out[j, i] *= -1
    return out


def audit_matrix(A: np.ndarray, name: str, beta: float) -> None:
    A = np.asarray(A, dtype=np.int64)
    r = len(A)
    h = r - 1
    rho = math.tanh(2 * beta)
    q = qnorm(A)
    alpha = q / r**1.5
    threshold = 9 * alpha * alpha * r / 16

    parent_X, parent_sigma, parent_energy = oriented_states(A)
    parent_weights = np.exp(beta * parent_energy - np.max(beta * parent_energy))
    parent_weights /= parent_weights.sum()

    deletion_polynomials = []
    cavity_moments = []
    signed_degree_two = []
    competitor_floor_bounds = []
    full_q_floor_bounds = []
    floor_ratios = []
    derivative_envelopes = []
    extension_slacks = []
    profile_upper_gaps = []

    for i in range(r):
        data = child_data(A, i, beta)
        P_child = cycle_polynomial(data["C"], beta)
        deletion_polynomials.append(P_child)
        cavity_moments.append(data["v"])
        derivative_envelopes.append(data["derivative_envelope"])
        extension_slacks.append(data["extension_product_mean"])
        profile_upper_gaps.append(data["profile_upper"] - data["kappa"])

        # Exact parent-to-child change of measure.  The coordinate flip is
        # implemented in projective representatives by re-gauging if needed;
        # h_i itself is invariant under the global re-gauging.
        local_field = (
            parent_sigma
            * parent_X[:, i]
            * (parent_X @ A[i])
        )
        affinity = float(parent_weights @ np.exp(-2 * beta * local_field))
        tilted_square = float(
            parent_weights @ (local_field**2 * np.exp(-2 * beta * local_field))
        )
        assert math.isclose(tilted_square, affinity * data["v"], rel_tol=2e-10)

        P_parent = cycle_polynomial(A, beta)
        affinity_from_cycles = P_child / (math.cosh(2 * beta) ** h * P_parent)
        assert math.isclose(affinity, affinity_from_cycles, rel_tol=2e-10)

        # Enumerate the actual complete-signing star competitors.  Their
        # degree-two Fourier sum is the signed deg_i(F)=2 cycle deck.
        twists = [star_twist(A, i, mask) for mask in range(1 << h)]
        f = np.asarray([cycle_polynomial(twist, beta) for twist in twists])
        mean_f = float(np.mean(f))
        assert math.isclose(mean_f, P_child, rel_tol=2e-10, abs_tol=1e-14)
        k2 = np.empty(1 << h)
        for mask in range(1 << h):
            k = bin(mask).count("1")
            k2[mask] = ((h - 2 * k) ** 2 - h) / 2
        a2 = float(np.mean(f * k2))
        signed_degree_two.append(a2)
        assert math.isclose(
            P_child * (data["v"] - h),
            2 * a2 / (rho * rho),
            rel_tol=4e-9,
            abs_tol=2e-12,
        )

        # Exact Q-minimality floor for every actual star competitor.  The
        # normalized oriented-projective partition has 2^r atoms.
        universal_floor = 2 * math.cosh(beta * q) / (
            2**r * math.cosh(2 * beta) ** (r * h // 2)
        )
        assert np.min(f) + 1e-14 >= universal_floor
        floor_ratio = universal_floor / P_child
        floor_ratios.append(floor_ratio)
        krawtchouk_bound = h - (
            2 * (h // 2) / (rho * rho)
        ) * (1 - floor_ratio)
        competitor_floor_bounds.append(krawtchouk_bound)
        assert data["v"] + 2e-9 >= krawtchouk_bound

        # Retaining each competitor's actual norm is a strictly stronger
        # pointwise floor and still uses no unsigned-cycle majorant.
        competitor_q = np.asarray([qnorm(twist) for twist in twists])
        q_floors = 2 * np.cosh(beta * competitor_q) / (
            2**r * math.cosh(2 * beta) ** (r * h // 2)
        )
        assert np.min(f - q_floors) >= -1e-13
        kmin = -(h // 2)
        full_a2_lower = kmin * P_child + float(np.mean(q_floors * (k2 - kmin)))
        full_q_bound = h + 2 * full_a2_lower / (rho * rho * P_child)
        full_q_floor_bounds.append(full_q_bound)
        assert data["v"] + 2e-9 >= full_q_bound

    weights = np.asarray(deletion_polynomials)
    moments = np.asarray(cavity_moments)
    a2s = np.asarray(signed_degree_two)
    weighted_v = float(weights @ moments / weights.sum())
    signed_ratio = float(a2s.sum() / weights.sum())
    assert math.isclose(
        weighted_v,
        h + 2 * signed_ratio / (rho * rho),
        rel_tol=4e-9,
        abs_tol=2e-12,
    )
    signed_deck_sufficient_rhs = (rho * rho / 2) * (threshold - h)

    print(
        name,
        f"beta={beta}",
        {
            "r": r,
            "q": q,
            "moment_threshold": threshold,
            "v_range": (min(cavity_moments), max(cavity_moments)),
            "P_deleted_weighted_v": weighted_v,
            "N2_over_N0": signed_ratio,
            "deck_sufficient_rhs": signed_deck_sufficient_rhs,
            "deck_sufficient": signed_ratio + 1e-10 >= signed_deck_sufficient_rhs,
            "competitor_floor_bound_range": (
                min(competitor_floor_bounds),
                max(competitor_floor_bounds),
            ),
            "competitor_floor_ratio_range": (min(floor_ratios), max(floor_ratios)),
            "full_competitor_Q_floor_bound_range": (
                min(full_q_floor_bounds),
                max(full_q_floor_bounds),
            ),
            "4v_over_derivative_envelope_range": (
                min(4 * v / d for v, d in zip(cavity_moments, derivative_envelopes)),
                max(4 * v / d for v, d in zip(cavity_moments, derivative_envelopes)),
            ),
            "extension_product_mean_range": (min(extension_slacks), max(extension_slacks)),
            "child_profile_kappa_upper_gap_range": (
                min(profile_upper_gaps),
                max(profile_upper_gaps),
            ),
        },
    )


def main() -> None:
    for beta in (0.25, 0.5, 1.0):
        for A, name in ((A5, "A5"), (A6, "A6"), (A8, "A8"), (A9, "A9")):
            audit_matrix(A, name, beta)
    print("PASS: signed deletion deck, child tilt, deficit square, and competitor floor")


if __name__ == "__main__":
    main()
