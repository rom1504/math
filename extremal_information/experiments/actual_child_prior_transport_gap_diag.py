#!/usr/bin/env python3
"""Finite PT.1 candidate diagnostic on certified actual children.

This is deliberately a candidate-family test, not minimization over the full
child-prior convex hull.  It evaluates covariance-top radial candidates and
posterior barycenters selected by collision severity on complete bridge cubes.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy
from scipy.optimize import linprog


ROOT = Path(__file__).resolve().parents[2]
EXP = ROOT / "extremal_information" / "experiments"
sys.path.insert(0, str(EXP))

import actual_child_bridge_law_exact as exact  # noqa: E402
import actual_child_latent_collision_diag as collision  # noqa: E402
import actual_child_escort_low_degree_falsifier as low_degree  # noqa: E402


BETA = 4.0
RAY_FRACTIONS = (0.01, 0.05, 0.10, 0.25, 0.50, 0.90)
TOP_COUNT = 24


def latent_prior(left, right, total_order, orientation):
    m, n = len(left), len(right)
    t = BETA / math.sqrt(total_order)
    x = exact.projective_spins(m).astype(np.int16)
    y = exact.projective_spins(n).astype(np.int16)
    ex = exact.energies_for_matrix(left, x)
    ey = exact.energies_for_matrix(right, y)
    atoms = []
    weights = []
    for xi, exi in zip(x, ex):
        for yj, eyj in zip(y, ey):
            atoms.append((xi[:, None] * yj[None, :]).reshape(-1))
            weights.append(math.cosh(t * float(exi + orientation * eyj)))
    atoms = np.asarray(atoms, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)
    weights /= float(np.sum(weights))
    return atoms, weights


def normalized_gap(A, atoms, weights):
    d = atoms.shape[1]
    norm_square = float(np.dot(A, A))
    scores = atoms @ A / (2.0 * math.sqrt(d))
    # The prior is centrally symmetrized, hence the cosh MGF.
    # Every candidate lies in conv(+-atoms), so |score|<=sqrt(d)/2 at the
    # tested orders and direct evaluation is safely in range.
    log_mgf = math.log(float(np.dot(weights, np.cosh(scores))))
    return {
        "frobenius_norm_square": norm_square,
        "frobenius_norm_square_over_d": norm_square / d,
        "log_mgf": log_mgf,
        "normalized_transport_gap": (
            0.5 - math.sqrt(d) * log_mgf / norm_square
        ),
    }


def covariance_ray(atoms, weights):
    covariance = atoms.T @ (weights[:, None] * atoms)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    top_value = float(eigenvalues[-1])
    direction = eigenvectors[:, -1]
    direction /= np.linalg.norm(direction)

    # Atomic norm of the top direction relative to conv(+-atoms).  If
    # min ||u||_1 subject to atoms^T u=direction is tau, then the radial
    # boundary is direction/tau.
    k = len(atoms)
    objective = np.ones(2 * k, dtype=np.float64)
    equality = np.concatenate([atoms.T, -atoms.T], axis=1)
    result = linprog(
        objective,
        A_eq=equality,
        b_eq=direction,
        bounds=(0.0, None),
        method="highs",
    )
    if not result.success:
        raise RuntimeError(result.message)
    atomic_norm = float(result.fun)
    radius = 1.0 / atomic_norm
    reconstruction = equality @ result.x
    reconstruction_error = float(np.max(np.abs(reconstruction - direction)))
    tangent_gap = 0.5 - top_value / (8.0 * math.sqrt(atoms.shape[1]))
    ray = []
    for fraction in RAY_FRACTIONS:
        candidate = fraction * radius * direction
        row = normalized_gap(candidate, atoms, weights)
        row["fraction_of_atomic_radial_boundary"] = fraction
        ray.append(row)
    return {
        "top_covariance_eigenvalue": top_value,
        "tangent_normalized_gap": tangent_gap,
        "top_direction_atomic_norm": atomic_norm,
        "top_direction_radial_boundary": radius,
        "LP_reconstruction_max_error": reconstruction_error,
        "ray_candidates": ray,
    }


def posterior_barycenter(index, log_z_t, raw_t, dimension):
    rho = math.tanh(raw_t)
    answer = np.empty(dimension, dtype=np.float64)
    for edge in range(dimension):
        bit = 1.0 - 2.0 * float((index >> edge) & 1)
        flipped = index ^ (1 << edge)
        r = (
            bit
            * math.tanh(0.5 * (log_z_t[index] - log_z_t[flipped]))
            / rho
        )
        answer[edge] = (r + rho * bit) / (1.0 + rho * bit * r)
    return answer


def run_record(record, child_cache):
    total_order = record["N"]
    left_order, right_order = record["split"]
    orientation = record["orientation"]
    left, left_certificate = low_degree.child_record(
        left_order, total_order, child_cache
    )
    right, right_certificate = low_degree.child_record(
        right_order, total_order, child_cache
    )
    atoms, weights = latent_prior(left, right, total_order, orientation)
    covariance = covariance_ray(atoms, weights)

    z_t, z_2t, z_0, audit = collision.channel_arrays(
        left, right, total_order, orientation
    )
    log_z_t = np.log(z_t).astype(np.float64)
    log_k = (
        float(np.log(z_0))
        + np.log(z_2t).astype(np.float64)
        - 2.0 * log_z_t
    )
    shifted = -log_z_t
    shifted -= float(np.max(shifted))
    q = np.exp(shifted)
    q /= float(np.sum(q))
    complement_mask = len(q) - 1
    canonical = np.arange(len(q), dtype=np.int64)
    canonical = canonical[canonical < (canonical ^ complement_mask)]
    top_collision = canonical[
        np.argsort(-log_k[canonical], kind="stable")[:TOP_COUNT]
    ]
    log_contribution = np.log(q) + log_k
    top_contribution = canonical[
        np.argsort(-log_contribution[canonical], kind="stable")[:TOP_COUNT]
    ]
    chosen = sorted(set(top_collision.tolist()) | set(top_contribution.tolist()))

    barycenters = []
    raw_t = BETA / math.sqrt(total_order)
    maximum_direct_posterior_error = 0.0
    for index in chosen:
        A = posterior_barycenter(index, log_z_t, raw_t, atoms.shape[1])
        bridge = 1.0 - 2.0 * (
            (index >> np.arange(atoms.shape[1], dtype=np.int64)) & 1
        )
        latent_field = raw_t * (atoms @ bridge)
        direct_A = (
            (weights * np.sinh(latent_field)) @ atoms
            / float(np.dot(weights, np.cosh(latent_field)))
        )
        maximum_direct_posterior_error = max(
            maximum_direct_posterior_error,
            float(np.max(np.abs(A - direct_A))),
        )
        row = normalized_gap(A, atoms, weights)
        row.update(
            {
                "bridge_mask": int(index),
                "log_collision_over_N": float(log_k[index] / total_order),
                "q_mass": float(q[index] + q[index ^ complement_mask]),
                "selected_by_top_collision": bool(index in set(top_collision)),
                "selected_by_top_qK_contribution": bool(
                    index in set(top_contribution)
                ),
            }
        )
        barycenters.append(row)
    barycenters.sort(key=lambda row: row["normalized_transport_gap"])
    if maximum_direct_posterior_error > 2e-10:
        raise AssertionError("flip posterior did not match direct Gibbs posterior")
    return {
        "N": total_order,
        "split": [left_order, right_order],
        "orientation": orientation,
        "child_selection": {
            "left": left_certificate,
            "right": right_certificate,
        },
        "latent_projective_atoms_before_central_symmetrization": len(atoms),
        "channel_audit": audit,
        "covariance_top_direction": covariance,
        "posterior_candidate_count": len(barycenters),
        "maximum_flip_vs_direct_posterior_coordinate_error": (
            maximum_direct_posterior_error
        ),
        "minimum_posterior_barycenter_gap": barycenters[0],
        "posterior_barycenter_candidates": barycenters,
    }


def main():
    mp.mp.dps = 80
    records = []
    child_cache = {}
    for total_order, left_order, right_order in collision.PLANS:
        for orientation in collision.ORIENTATIONS:
            source_record = {
                "N": total_order,
                "split": [left_order, right_order],
                "orientation": orientation,
            }
            result = run_record(source_record, child_cache)
            records.append(result)
            cov = result["covariance_top_direction"]
            post = result["minimum_posterior_barycenter_gap"]
            print(
                f"N={result['N']} eps={result['orientation']:+d} "
                f"cov_tan={cov['tangent_normalized_gap']:.9g} "
                f"cov_ray_min={min(x['normalized_transport_gap'] for x in cov['ray_candidates']):.9g} "
                f"post_min={post['normalized_transport_gap']:.9g} "
                f"post_norm2/d={post['frobenius_norm_square_over_d']:.6g}",
                flush=True,
            )
    output = {
        "schema": "actual-child-prior-transport-gap-candidates-v1",
        "classification": (
            "exhaustive child/signing and complete bridge enumeration; "
            "interval child-pressure certification; numerical covariance, "
            "linear-program, posterior, and MGF candidate evaluation"
        ),
        "parameters": {
            "beta": BETA,
            "lambda_for_bridge_candidate_ranking": collision.LAMBDA,
            "plans": [list(plan) for plan in collision.PLANS],
            "orientations": list(collision.ORIENTATIONS),
            "covariance_ray_fractions": list(RAY_FRACTIONS),
            "top_bridge_pairs_per_ranking": TOP_COUNT,
            "mp_dps_for_child_selection": mp.mp.dps,
            "numpy_version": np.__version__,
            "scipy_version": scipy.__version__,
        },
        "PT1_normalized_gap": (
            "c(A)=1/2-sqrt(d) log E_mu exp(<A,Q>/(2sqrt(d)))/||A||_F^2; "
            "PT.1 with gap c requires c(A)>=c for every nonzero A in conv(supp mu)"
        ),
        "candidate_family": {
            "covariance": (
                "top eigenvector of E_mu vec(Q)vec(Q)^T, evaluated along fixed "
                "fractions of its exact LP-computed radial interval in conv(+-supp mu)"
            ),
            "posterior": (
                "full posterior barycenters from the 24 largest K0 bridge pairs "
                "and 24 largest per-word qK0 contribution bridge pairs"
            ),
        },
        "records": records,
        "scope": {
            "exhaustive": [
                "child signing and exact energy-histogram enumeration",
                "both declared orientations of every complete bridge cube",
                "covariance matrix of every latent prior atom",
            ],
            "candidate_only": [
                "six radial scales on the covariance top eigendirection",
                "24 antipodal bridge pairs with largest K0",
                "24 antipodal bridge pairs with largest per-word qK0 contribution",
            ],
            "not_claimed": [
                "minimization of the PT.1 gap over the complete convex hull",
                "a positive all-candidate or all-order transport gap",
                "an asymptotic lower bound",
                "a Level-6 recurrence implication",
            ],
        },
    }
    output_path = (
        ROOT
        / "computations/results/actual_child_prior_transport_gap_diag.json"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
