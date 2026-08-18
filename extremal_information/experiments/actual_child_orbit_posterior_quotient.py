#!/usr/bin/env python3
"""Exact finite audit of the actual-child latent orbit quotient.

The two order-eight matrices are the certified pressure-minimizer classes at
raw temperature ``t >= 3``.  We attach the unique order-two child and use the
negative disorder exponent ``lambda = 1``.  All structural enumerations are
integer exact.  The posterior-response comparison is evaluated in the finite
field F_p at the formal value ``z = 2``.  Nonzero denominators and distinct
residues prove distinct rational response functions; transcendence of
``exp(3)`` then proves distinction at the actual raw temperature ``t = 3``.

This is a finite theorem/falsifier, not an asymptotic claim.
"""

from __future__ import annotations

import hashlib
import itertools
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from actual_child_radial_ceiling_witness import A0, A1


PRIME = 1_000_003
RAW_T = 3
DISORDER_LAMBDA = 1
RIGHT = ((0, 1), (1, 0))


def projective_spins(n: int) -> np.ndarray:
    return np.asarray(
        [(1,) + tail for tail in itertools.product((-1, 1), repeat=n - 1)],
        dtype=np.int16,
    )


def energies(matrix: np.ndarray, spins: np.ndarray) -> np.ndarray:
    return (
        np.einsum("bi,ij,bj->b", spins, matrix, spins, optimize=True) // 2
    ).astype(np.int16)


def cosh_log_two_mod(k: int) -> int:
    """Return cosh(k log 2) in F_PRIME."""

    k = abs(int(k))
    two_k = pow(2, k, PRIME)
    return (
        (two_k + pow(two_k, PRIME - 2, PRIME))
        * pow(2, PRIME - 2, PRIME)
    ) % PRIME


def signed_similarities(
    matrix: np.ndarray, delta: int = 1
) -> list[tuple[tuple[int, ...], tuple[int, ...]]]:
    """Signed permutations U with U^T A U = delta*A, modulo global sign."""

    n = len(matrix)
    answer: list[tuple[tuple[int, ...], tuple[int, ...]]] = []
    for permutation in itertools.permutations(range(n)):
        signs = [1] * n
        for j in range(1, n):
            signs[j] = int(
                delta
                * matrix[0, j]
                * matrix[permutation[0], permutation[j]]
            )
        if all(
            signs[i]
            * signs[j]
            * matrix[permutation[i], permutation[j]]
            == delta * matrix[i, j]
            for i in range(n)
            for j in range(i + 1, n)
        ):
            answer.append((tuple(permutation), tuple(signs)))
    return answer


def signed_automorphisms(
    matrix: np.ndarray,
) -> list[tuple[tuple[int, ...], tuple[int, ...]]]:
    return signed_similarities(matrix, 1)


def spin_orbits(
    spins: np.ndarray,
    automorphisms: list[tuple[tuple[int, ...], tuple[int, ...]]],
) -> tuple[dict[tuple[int, ...], int], list[int]]:
    def act(
        spin: tuple[int, ...],
        element: tuple[tuple[int, ...], tuple[int, ...]],
    ) -> tuple[int, ...]:
        permutation, signs = element
        image = tuple(
            signs[i] * spin[permutation[i]] for i in range(len(spin))
        )
        return image if image[0] == 1 else tuple(-value for value in image)

    unseen = {tuple(map(int, spin)) for spin in spins}
    orbit_id: dict[tuple[int, ...], int] = {}
    sizes: list[int] = []
    while unseen:
        representative = next(iter(unseen))
        orbit = {act(representative, element) for element in automorphisms}
        unseen.difference_update(orbit)
        index = len(sizes)
        for spin in orbit:
            orbit_id[spin] = index
        sizes.append(len(orbit))
    return orbit_id, sizes


def simultaneous_similarity_orbits(
    left_spins: np.ndarray,
    right_spins: np.ndarray,
    left: np.ndarray,
    right: np.ndarray,
) -> tuple[dict[tuple[tuple[int, ...], tuple[int, ...]], int], list[int]]:
    """Orbits under simultaneous automorphisms or simultaneous anti-isomorphisms."""

    def act_spin(
        spin: tuple[int, ...],
        element: tuple[tuple[int, ...], tuple[int, ...]],
    ) -> tuple[int, ...]:
        permutation, signs = element
        image = tuple(
            signs[i] * spin[permutation[i]] for i in range(len(spin))
        )
        return image if image[0] == 1 else tuple(-value for value in image)

    left_states = [tuple(map(int, spin)) for spin in left_spins]
    right_states = [tuple(map(int, spin)) for spin in right_spins]
    state_set = set(itertools.product(left_states, right_states))
    transformations = []
    for delta in (1, -1):
        left_maps = signed_similarities(left, delta)
        right_maps = signed_similarities(right, delta)
        transformations.extend(itertools.product(left_maps, right_maps))
    unseen = set(state_set)
    orbit_id: dict[tuple[tuple[int, ...], tuple[int, ...]], int] = {}
    sizes: list[int] = []
    while unseen:
        x, y = next(iter(unseen))
        orbit = {
            (act_spin(x, left_map), act_spin(y, right_map))
            for left_map, right_map in transformations
        }
        if not orbit <= state_set:
            raise AssertionError("similarity action left the projective state set")
        unseen.difference_update(orbit)
        index = len(sizes)
        for state in orbit:
            orbit_id[state] = index
        sizes.append(len(orbit))
    return orbit_id, sizes


def rooted_energy_overlap_profiles(
    spins: np.ndarray, energy: np.ndarray
) -> list[tuple[tuple[tuple[int, int], int], ...]]:
    profiles = []
    for root in spins:
        profiles.append(
            tuple(
                sorted(
                    Counter(
                        (int(value), abs(int(root @ other)))
                        for value, other in zip(energy, spins)
                    ).items()
                )
            )
        )
    return profiles


def bridge_block(start: int, stop: int, dimension: int) -> np.ndarray:
    indices = np.arange(start, stop, dtype=np.uint64)
    shifts = np.arange(dimension, dtype=np.uint64)
    bits = ((indices[:, None] >> shifts) & 1).astype(np.int16)
    return (1 - 2 * bits).astype(np.int16)


def response_residues(
    words: np.ndarray, combined_energy: np.ndarray
) -> tuple[np.ndarray, int]:
    """Return the Q-dependent part of bar(mu)(Q)/mu(Q) modulo PRIME.

    At lambda=1 and z=2, if

        O(B) = sum_Q cosh(E(Q)log2) cosh(<B,Q>log2),

    then equality of posterior/prior ratios is equality of

        N(Q) = sum_B cosh(<B,Q>log2) / O(B)^2.

    The omitted normalizer is common to every Q.
    """

    dimension = words.shape[1]
    prior_weight = np.asarray(
        [cosh_log_two_mod(value) for value in combined_energy],
        dtype=np.int64,
    )
    answer = np.zeros(len(words), dtype=np.int64)
    inverse_normalizer = 0
    for start in range(0, 1 << dimension, 512):
        bridge = bridge_block(start, min(start + 512, 1 << dimension), dimension)
        correlation = bridge @ words.T
        kernel = np.empty_like(correlation, dtype=np.int64)
        for value in range(-dimension, dimension + 1):
            kernel[correlation == value] = cosh_log_two_mod(value)
        output = np.remainder(kernel @ prior_weight, PRIME)
        if np.any(output == 0):
            raise AssertionError("chosen finite-field prime kills an output denominator")
        inverse = np.asarray(
            [pow(int(value), PRIME - 2, PRIME) for value in output],
            dtype=np.int64,
        )
        inverse_normalizer = (
            inverse_normalizer + int(np.sum(inverse) % PRIME)
        ) % PRIME
        inverse_square = np.remainder(inverse * inverse, PRIME)
        # 512 * PRIME^2 is safely below the signed int64 ceiling.
        answer = np.remainder(
            answer + np.remainder(kernel.T @ inverse_square, PRIME), PRIME
        )
    if inverse_normalizer == 0:
        raise AssertionError("chosen finite-field prime kills the escort normalizer")
    return answer, inverse_normalizer


def digest(values: list[int]) -> str:
    return hashlib.sha256(
        ",".join(map(str, values)).encode("ascii")
    ).hexdigest()


def audit_class(class_id: int, raw_matrix: tuple[tuple[int, ...], ...]) -> dict:
    left = np.asarray(raw_matrix, dtype=np.int16)
    right = np.asarray(RIGHT, dtype=np.int16)
    left_spins = projective_spins(len(left))
    right_spins = projective_spins(len(right))
    left_energy = energies(left, left_spins)
    right_energy = energies(right, right_spins)

    words = np.asarray(
        [
            (x[:, None] * y[None, :]).reshape(-1)
            for x in left_spins
            for y in right_spins
        ],
        dtype=np.int16,
    )
    combined_energy = np.asarray(
        [int(ex + ey) for ex in left_energy for ey in right_energy],
        dtype=np.int16,
    )
    residues, inverse_normalizer = response_residues(words, combined_energy)

    automorphisms = signed_automorphisms(left)
    orbit_id, orbit_sizes = spin_orbits(left_spins, automorphisms)
    simultaneous_orbit_id, simultaneous_orbit_sizes = simultaneous_similarity_orbits(
        left_spins, right_spins, left, right
    )
    profiles = rooted_energy_overlap_profiles(left_spins, left_energy)

    profile_orbits: defaultdict[object, set[int]] = defaultdict(set)
    for spin, profile in zip(left_spins, profiles):
        profile_orbits[profile].add(orbit_id[tuple(map(int, spin))])
    if len(profile_orbits) != len(orbit_sizes) or any(
        len(ids) != 1 for ids in profile_orbits.values()
    ):
        raise AssertionError("rooted profiles do not equal automorphism orbits")

    joint_orbit_values: defaultdict[tuple[int, int], list[int]] = defaultdict(list)
    simultaneous_orbit_values: defaultdict[int, list[int]] = defaultdict(list)
    shell_values: defaultdict[int, set[int]] = defaultdict(set)
    index = 0
    for x in left_spins:
        left_id = orbit_id[tuple(map(int, x))]
        for ey in right_energy:
            value = int(residues[index])
            joint_orbit_values[(left_id, int(ey))].append(value)
            pair_state = (
                tuple(map(int, x)),
                tuple(map(int, right_spins[index % len(right_spins)])),
            )
            simultaneous_orbit_values[simultaneous_orbit_id[pair_state]].append(value)
            shell_values[int(combined_energy[index])].add(value)
            index += 1
    if any(len(set(values)) != 1 for values in joint_orbit_values.values()):
        raise AssertionError("posterior response is not constant on a joint orbit")
    if any(len(set(values)) != 1 for values in simultaneous_orbit_values.values()):
        raise AssertionError(
            "posterior response is not constant on a simultaneous-similarity orbit"
        )

    orbit_residues = sorted(
        {values[0] for values in joint_orbit_values.values()}
    )
    # The response uses as many exact values as there are left orbit/profile
    # cells.  Distinct residues prove distinct rational functions.
    if len(orbit_residues) != len(orbit_sizes):
        raise AssertionError("unexpected collision in the exact response quotient")
    if len(orbit_residues) != len(simultaneous_orbit_sizes):
        raise AssertionError("response did not saturate simultaneous-similarity quotient")

    return {
        "class_id": class_id,
        "left_signed_automorphism_group_order_projective": len(automorphisms),
        "left_spin_orbit_count": len(orbit_sizes),
        "left_spin_orbit_sizes": sorted(orbit_sizes),
        "rooted_energy_absolute_overlap_profile_count": len(profile_orbits),
        "joint_factor_orbit_count": len(joint_orbit_values),
        "simultaneous_similarity_orbit_count": len(simultaneous_orbit_sizes),
        "simultaneous_similarity_orbit_sizes": sorted(simultaneous_orbit_sizes),
        "distinct_posterior_response_residues": len(orbit_residues),
        "combined_energy_shell_count": len(shell_values),
        "response_values_per_energy_shell": {
            str(shell): len(values)
            for shell, values in sorted(shell_values.items())
        },
        "inverse_escort_normalizer_residue": inverse_normalizer,
        "posterior_response_residue_sha256": digest(orbit_residues),
    }


def main() -> None:
    records = [audit_class(0, A0), audit_class(1, A1)]
    result = {
        "schema": "actual-child-orbit-posterior-quotient-v1",
        "classification": "exact finite theorem and modular certificate",
        "raw_actual_temperature": RAW_T,
        "negative_disorder_lambda": DISORDER_LAMBDA,
        "left_order": 8,
        "right_order": 2,
        "finite_field_prime": PRIME,
        "evaluation_point": "formal z evaluated at 2 in F_p",
        "records": records,
        "logical_certificate": (
            "Distinct finite-field residues imply distinct Q(z)-valued "
            "posterior-response rational functions. Since exp(3) is "
            "transcendental, those functions remain distinct at the actual "
            "raw temperature t=3."
        ),
        "scope": (
            "The finite actual optimizers admit an exact rooted-profile/orbit "
            "quotient, but this gives no all-order bound on orbit count and no "
            "row-product-regret theorem."
        ),
    }
    output = Path(
        "computations/results/actual_child_orbit_posterior_quotient.json"
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
