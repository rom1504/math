#!/usr/bin/env python3
"""Exact finite scaling audit for the actual-child orbit quotient.

This script uses only the stored exhaustive cap-minimizer classifications at
orders 3 through 8.  At raw temperature ``t=3`` those classifications also
certify the thermal pressure minimizers (see ``pressure_certificate`` below).
Each eligible left child is paired with the unique order-two child and the
negative-disorder exponent is ``lambda=1``.

All orbit/profile enumerations are integer exact.  Posterior-response
functions are evaluated at ``z=2`` in a prime field.  Pairwise-distinct
residues prove distinct rational functions over Q(z), and hence distinction
at the actual value ``z=exp(3)`` by transcendence.  The sole modular collision,
at order three, is proved to be an identity by an exact denominator-type
cancellation.

No order-nine pressure classification is available in the repository.  The
stored order-nine matrix is therefore reported only as an explicitly
ineligible ground-state diagnostic; it is not used as an actual child.
"""

from __future__ import annotations

import itertools
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

import actual_child_orbit_posterior_quotient as orbit


ACTUAL_CLASS_IDS = {
    3: (0,),
    4: (0,),
    5: (0,),
    6: (0,),
    7: (2,),
    8: (0, 1),
}


def load_classes(order: int) -> list[dict]:
    path = Path(f"computations/results/m{order}_minimizer_orbits.json")
    return json.loads(path.read_text())["classes"]


def class_histogram(row: dict) -> dict[int, int]:
    return {
        int(energy): int(count)
        for energy, count in row["representative_profile"][
            "energy_histogram"
        ].items()
    }


def absolute_histogram(row: dict) -> dict[int, int]:
    answer: Counter[int] = Counter()
    for energy, count in class_histogram(row).items():
        answer[abs(energy)] += count
    return dict(sorted(answer.items()))


def pressure_certificate() -> dict:
    """Record the exact finite argument selecting the actual classes.

    Every non-cap-minimizing signing has cap at least ``M_n+2``.  For
    ``3<=n<=8`` and ``t>=3``,

        cosh((M_n+2)t) / (2^(n-1) cosh(M_n t))
            >= exp(2t) / 2^n > 1.

    Hence it remains only to compare the exhaustively classified cap classes.
    There is one class through order six; the two order-eight classes have the
    same absolute-energy histogram.  At order seven the exact differences
    between the projective partition sums of classes 0/1 and class 2 are the
    displayed positive cosh combinations.
    """

    classes7 = load_classes(7)
    histograms7 = [absolute_histogram(row) for row in classes7]
    assert histograms7 == [
        {1: 15, 3: 21, 5: 15, 7: 9, 9: 4},
        {1: 21, 3: 21, 5: 7, 7: 8, 9: 7},
        {1: 21, 3: 13, 5: 15, 7: 12, 9: 3},
    ]
    classes8 = load_classes(8)
    assert absolute_histogram(classes8[0]) == absolute_histogram(classes8[1])
    return {
        "raw_temperature_range": "t>=3",
        "nonminimal_cap_exclusion": (
            "cosh((M+2)t)/(2^(n-1)cosh(Mt)) >= exp(2t)/2^n > 1 "
            "for n<=8 and t>=3"
        ),
        "order_7_projective_absolute_histograms": histograms7,
        "order_7_class_0_minus_class_2_partition_sum": (
            "cosh(9t)-3cosh(7t)+8cosh(3t)-6cosh(t)>0 for t>=3"
        ),
        "order_7_class_1_minus_class_2_partition_sum": (
            "4cosh(9t)-4cosh(7t)-8cosh(5t)+8cosh(3t)>0 for t>=3"
        ),
        "positivity_reason": (
            "cosh(9t)/cosh(7t)>=exp(2t)/2>9 and cosh(7t)>=cosh(5t)>=cosh(t)"
        ),
        "actual_class_ids": {
            str(order): list(ids) for order, ids in ACTUAL_CLASS_IDS.items()
        },
    }


def exact_order_three_collision(
    words: np.ndarray,
    combined_energy: np.ndarray,
    response: np.ndarray,
    simultaneous_id: dict,
    left_spins: np.ndarray,
    right_spins: np.ndarray,
) -> dict:
    """Prove the unique order-three response collision for every z.

    For a bridge B, its output denominator is determined by the histogram of
    pairs ``(|E(Q)|, |<B,Q>|)``.  Within each denominator type, the coefficient
    of every kernel ``C_k(z)`` in the difference of the two colliding response
    numerators vanishes.  This is an exact rational-function identity.
    """

    orbit_indices: defaultdict[int, list[int]] = defaultdict(list)
    index = 0
    for x in left_spins:
        for y in right_spins:
            state = (tuple(map(int, x)), tuple(map(int, y)))
            orbit_indices[simultaneous_id[state]].append(index)
            index += 1
    representative = {key: values[0] for key, values in orbit_indices.items()}
    residue_to_orbits: defaultdict[int, list[int]] = defaultdict(list)
    for key, index in representative.items():
        residue_to_orbits[int(response[index])].append(key)
    collision = [keys for keys in residue_to_orbits.values() if len(keys) > 1]
    assert len(collision) == 1 and len(collision[0]) == 2
    first, second = [representative[key] for key in collision[0]]

    net_by_denominator: defaultdict[tuple, Counter[int]] = defaultdict(Counter)
    for raw_bridge in itertools.product((-1, 1), repeat=words.shape[1]):
        bridge = np.asarray(raw_bridge, dtype=np.int16)
        correlations = [abs(int(bridge @ word)) for word in words]
        denominator_type = tuple(
            sorted(
                Counter(
                    zip(map(abs, map(int, combined_energy)), correlations)
                ).items()
            )
        )
        net_by_denominator[denominator_type][correlations[first]] += 1
        net_by_denominator[denominator_type][correlations[second]] -= 1
    assert all(not +counter for counter in net_by_denominator.values())
    return {
        "simultaneous_orbit_representative_indices": [first, second],
        "denominator_types": len(net_by_denominator),
        "certificate": (
            "For every output-denominator type, every C_k(z) coefficient "
            "in the numerator difference is exactly zero. Therefore the two "
            "response rational functions agree for every z and every t."
        ),
    }


def audit_actual_class(order: int, class_row: dict) -> dict:
    left = np.asarray(class_row["representative_matrix"], dtype=np.int16)
    right = np.asarray(orbit.RIGHT, dtype=np.int16)
    left_spins = orbit.projective_spins(order)
    right_spins = orbit.projective_spins(2)
    left_energy = orbit.energies(left, left_spins)
    right_energy = orbit.energies(right, right_spins)
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
    response, inverse_normalizer = orbit.response_residues(
        words, combined_energy
    )

    automorphisms = orbit.signed_automorphisms(left)
    spin_id, spin_sizes = orbit.spin_orbits(left_spins, automorphisms)
    simultaneous_id, simultaneous_sizes = orbit.simultaneous_similarity_orbits(
        left_spins, right_spins, left, right
    )
    profiles = orbit.rooted_energy_overlap_profiles(left_spins, left_energy)

    profile_ids = {profile: i for i, profile in enumerate(dict.fromkeys(profiles))}
    profile_to_spin_orbits: defaultdict[int, set[int]] = defaultdict(set)
    spin_orbit_to_profiles: defaultdict[int, set[int]] = defaultdict(set)
    simultaneous_response: defaultdict[int, set[int]] = defaultdict(set)
    index = 0
    for x_index, x in enumerate(left_spins):
        spin_cell = spin_id[tuple(map(int, x))]
        profile_cell = profile_ids[profiles[x_index]]
        profile_to_spin_orbits[profile_cell].add(spin_cell)
        spin_orbit_to_profiles[spin_cell].add(profile_cell)
        for y in right_spins:
            state = (tuple(map(int, x)), tuple(map(int, y)))
            simultaneous_response[simultaneous_id[state]].add(
                int(response[index])
            )
            index += 1
    assert all(len(values) == 1 for values in simultaneous_response.values())

    distinct_response = len(set(map(int, response)))
    simultaneous_count = len(simultaneous_sizes)
    record = {
        "order": order,
        "class_id": int(class_row["class"]),
        "projective_spin_count": len(left_spins),
        "signed_automorphism_group_order_projective": len(automorphisms),
        "spin_automorphism_orbit_count": len(spin_sizes),
        "rooted_energy_absolute_overlap_profile_count": len(profile_ids),
        "rooted_profiles_equal_spin_orbits": (
            len(profile_ids) == len(spin_sizes)
            and all(len(values) == 1 for values in profile_to_spin_orbits.values())
            and all(len(values) == 1 for values in spin_orbit_to_profiles.values())
        ),
        "simultaneous_signed_similarity_orbit_count": simultaneous_count,
        "distinct_negative_posterior_response_functions": distinct_response,
        "response_saturates_simultaneous_quotient": (
            distinct_response == simultaneous_count
        ),
        "inverse_escort_normalizer_residue": inverse_normalizer,
        "response_residue_sha256": orbit.digest(sorted(set(map(int, response)))),
    }
    if order == 3:
        record["exact_nonsaturation_identity"] = exact_order_three_collision(
            words,
            combined_energy,
            response,
            simultaneous_id,
            left_spins,
            right_spins,
        )
        assert distinct_response == 3 and simultaneous_count == 4
    else:
        assert distinct_response == simultaneous_count
    return record


def order_two_record() -> dict:
    matrix = np.asarray(orbit.RIGHT, dtype=np.int16)
    row = {
        "class": 0,
        "representative_matrix": matrix.tolist(),
    }
    return audit_actual_class(2, row)


def ineligible_order_nine_diagnostic() -> dict:
    payload = json.loads(Path("computations/results/exact_m9.json").read_text())
    matrix = np.asarray(payload["matrix"], dtype=np.int16)
    spins = orbit.projective_spins(9)
    energy = orbit.energies(matrix, spins)
    automorphisms = orbit.signed_automorphisms(matrix)
    spin_id, spin_sizes = orbit.spin_orbits(spins, automorphisms)
    profiles = orbit.rooted_energy_overlap_profiles(spins, energy)
    profile_to_orbits: defaultdict[object, set[int]] = defaultdict(set)
    for spin, profile in zip(spins, profiles):
        profile_to_orbits[profile].add(spin_id[tuple(map(int, spin))])
    return {
        "eligible_as_actual_pressure_child": False,
        "reason": (
            "The repository stores one certified cap minimizer, but no "
            "exhaustive thermal-pressure-minimizer classification at order 9."
        ),
        "projective_spin_count": len(spins),
        "signed_automorphism_group_order_projective": len(automorphisms),
        "spin_automorphism_orbit_count": len(spin_sizes),
        "rooted_energy_absolute_overlap_profile_count": len(set(profiles)),
        "rooted_profiles_equal_spin_orbits": (
            len(set(profiles)) == len(spin_sizes)
            and all(len(values) == 1 for values in profile_to_orbits.values())
        ),
    }


def main() -> None:
    records = [order_two_record()]
    for order, class_ids in ACTUAL_CLASS_IDS.items():
        rows = load_classes(order)
        for class_id in class_ids:
            records.append(audit_actual_class(order, rows[class_id]))
    result = {
        "schema": "actual-child-orbit-scaling-v1",
        "classification": (
            "exact finite actual-pressure theorem and modular certificate"
        ),
        "raw_actual_temperature": 3,
        "negative_disorder_lambda": 1,
        "right_child_order": 2,
        "finite_field_prime": orbit.PRIME,
        "finite_field_evaluation": "formal z evaluated at 2 in F_p",
        "pressure_certificate": pressure_certificate(),
        "records": records,
        "order_9_ground_state_only_diagnostic": ineligible_order_nine_diagnostic(),
        "logical_certificate": (
            "Distinct residues at z=2 imply distinct Q(z)-valued response "
            "functions; a nonzero rational function cannot vanish at the "
            "transcendental actual point z=exp(3). Order 3 is handled by a "
            "separate exact all-z cancellation certificate."
        ),
        "scope": (
            "The actual finite response usually saturates its full available "
            "simultaneous symmetry quotient. This is not an all-order orbit "
            "bound and gives no row-product-regret or recurrence theorem."
        ),
    }
    output = Path("computations/results/actual_child_orbit_scaling.json")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
