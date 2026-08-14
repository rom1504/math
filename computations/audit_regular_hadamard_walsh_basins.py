#!/usr/bin/env python3
"""Walsh-coordinate and basin audit for the regular-Hadamard obstruction.

The audit deliberately separates two experiments.

``local`` starts a Boolean spin uniformly at random and applies the prescribed
largest-gain/least-index one-flip ascent to ``-r^T (K-I) r``.  Its bad-basin
event means that the terminal spin belongs to the affine-orthogonal orbit of
the explicit codimension-two plateaued spin.

``full`` samples the row-sign input of the two-block signing from the
obstruction note.  It constructs both anchored augmented problems in switched
coordinates, applies the same deterministic ascent, and records the better
defect and unmatched-core number.  This is the probability relevant to the
expected row-sign certificate; it is not inferred from the local experiment.

All randomness uses numpy's PCG64 generator with the seed stored in the JSON.
Temporary outputs should be placed below /home/math/quadra/tmp.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter, deque
from functools import lru_cache
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np


def parity_table(n: int) -> np.ndarray:
    return np.fromiter(
        (bin(i).count("1") & 1 for i in range(n)), dtype=np.int8, count=n
    )


@lru_cache(maxsize=None)
def quadratic_data(k: int) -> Dict[str, np.ndarray]:
    """Return f, K, C, symplectic Walsh H, and the explicit bad spin."""
    s = 1 << k
    m = s * s
    values = np.arange(m, dtype=np.int64)
    a = values & (s - 1)
    b = values >> k
    parity = parity_table(s)
    f = (1 - 2 * parity[np.bitwise_and(a, b)]).astype(np.int32)
    k_matrix = f[np.bitwise_xor(values[:, None], values[None, :])]
    c_matrix = k_matrix - np.eye(m, dtype=np.int32)
    h_matrix = f[:, None] * k_matrix * f[None, :]
    in_subspace = (
        (((b >> 1) & 1) == ((a >> 1) & 1))
        & ((b & 1) == ((a & 1) ^ ((a >> 1) & 1)))
    )
    bad_spin = np.where(in_subspace, -1, 1).astype(np.int32)
    ones = np.ones(m, dtype=np.int32)
    assert np.array_equal(k_matrix @ k_matrix, m * np.eye(m, dtype=np.int32))
    assert np.array_equal(k_matrix @ ones, s * ones)
    assert np.array_equal(h_matrix @ f, s * f)
    assert np.array_equal(f[:, None] * h_matrix * f[None, :], k_matrix)
    return {
        "f": f,
        "K": k_matrix,
        "C": c_matrix,
        "H": h_matrix,
        "bad_spin": bad_spin,
        "in_subspace": in_subspace,
    }


def spins_from_integer(value: int, order: int) -> np.ndarray:
    return np.fromiter(
        (-1 if ((value >> i) & 1) else 1 for i in range(order)),
        dtype=np.int32,
        count=order,
    )


def negative_mask(spin: np.ndarray) -> int:
    answer = 0
    for i in np.flatnonzero(spin < 0):
        answer |= 1 << int(i)
    return answer


def popcount(value: int) -> int:
    # Python 3.8 compatibility for the repository environment.
    return bin(value).count("1")


def q_value(u: int, k: int) -> int:
    mask = (1 << k) - 1
    return popcount((u & mask) & (u >> k)) & 1


def beta_value(u: int, v: int, k: int) -> int:
    mask = (1 << k) - 1
    a, b = u & mask, u >> k
    c, d = v & mask, v >> k
    return (popcount(a & d) + popcount(c & b)) & 1


def permutation_translation(order: int, t: int) -> Tuple[int, ...]:
    return tuple(u ^ t for u in range(order))


def permutation_orthogonal_transvection(
    order: int, k: int, v: int
) -> Tuple[int, ...]:
    # q(v)=1 makes u -> u+beta(u,v)v an isometry of q in characteristic 2.
    assert q_value(v, k) == 1
    return tuple(u ^ (v if beta_value(u, v, k) else 0) for u in range(order))


def permutation_swap_halves(order: int, k: int) -> Tuple[int, ...]:
    """The q-isometry (a,b)->(b,a), outside the transvection subgroup at k=2."""
    mask = (1 << k) - 1
    return tuple(((u & mask) << k) | (u >> k) for u in range(order))


def permutation_orthogonal_shear(order: int, k: int) -> Tuple[int, ...]:
    """The q-isometry adding (b_1,b_0,0,...) to a."""
    if k < 2:
        raise ValueError("the shear needs two coordinates")
    mask = (1 << k) - 1
    answer = []
    for u in range(order):
        a, b = u & mask, u >> k
        shear = ((b >> 1) & 1) | ((b & 1) << 1)
        answer.append((a ^ shear) | (b << k))
    return tuple(answer)


def permute_mask(mask: int, permutation: Sequence[int]) -> int:
    answer = 0
    for old, new in enumerate(permutation):
        if (mask >> old) & 1:
            answer |= 1 << int(new)
    return answer


@lru_cache(maxsize=None)
def affine_orthogonal_bad_orbit(k: int) -> Set[int]:
    """Orbit under translations and orthogonal transvections, including sign."""
    data = quadratic_data(k)
    order = len(data["f"])
    start = negative_mask(data["bad_spin"])
    generators: List[Tuple[int, ...]] = [
        permutation_translation(order, 1 << i) for i in range(2 * k)
    ]
    generators.append(permutation_swap_halves(order, k))
    generators.append(permutation_orthogonal_shear(order, k))
    generators.extend(
        permutation_orthogonal_transvection(order, k, v)
        for v in range(1, order)
        if q_value(v, k) == 1
    )
    orbit = {start}
    queue = deque([start])
    while queue:
        state = queue.popleft()
        for permutation in generators:
            image = permute_mask(state, permutation)
            if image not in orbit:
                orbit.add(image)
                queue.append(image)
    full_mask = (1 << order) - 1
    orbit.update(full_mask ^ state for state in tuple(orbit))

    # Directly check that every enumerated spin has the claimed local fields.
    c_matrix = data["C"]
    s = 1 << k
    for state in orbit:
        spin = spins_from_integer(state, order)
        products = spin * (c_matrix @ spin)
        if not (
            int(-(spin @ c_matrix @ spin)) == order * (s // 2 + 1)
            and np.all(products <= -1)
        ):
            raise AssertionError("affine-orthogonal image lost bad stability")
    return orbit


def k2_symmetry_invariant_tie_audit() -> Dict[str, object]:
    """Exact projective quotient for the symmetry-invariant random tie rule.

    The repository's least-index rule is not permutation equivariant.  This
    auxiliary calculation replaces a tied largest gain by a uniform choice;
    it is therefore an honest quotient under the affine-orthogonal subgroup.
    It is reported only as a tie-sensitivity control, never as the row-sign
    basin probability.
    """
    k = 2
    data = quadratic_data(k)
    matrix = data["C"]
    order = len(matrix)
    total_states = 1 << order
    full_mask = total_states - 1
    bad_orbit = affine_orthogonal_bad_orbit(k)
    generators: List[Tuple[int, ...]] = [
        permutation_translation(order, 1 << i) for i in range(2 * k)
    ]
    generators.append(permutation_swap_halves(order, k))
    generators.append(permutation_orthogonal_shear(order, k))
    generators.extend(
        permutation_orthogonal_transvection(order, k, v)
        for v in range(1, order)
        if q_value(v, k) == 1
    )

    @lru_cache(maxsize=None)
    def outcome(state: int) -> Tuple[float, bool, bool]:
        spin = spins_from_integer(state, order)
        gains = 4 * spin * (matrix @ spin)
        best = int(gains.max())
        if best <= 0:
            hit = state in bad_orbit
            return float(hit), hit, hit
        moves = np.flatnonzero(gains == best)
        children = [outcome(state ^ (1 << int(i))) for i in moves]
        return (
            sum(child[0] for child in children) / len(children),
            any(child[1] for child in children),
            all(child[2] for child in children),
        )

    unseen = set(range(total_states))
    orbit_representatives: List[Tuple[int, int]] = []
    while unseen:
        representative = min(unseen)
        orbit = {representative, full_mask ^ representative}
        queue = deque(orbit)
        while queue:
            state = queue.popleft()
            for permutation in generators:
                image = permute_mask(state, permutation)
                for candidate in (image, full_mask ^ image):
                    if candidate not in orbit:
                        orbit.add(candidate)
                        queue.append(candidate)
        unseen.difference_update(orbit)
        orbit_representatives.append((representative, len(orbit)))

    weighted_probability = 0.0
    weighted_possible = 0
    weighted_forced = 0
    maximum_invariance_error = 0.0
    for representative, size in orbit_representatives:
        probability, possible, forced = outcome(representative)
        weighted_probability += size * probability
        weighted_possible += size * int(possible)
        weighted_forced += size * int(forced)
        # Generator checks are sufficient because the random-tie recurrence is
        # invariant and the subgroup is generated by these maps.
        for permutation in generators:
            image_probability = outcome(permute_mask(representative, permutation))[0]
            maximum_invariance_error = max(
                maximum_invariance_error, abs(probability - image_probability)
            )
    if maximum_invariance_error > 1e-12:
        raise AssertionError("random-tie basin failed subgroup invariance")
    return {
        "k": 2,
        "m": order,
        "subgroup": (
            "translations semidirect q-isometries generated by nonsingular "
            "transvections, half-swap, and an orthogonal shear; modulo global sign"
        ),
        "projective_input_orbit_count": len(orbit_representatives),
        "uniform_random_largest_gain_tie_probability": weighted_probability
        / total_states,
        "some_largest_gain_path_probability": weighted_possible / total_states,
        "every_largest_gain_path_probability": weighted_forced / total_states,
        "maximum_generator_invariance_error": maximum_invariance_error,
        "warning": (
            "auxiliary symmetry-invariant tie rule; the prescribed least-index "
            "probability is the exact local k=2 entry"
        ),
    }


def best_improvement(
    matrix: np.ndarray, spin: np.ndarray, target_sign: int
) -> Tuple[np.ndarray, int, int]:
    """Largest positive gain, least coordinate on ties, exact integer update."""
    spin = spin.astype(np.int32, copy=True)
    field = matrix @ spin
    flips = 0
    edge_mass = int(np.abs(np.triu(matrix, 1)).sum())
    while flips <= edge_mass:
        gains = -4 * target_sign * spin * field
        coordinate = int(np.argmax(gains))
        if int(gains[coordinate]) <= 0:
            energy = int(target_sign * (spin @ matrix @ spin))
            assert np.all(target_sign * spin * field >= 0)
            return spin, energy, flips
        old = int(spin[coordinate])
        spin[coordinate] = -old
        field -= 2 * old * matrix[:, coordinate]
        flips += 1
    raise AssertionError("strict ascent exceeded integral edge-mass bound")


def walsh_metrics(
    spin_in_translation_coordinates: np.ndarray,
    data: Dict[str, np.ndarray],
    bad_orbit: Optional[Set[int]],
) -> Dict[str, float]:
    """Metrics of z=f*x in the symplectic Walsh coordinates H."""
    f = data["f"]
    h_matrix = data["H"]
    order = len(f)
    s = int(round(math.sqrt(order)))
    z = f * spin_in_translation_coordinates
    spectrum = h_matrix @ z
    abs_spectrum = np.abs(spectrum)
    bent_affine_spectrum = h_matrix @ (f * z)
    result = {
        "max_abs_over_sqrt_m": float(abs_spectrum.max() / s),
        "support_fraction": float(np.count_nonzero(spectrum) / order),
        "fourth_ratio": float(np.sum(spectrum.astype(np.float64) ** 4) / order**3),
        "rm1_distance_fraction": float((order - abs_spectrum.max()) / (2 * order)),
        "distance_to_explicit_bent_affine_orbit_fraction": float(
            (order - np.max(np.abs(bent_affine_spectrum))) / (2 * order)
        ),
    }
    if bad_orbit is not None:
        state = negative_mask(spin_in_translation_coordinates)
        result["distance_to_bad_plateaued_orbit_fraction"] = min(
            popcount(state ^ candidate) for candidate in bad_orbit
        ) / order
    return result


class MetricAccumulator:
    def __init__(self) -> None:
        self.count = 0
        self.sums: Counter[str] = Counter()
        self.sumsq: Counter[str] = Counter()
        self.minimum: Dict[str, float] = {}
        self.maximum: Dict[str, float] = {}

    def add(self, values: Dict[str, float]) -> None:
        self.count += 1
        for key, raw in values.items():
            value = float(raw)
            self.sums[key] += value
            self.sumsq[key] += value * value
            self.minimum[key] = min(self.minimum.get(key, value), value)
            self.maximum[key] = max(self.maximum.get(key, value), value)

    def report(self) -> Dict[str, object]:
        if not self.count:
            return {"count": 0}
        answer: Dict[str, object] = {"count": self.count}
        for key in sorted(self.sums):
            mean = self.sums[key] / self.count
            variance = max(0.0, self.sumsq[key] / self.count - mean * mean)
            answer[key] = {
                "mean": mean,
                "standard_deviation": math.sqrt(variance),
                "minimum": self.minimum[key],
                "maximum": self.maximum[key],
            }
        return answer


def run_local(k: int, samples: int, seed: int, exact: bool) -> Dict[str, object]:
    data = quadratic_data(k)
    order = len(data["f"])
    c_matrix = data["C"]
    bad_orbit = affine_orthogonal_bad_orbit(k) if k <= 3 else None
    rng = np.random.default_rng(seed)
    if exact:
        if order > 20:
            raise ValueError("exact local enumeration is intentionally capped at 20 bits")
        # Global sign commutes with the dynamics, so fix coordinate zero to +1.
        starts: Iterable[np.ndarray] = (
            spins_from_integer(state << 1, order) for state in range(1 << (order - 1))
        )
        sample_count = 1 << (order - 1)
    else:
        def random_starts() -> Iterable[np.ndarray]:
            for _ in range(samples):
                spin = rng.choice(np.array([-1, 1], dtype=np.int32), size=order)
                spin[0] = 1  # exact global-sign quotient
                yield spin

        starts = random_starts()
        sample_count = samples

    all_metrics = MetricAccumulator()
    bad_metrics = MetricAccumulator()
    nonbad_metrics = MetricAccumulator()
    energy_histogram: Counter[int] = Counter()
    flip_histogram: Counter[int] = Counter()
    bad_count = 0
    bad_projective_endpoint_histogram: Counter[int] = Counter()
    bad_energy_shell_count = 0
    anti_self_dual_plateaued_count = 0
    bent_count = 0
    bent_and_bad_count = 0
    near_bent_count = 0
    near_bent_and_bad_count = 0
    bad_energy = order * ((1 << k) // 2 + 1)
    explicit_state = negative_mask(data["bad_spin"])
    explicit_projective_state = min(
        explicit_state, ((1 << order) - 1) ^ explicit_state
    )
    for start in starts:
        terminal, energy, flips = best_improvement(c_matrix, start, -1)
        metrics = walsh_metrics(start, data, bad_orbit if k <= 3 else None)
        all_metrics.add(metrics)
        energy_histogram[energy] += 1
        flip_histogram[flips] += 1
        state = negative_mask(terminal)
        in_bad_orbit = bad_orbit is not None and state in bad_orbit
        terminal_walsh_spin = data["f"] * terminal
        terminal_spectrum = data["H"] @ terminal_walsh_spin
        anti_self_dual_plateaued = bool(
            np.count_nonzero(terminal_spectrum) == order // 4
            and np.all(
                (terminal_spectrum == 0)
                | (terminal_spectrum == -2 * (1 << k) * terminal_walsh_spin)
            )
        )
        anti_self_dual_plateaued_count += int(anti_self_dual_plateaued)
        bent = abs(metrics["max_abs_over_sqrt_m"] - 1.0) < 1e-12
        near_bent = metrics["max_abs_over_sqrt_m"] <= 1.5 + 1e-12
        bent_count += int(bent)
        near_bent_count += int(near_bent)
        if energy == bad_energy:
            bad_energy_shell_count += 1
        if in_bad_orbit:
            bad_count += 1
            full_mask = (1 << order) - 1
            bad_projective_endpoint_histogram[min(state, full_mask ^ state)] += 1
            bent_and_bad_count += int(bent)
            near_bent_and_bad_count += int(near_bent)
            bad_metrics.add(metrics)
        else:
            nonbad_metrics.add(metrics)

    return {
        "k": k,
        "m": order,
        "mode": "exact_global_sign_quotient" if exact else "monte_carlo_global_sign_quotient",
        "seed": None if exact else seed,
        "samples_after_quotient": sample_count,
        "affine_orthogonal_bad_orbit_size_including_global_sign": (
            len(bad_orbit) if bad_orbit is not None else None
        ),
        "bad_orbit_landings": bad_count if bad_orbit is not None else None,
        "bad_orbit_probability": (
            bad_count / sample_count if bad_orbit is not None else None
        ),
        "bad_orbit_zero_hit_upper_95": (
            1 - 0.05 ** (1 / sample_count)
            if bad_orbit is not None and bad_count == 0
            else None
        ),
        "bad_projective_endpoint_landing_histogram": {
            hex(key): value
            for key, value in sorted(bad_projective_endpoint_histogram.items())
        },
        "explicit_bad_point_projective_mask": hex(explicit_projective_state),
        "explicit_bad_point_basin_count": bad_projective_endpoint_histogram.get(
            explicit_projective_state, 0
        ),
        "explicit_bad_point_basin_probability": (
            bad_projective_endpoint_histogram.get(explicit_projective_state, 0)
            / sample_count
        ),
        "bad_energy_shell_probability": bad_energy_shell_count / sample_count,
        "anti_self_dual_2_plateaued_terminal_count": anti_self_dual_plateaued_count,
        "anti_self_dual_2_plateaued_terminal_probability": (
            anti_self_dual_plateaued_count / sample_count
        ),
        "bent_initial_count": bent_count,
        "bad_landings_among_bent_initial": bent_and_bad_count,
        "near_bent_initial_count_max_coefficient_at_most_1.5_sqrt_m": near_bent_count,
        "bad_landings_among_near_bent_initial": near_bent_and_bad_count,
        "terminal_energy_histogram": {
            str(key): value for key, value in sorted(energy_histogram.items())
        },
        "flip_count_histogram": {
            str(key): value for key, value in sorted(flip_histogram.items())
        },
        "initial_walsh_metrics_all": all_metrics.report(),
        "initial_walsh_metrics_bad_basin": bad_metrics.report(),
        "initial_walsh_metrics_other_basins": nonbad_metrics.report(),
    }


def run_regularizing_base(k: int, samples: int, seed: int) -> Dict[str, object]:
    """Separate one-block diagnostic: +C regularization followed by -C descent."""
    data = quadratic_data(k)
    matrix = data["C"]
    h_matrix = data["H"]
    f = data["f"]
    order = len(matrix)
    s = 1 << k
    rng = np.random.default_rng(seed)
    bad_orbit = affine_orthogonal_bad_orbit(k) if k <= 3 else None
    all_metrics = MetricAccumulator()
    positive_metrics = MetricAccumulator()
    linear_metrics = MetricAccumulator()
    positive_count = 0
    linear_count = 0
    bent_base_count = 0
    plateaued_terminal_count = 0
    bent_to_plateaued_count = 0
    linear_bent_to_plateaued_count = 0
    linear_bent_base_count = 0
    explicit_terminal_count = 0
    threshold = 0.05

    for _ in range(samples):
        start = rng.choice(np.array([-1, 1], dtype=np.int32), size=order)
        start[0] = 1
        base, base_energy, plus_flips = best_improvement(matrix, start, +1)
        local_margins = base * (matrix @ base)
        if np.any(local_margins < 0):
            raise AssertionError("+C regularizing base is not one-flip stable")
        terminal, descending_cap, minus_flips = best_improvement(matrix, base, -1)
        relative_terminal = base * terminal
        switched_target = -base[:, None] * matrix * base[None, :]
        weights = local_margins.astype(np.int64)
        selected = relative_terminal < 0
        outside = ~selected
        if np.any(selected):
            induced_sums = switched_target[np.ix_(selected, selected)].sum(axis=1)
            p_sum = int(weights[selected].sum() + 2 * induced_sums.sum())
        else:
            p_sum = 0
        outside_mass = int(weights[outside].sum())
        unpaid = max(0, outside_mass - p_sum)
        delta = max(0, base_energy - descending_cap)
        if delta != 2 * unpaid:
            raise AssertionError((delta, 2 * unpaid, base_energy, descending_cap))
        if unpaid:
            ordered = np.sort(weights[outside])[::-1]
            kappa = int(np.searchsorted(np.cumsum(ordered), unpaid, side="left") + 1)
        else:
            kappa = 0

        g = f * base
        z = f * terminal
        g_spectrum = h_matrix @ g
        z_spectrum = h_matrix @ z
        # Independent Walsh-coordinate check of the energy/defect formula.
        delta_walsh = max(0, int(g @ g_spectrum + z @ z_spectrum - 2 * order))
        if delta != delta_walsh:
            raise AssertionError((delta, delta_walsh))
        bent_base = bool(np.all(np.abs(g_spectrum) == s))
        plateaued_terminal = bool(
            np.count_nonzero(z_spectrum) == order // 4
            and np.all((z_spectrum == 0) | (z_spectrum == -2 * s * z))
        )
        explicit_terminal = (
            bad_orbit is not None and negative_mask(terminal) in bad_orbit
        )
        distance_orbit = bad_orbit if k <= 3 else None
        base_walsh = walsh_metrics(base, data, distance_orbit)
        terminal_walsh = walsh_metrics(terminal, data, distance_orbit)
        metrics: Dict[str, float] = {
            "defect_over_m32": delta / order**1.5,
            "kappa_over_m": kappa / order,
            "selected_fraction": float(selected.mean()),
            "outside_fraction": float(outside.mean()),
            "plus_flips_over_m": plus_flips / order,
            "minus_flips_over_m": minus_flips / order,
            "base_bent": float(bent_base),
            "terminal_anti_self_dual_2_plateaued": float(plateaued_terminal),
            "terminal_in_explicit_affine_orthogonal_orbit": float(explicit_terminal),
            "base_bent_residual": float(
                np.sum((np.abs(g_spectrum) - s) ** 2) / order**2
            ),
            "terminal_anti_2_plateaued_residual": float(
                np.sum(
                    np.minimum(
                        z_spectrum.astype(np.int64) ** 2,
                        (z_spectrum.astype(np.int64) + 2 * s * z) ** 2,
                    )
                )
                / order**2
            ),
        }
        metrics.update({f"base_{key}": value for key, value in base_walsh.items()})
        metrics.update(
            {f"terminal_{key}": value for key, value in terminal_walsh.items()}
        )
        all_metrics.add(metrics)
        if delta > 0:
            positive_count += 1
            positive_metrics.add(metrics)
        if kappa / order >= threshold:
            linear_count += 1
            linear_metrics.add(metrics)
        bent_base_count += int(bent_base)
        plateaued_terminal_count += int(plateaued_terminal)
        bent_to_plateaued_count += int(bent_base and plateaued_terminal)
        linear_bent_to_plateaued_count += int(
            kappa / order >= threshold and bent_base and plateaued_terminal
        )
        linear_bent_base_count += int(kappa / order >= threshold and bent_base)
        explicit_terminal_count += int(explicit_terminal)

    return {
        "k": k,
        "m": order,
        "mode": "+C largest-gain regularization then -C largest-gain descent",
        "seed": seed,
        "samples_after_global_sign_quotient": samples,
        "linear_core_threshold_kappa_over_m": threshold,
        "positive_defect_count": positive_count,
        "positive_defect_probability": positive_count / samples,
        "linear_core_count": linear_count,
        "linear_core_probability": linear_count / samples,
        "bent_base_count": bent_base_count,
        "anti_self_dual_2_plateaued_terminal_count": plateaued_terminal_count,
        "bent_base_to_anti_self_dual_2_plateaued_terminal_count": (
            bent_to_plateaued_count
        ),
        "linear_core_and_bent_to_plateaued_count": linear_bent_to_plateaued_count,
        "linear_core_and_bent_base_count": linear_bent_base_count,
        "explicit_bad_orbit_terminal_count": explicit_terminal_count,
        "all": all_metrics.report(),
        "positive_defect_only": positive_metrics.report(),
        "linear_core_only": linear_metrics.report(),
    }


def sylvester(order: int) -> np.ndarray:
    values = np.arange(order, dtype=np.int64)
    parity = parity_table(order)
    return (1 - 2 * parity[np.bitwise_and(values[:, None], values[None, :])]).astype(
        np.int32
    )


def full_obstruction_matrix(k: int) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    data = quadratic_data(k)
    c_matrix = data["C"]
    order = len(c_matrix)
    l_matrix = np.array([[1, -1], [-1, 1]], dtype=np.int32)
    bridge = np.kron(sylvester(order // 2), l_matrix).astype(np.int32)
    assert np.array_equal(bridge @ np.ones(order, dtype=np.int32), np.zeros(order))
    full = np.block([[c_matrix, bridge], [bridge.T, -c_matrix]]).astype(np.int32)
    assert np.all(np.diag(full) == 0)
    assert np.all(np.abs(full + np.eye(2 * order, dtype=np.int32)) == 1)
    data = dict(data)
    data["bridge"] = bridge
    return full, data


def augmented_branch(
    matrix: np.ndarray,
    switched_field: np.ndarray,
    input_spin: np.ndarray,
    anchor: np.ndarray,
    anchor_energy: int,
    free_energy: int,
) -> Dict[str, object]:
    free = ~anchor
    anchor_vertices = np.flatnonzero(anchor)
    free_vertices = np.flatnonzero(free)
    x_free = input_spin[free_vertices]
    x_anchor = input_spin[anchor_vertices]
    free_block = (
        x_free[:, None]
        * matrix[np.ix_(free_vertices, free_vertices)]
        * x_free[None, :]
    ).astype(np.int32)
    cross = (
        x_free
        * (matrix[np.ix_(free_vertices, anchor_vertices)] @ x_anchor)
    ).astype(np.int32)
    sigma = 1 if anchor_energy > 0 else -1
    augmented = np.zeros((len(free_vertices) + 1,) * 2, dtype=np.int32)
    augmented[:-1, :-1] = free_block
    augmented[:-1, -1] = cross
    augmented[-1, :-1] = cross
    initial_free = np.where(
        sigma * cross > 0, 1, np.where(sigma * cross < 0, -1, 1)
    ).astype(np.int32)
    initial = np.concatenate((initial_free, np.ones(1, dtype=np.int32)))
    terminal, cap, flips = best_improvement(augmented, initial, sigma)
    effective = (terminal[-1] * terminal[:-1]).astype(np.int32)
    target_block = sigma * free_block
    weights = (-sigma * switched_field[free_vertices]).astype(np.int64)
    if np.any(weights < 0):
        raise AssertionError("hard-branch free weights have wrong sign")
    selected = effective < 0
    outside = ~selected
    if np.any(selected):
        induced_sums = target_block[np.ix_(selected, selected)].sum(axis=1)
        p_sum = int(weights[selected].sum() + 2 * induced_sums.sum())
    else:
        p_sum = 0
    outside_mass = int(weights[outside].sum())
    unpaid = max(0, outside_mass - p_sum)
    defect_formula = 2 * unpaid
    defect = max(0, abs(free_energy) - cap)
    if defect != defect_formula:
        raise AssertionError(
            (defect, defect_formula, free_energy, cap, outside_mass, p_sum)
        )
    if unpaid:
        ordered = np.sort(weights[outside])[::-1]
        kappa = int(np.searchsorted(np.cumsum(ordered), unpaid, side="left") + 1)
    else:
        kappa = 0
    return {
        "cap": cap,
        "defect": defect,
        "flips": flips,
        "kappa": kappa,
        "selected_size": int(selected.sum()),
        "outside_size": int(outside.sum()),
        "unpaid_mass": unpaid,
        "terminal_effective": effective,
        "free_vertices": free_vertices,
    }


def full_sample(
    matrix: np.ndarray,
    data: Dict[str, np.ndarray],
    input_spin: np.ndarray,
    bad_orbit: Optional[Set[int]],
) -> Dict[str, object]:
    order = len(matrix)
    half = order // 2
    field = matrix @ input_spin
    switched_field = input_spin * field
    agreement = switched_field >= 0
    x_i = input_spin[agreement]
    x_j = input_spin[~agreement]
    p = int(x_i @ matrix[np.ix_(agreement, agreement)] @ x_i)
    r = int(x_j @ matrix[np.ix_(~agreement, ~agreement)] @ x_j)
    response = p - r
    if p * r >= 0:
        branch_metrics = {
            "defect_over_N32": 0.0,
            "kappa_over_N": 0.0,
            "kappa_branch_selected_over_N": 0.0,
            "kappa_branch_outside_over_N": 0.0,
            "defect_branch_selected_over_N": 0.0,
            "defect_branch_outside_over_N": 0.0,
            "hard_branch": 0.0,
            "explicit_native_bad_core": 0.0,
        }
        return {"metrics": branch_metrics, "branches": None, "response": response}
    if not (p > 0 > r):
        raise AssertionError((p, r))
    branch_j = augmented_branch(
        matrix, switched_field, input_spin, agreement, p, r
    )
    branch_i = augmented_branch(
        matrix, switched_field, input_spin, ~agreement, r, p
    )
    branches = [branch_j, branch_i]
    chosen_defect = min(
        branches, key=lambda item: (int(item["defect"]), int(item["kappa"]))
    )
    chosen_kappa = min(
        branches, key=lambda item: (int(item["kappa"]), int(item["defect"]))
    )
    defect = min(int(branch_j["defect"]), int(branch_i["defect"]))
    if defect != max(0, response - max(p + int(branch_j["cap"]), -r + int(branch_i["cap"]))):
        raise AssertionError("combined certificate identity failed")

    explicit = False
    if bad_orbit is not None:
        native_top = set(range(half))
        native_bottom = set(range(half, order))
        for branch in branches:
            vertices = [int(v) for v in branch["free_vertices"]]
            vertex_set = set(vertices)
            if vertex_set == native_top or vertex_set == native_bottom:
                effective = np.asarray(branch["terminal_effective"])
                # ``effective`` is the terminal spin relative to the row-sign
                # input.  The affine-O(q) orbit lives in the original
                # physical C coordinates, so undo that switch before testing.
                physical_terminal = input_spin[vertices] * effective
                explicit = explicit or negative_mask(physical_terminal) in bad_orbit

    metrics = {
        "defect_over_N32": defect / (order**1.5),
        # The diffuse-core theorem uses kappa_*=min(kappa_I,kappa_J),
        # independently of which branch minimizes the numerical defect.
        "kappa_over_N": int(chosen_kappa["kappa"]) / order,
        "kappa_branch_selected_over_N": int(chosen_kappa["selected_size"]) / order,
        "kappa_branch_outside_over_N": int(chosen_kappa["outside_size"]) / order,
        "defect_branch_selected_over_N": int(chosen_defect["selected_size"]) / order,
        "defect_branch_outside_over_N": int(chosen_defect["outside_size"]) / order,
        "hard_branch": 1.0,
        "explicit_native_bad_core": float(explicit),
    }
    return {"metrics": metrics, "branches": branches, "response": response}


def run_full(
    k: int,
    samples: int,
    seed: int,
    diagonal_slice_exact: bool,
) -> Dict[str, object]:
    matrix, data = full_obstruction_matrix(k)
    half = len(matrix) // 2
    bad_orbit = affine_orthogonal_bad_orbit(k) if k <= 3 else None
    rng = np.random.default_rng(seed)
    if diagonal_slice_exact:
        if half > 20:
            raise ValueError("exact diagonal slice is intentionally capped at 20 bits")
        iterator: Iterable[np.ndarray] = (
            np.concatenate((x, x))
            for x in (
                spins_from_integer(state << 1, half)
                for state in range(1 << (half - 1))
            )
        )
        sample_count = 1 << (half - 1)
        mode = "exact_X2_equals_X1_slice_global_sign_quotient"
    else:
        def random_inputs() -> Iterable[np.ndarray]:
            for _ in range(samples):
                spin = rng.choice(np.array([-1, 1], dtype=np.int32), size=2 * half)
                spin[0] = 1
                yield spin

        iterator = random_inputs()
        sample_count = samples
        mode = "uniform_monte_carlo_global_sign_quotient"

    all_metrics = MetricAccumulator()
    hard_metrics = MetricAccumulator()
    diffuse_metrics = MetricAccumulator()
    nondiffuse_metrics = MetricAccumulator()
    bent_input_metrics = MetricAccumulator()
    explicit_count = 0
    diffuse_count = 0
    positive_defect_count = 0
    bent_count = 0
    diffuse_and_bent_count = 0
    near_bent_count = 0
    diffuse_and_near_bent_count = 0
    kappa_histogram: Counter[int] = Counter()
    threshold = 0.05
    for input_spin in iterator:
        sample = full_sample(matrix, data, input_spin, bad_orbit)
        metrics = dict(sample["metrics"])
        block_metrics = [
            walsh_metrics(input_spin[:half], data, bad_orbit),
            walsh_metrics(input_spin[half:], data, bad_orbit),
        ]
        for key in block_metrics[0]:
            metrics[f"input_block_mean_{key}"] = 0.5 * (
                block_metrics[0][key] + block_metrics[1][key]
            )
        all_metrics.add(metrics)
        if metrics["hard_branch"]:
            hard_metrics.add(metrics)
        if metrics["kappa_over_N"] >= threshold:
            diffuse_count += 1
            diffuse_metrics.add(metrics)
        else:
            nondiffuse_metrics.add(metrics)
        bent = all(
            abs(block["max_abs_over_sqrt_m"] - 1.0) < 1e-12
            for block in block_metrics
        )
        near_bent = all(
            block["max_abs_over_sqrt_m"] <= 1.5 + 1e-12
            for block in block_metrics
        )
        if bent:
            bent_input_metrics.add(metrics)
        bent_count += int(bent)
        near_bent_count += int(near_bent)
        diffuse_and_bent_count += int(bent and metrics["kappa_over_N"] >= threshold)
        diffuse_and_near_bent_count += int(
            near_bent and metrics["kappa_over_N"] >= threshold
        )
        positive_defect_count += int(metrics["defect_over_N32"] > 0)
        kappa_histogram[int(round(metrics["kappa_over_N"] * (2 * half)))] += 1
        explicit_count += int(metrics["explicit_native_bad_core"])

    return {
        "k": k,
        "m_per_block": half,
        "N": 2 * half,
        "mode": mode,
        "seed": None if diagonal_slice_exact else seed,
        "samples_after_global_sign_quotient": sample_count,
        "diffuse_threshold_kappa_over_N": threshold,
        "diffuse_count": diffuse_count,
        "diffuse_probability": diffuse_count / sample_count,
        "positive_defect_count": positive_defect_count,
        "positive_defect_probability": positive_defect_count / sample_count,
        "both_blocks_bent_count": bent_count,
        "diffuse_and_both_blocks_bent_count": diffuse_and_bent_count,
        "both_blocks_near_bent_count": near_bent_count,
        "diffuse_and_both_blocks_near_bent_count": diffuse_and_near_bent_count,
        "kappa_histogram": {
            str(key): value for key, value in sorted(kappa_histogram.items())
        },
        "explicit_native_bad_core_count": explicit_count,
        "explicit_native_bad_core_probability": explicit_count / sample_count,
        "all": all_metrics.report(),
        "hard_branch_only": hard_metrics.report(),
        "diffuse_only": diffuse_metrics.report(),
        "nondiffuse_only": nondiffuse_metrics.report(),
        "both_input_blocks_bent": bent_input_metrics.report(),
    }


def fwht(values: np.ndarray) -> np.ndarray:
    """Unnormalized standard Walsh transform, exact int64 arithmetic."""
    transformed = values.astype(np.int64, copy=True)
    width = 1
    while width < len(transformed):
        blocks = transformed.reshape(-1, 2 * width)
        left = blocks[:, :width].copy()
        right = blocks[:, width:].copy()
        blocks[:, :width] = left + right
        blocks[:, width:] = left - right
        width *= 2
    return transformed


def symplectic_walsh(values: np.ndarray, k: int) -> np.ndarray:
    standard = fwht(values)
    mask = (1 << k) - 1
    indices = np.arange(len(values), dtype=np.int64)
    swapped = ((indices & mask) << k) | (indices >> k)
    return standard[swapped]


def structured_seed(k: int, base: np.ndarray) -> np.ndarray:
    """base(a0,a1,b0,b1) times f on the remaining coordinate pairs."""
    order = 1 << (2 * k)
    values = np.arange(order, dtype=np.int64)
    mask = (1 << k) - 1
    a, b = values & mask, values >> k
    base_index = (a & 3) | ((b & 3) << 2)
    tail_parity = np.fromiter(
        (bin(int(x)).count("1") & 1 for x in ((a >> 2) & (b >> 2))),
        dtype=np.int8,
        count=order,
    )
    tail_f = 1 - 2 * tail_parity
    return (base[base_index] * tail_f).astype(np.int64)


def walsh_best_improvement(
    initial_g: np.ndarray, k: int, target_sign: int
) -> Tuple[np.ndarray, np.ndarray, int, int]:
    """Best improvement for C=K-I, maintained in bent-Walsh coordinates."""
    g = initial_g.astype(np.int64, copy=True)
    hg = symplectic_walsh(g, k)
    order = len(g)
    s = 1 << k
    values = np.arange(order, dtype=np.int64)
    a, b = values & (s - 1), values >> k
    parity = parity_table(s)
    flips = 0
    while True:
        margins = g * hg - 1
        gains = -4 * target_sign * margins
        coordinate = int(np.argmax(gains))
        if int(gains[coordinate]) <= 0:
            energy = int(target_sign * (g @ hg - order))
            if np.any(target_sign * margins < 0):
                raise AssertionError("Walsh ascent terminal condition failed")
            if not np.array_equal(hg, symplectic_walsh(g, k)):
                raise AssertionError("incremental Walsh response drifted")
            return g, hg, energy, flips
        old = int(g[coordinate])
        g[coordinate] = -old
        c = coordinate & (s - 1)
        d = coordinate >> k
        beta_parity = parity[np.bitwise_xor(a & d, b & c)]
        column = 1 - 2 * beta_parity.astype(np.int64)
        hg -= 2 * old * column
        flips += 1


def core_from_walsh_pair(
    base_g: np.ndarray,
    base_hg: np.ndarray,
    terminal_g: np.ndarray,
    terminal_hg: np.ndarray,
) -> Dict[str, object]:
    order = len(base_g)
    weights = base_g * base_hg - 1
    if np.any(weights < 0):
        raise AssertionError("structured base is not regularizing")
    relative = base_g * terminal_g
    outside = relative > 0
    base_energy = int(base_g @ base_hg - order)
    terminal_target_energy = int(order - terminal_g @ terminal_hg)
    delta = max(0, base_energy - terminal_target_energy)
    if delta % 2:
        raise AssertionError("unpaid mass is not integral")
    unpaid = delta // 2
    if unpaid:
        ordered = np.sort(weights[outside])[::-1]
        kappa = int(np.searchsorted(np.cumsum(ordered), unpaid, side="left") + 1)
    else:
        kappa = 0
    return {
        "base_energy": base_energy,
        "terminal_target_energy": terminal_target_energy,
        "delta": delta,
        "delta_over_m32": delta / order**1.5,
        "unpaid_mass": unpaid,
        "kappa": kappa,
        "kappa_over_m": kappa / order,
        "selected_size": int(np.count_nonzero(relative < 0)),
        "outside_size": int(np.count_nonzero(outside)),
        "terminal_walsh_support": int(np.count_nonzero(terminal_hg)),
        "terminal_walsh_distinct_absolute_levels": sorted(
            int(value) for value in np.unique(np.abs(terminal_hg))
        ),
    }


def structured_counterfamily_audit(max_trajectory_k: int = 8) -> Dict[str, object]:
    """Verify the exact non-flat pair and its actual prescribed trajectory."""
    z_star = np.array(
        [1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1],
        dtype=np.int64,
    )
    w = np.array([1, -1, -1, -1] + [1] * 12, dtype=np.int64)
    v = np.array([-1] + [1] * 15, dtype=np.int64)
    exact_pairs: List[Dict[str, object]] = []
    trajectories: List[Dict[str, object]] = []
    for k in range(2, max_trajectory_k + 1):
        order = 1 << (2 * k)
        s = 1 << k
        z = structured_seed(k, z_star)
        g = structured_seed(k, w)
        comparison_v = structured_seed(k, v)
        hz = symplectic_walsh(z, k)
        hg = symplectic_walsh(g, k)
        h_comparison_v = symplectic_walsh(comparison_v, k)
        z_products = z * hz
        g_products = g * hg
        if set(np.unique(z_products)) != {-7 * s // 2, -s // 2}:
            raise AssertionError("z-star tensor product formula failed")
        if set(np.unique(g_products)) != {s // 2, 3 * s // 2, 5 * s // 2}:
            raise AssertionError("w tensor product formula failed")
        static_core = core_from_walsh_pair(g, hg, z, hz)
        expected_delta = order * (s // 8 - 2) if s > 16 else 0
        expected_kappa = (
            math.ceil(order * (s - 16) / (8 * (5 * s - 2))) if s > 16 else 0
        )
        if static_core["delta"] != expected_delta or static_core["kappa"] != expected_kappa:
            raise AssertionError((k, static_core, expected_delta, expected_kappa))
        comparison_core = core_from_walsh_pair(
            g, hg, comparison_v, h_comparison_v
        )
        expected_comparison_kappa = (
            math.ceil(order * (s - 16) / (8 * (s - 2))) if s > 16 else 0
        )
        if (
            comparison_core["delta"] != expected_delta
            or comparison_core["kappa"] != expected_comparison_kappa
        ):
            raise AssertionError(
                (
                    "comparison-v",
                    k,
                    comparison_core,
                    expected_delta,
                    expected_comparison_kappa,
                )
            )
        exact_pairs.append(
            {
                "k": k,
                "m": order,
                "base_signed_walsh_products": {
                    str(int(value)): int(np.count_nonzero(g_products == value))
                    for value in np.unique(g_products)
                },
                "terminal_signed_walsh_products": {
                    str(int(value)): int(np.count_nonzero(z_products == value))
                    for value in np.unique(z_products)
                },
                **static_core,
            }
        )

        trajectory_g, trajectory_hg, trajectory_cap, flips = walsh_best_improvement(
            g, k, -1
        )
        trajectory_core = core_from_walsh_pair(
            g, hg, trajectory_g, trajectory_hg
        )
        if trajectory_core["terminal_target_energy"] != trajectory_cap:
            raise AssertionError("trajectory cap mismatch")
        coordinates = np.arange(order, dtype=np.int64)
        mask = s - 1
        a, b = coordinates & mask, coordinates >> k
        base_index = (a & 3) | ((b & 3) << 2)
        selected = g * trajectory_g < 0
        deviation_from_v = trajectory_g != comparison_v
        selected_by_base_fibre = [
            int(np.count_nonzero(selected & (base_index == fibre)))
            for fibre in range(16)
        ]
        deviation_from_v_by_base_fibre = [
            int(np.count_nonzero(deviation_from_v & (base_index == fibre)))
            for fibre in range(16)
        ]
        if sum(selected_by_base_fibre) != trajectory_core["selected_size"]:
            raise AssertionError("selected fibre counts do not sum")
        if sum(deviation_from_v_by_base_fibre) != int(
            np.count_nonzero(deviation_from_v)
        ):
            raise AssertionError("comparison-v fibre counts do not sum")
        trajectories.append(
            {
                "k": k,
                "m": order,
                "flip_count": flips,
                "distance_to_comparison_v": int(
                    np.count_nonzero(deviation_from_v)
                ),
                "selected_by_base_fibre": selected_by_base_fibre,
                "deviation_from_comparison_v_by_base_fibre": (
                    deviation_from_v_by_base_fibre
                ),
                "comparison_v_static_core": comparison_core,
                **trajectory_core,
            }
        )
    return {
        "seed_coordinate_order": "j=a0+2*a1+4*b0+8*b1, then f on tail pairs",
        "z_star": z_star.tolist(),
        "w": w.tolist(),
        "v": v.tolist(),
        "exact_nonflat_pair": exact_pairs,
        "prescribed_descent_from_nonbent_regularizing_seed": trajectories,
    }


def canonical_payload_hash(payload: Dict[str, object]) -> str:
    cleaned = dict(payload)
    cleaned.pop("canonical_payload_sha256", None)
    encoded = json.dumps(cleaned, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("computations/results/regular_hadamard_walsh_basins.json"),
    )
    parser.add_argument("--seed", type=int, default=20260814)
    parser.add_argument("--local-k3-samples", type=int, default=50000)
    parser.add_argument("--local-k4-samples", type=int, default=10000)
    parser.add_argument("--full-k2-samples", type=int, default=20000)
    parser.add_argument("--full-k3-samples", type=int, default=3000)
    parser.add_argument("--regularizing-k3-samples", type=int, default=10000)
    parser.add_argument("--regularizing-k4-samples", type=int, default=3000)
    parser.add_argument(
        "--structured-max-k",
        type=int,
        default=8,
        help="largest exact structured-seed trajectory order (9 takes several minutes)",
    )
    args = parser.parse_args()

    payload: Dict[str, object] = {
        "schema": "regular-hadamard-walsh-basin-audit-v1",
        "normalization": (
            "Q(A)=max_z |z^T A z|; full two-block defects use N^(3/2), "
            "one-block and structured defects use m^(3/2)"
        ),
        "dynamics": (
            "largest strictly positive one-flip gain, least coordinate on ties; "
            "augmented zero cross-fields start at switched spin +1"
        ),
        "local_experiment": [
            run_local(2, 0, args.seed + 2, True),
            run_local(3, args.local_k3_samples, args.seed + 3, False),
            run_local(4, args.local_k4_samples, args.seed + 4, False),
        ],
        "symmetry_invariant_tie_control": k2_symmetry_invariant_tie_audit(),
        "generic_regularizing_base_experiment": [
            run_regularizing_base(3, args.regularizing_k3_samples, args.seed + 33),
            run_regularizing_base(4, args.regularizing_k4_samples, args.seed + 34),
        ],
        "structured_nonflat_counterfamily": structured_counterfamily_audit(
            args.structured_max_k
        ),
        "full_row_sign_experiment": [
            run_full(2, 0, args.seed + 12, True),
            run_full(2, args.full_k2_samples, args.seed + 22, False),
            run_full(3, args.full_k3_samples, args.seed + 23, False),
        ],
    }
    payload["canonical_payload_sha256"] = canonical_payload_hash(payload)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
