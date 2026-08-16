#!/usr/bin/env python3
"""Independent finite audit of vector BLR and the response constant 11.

The checks are deliberately independent of the originating synchronization
script.  No external solver is needed.
"""

from __future__ import annotations

import argparse
import itertools
import json
from collections import Counter, deque
from pathlib import Path


def popcount(x: int) -> int:
    return bin(x).count("1")


def xor_sum(values) -> int:
    total = 0
    for value in values:
        total ^= value
    return total


def bit_at(table: int, x: int) -> int:
    """Truth-table bit; the value at zero is fixed to zero."""

    return 0 if x == 0 else (table >> (x - 1)) & 1


def linear_tables(k: int) -> list[int]:
    n = 1 << k
    return [
        sum(((popcount(alpha & x) & 1) << (x - 1)) for x in range(1, n))
        for alpha in range(n)
    ]


def best_fourier_character(table: int, k: int) -> tuple[int, int]:
    """Return (alpha, correlation numerator) with deterministic tie-breaking."""

    n = 1 << k
    best_alpha = 0
    best_correlation = -n - 1
    for alpha in range(n):
        correlation = sum(
            1
            if bit_at(table, x) == (popcount(alpha & x) & 1)
            else -1
            for x in range(n)
        )
        if correlation > best_correlation:
            best_alpha = alpha
            best_correlation = correlation
    return best_alpha, best_correlation


def scalar_blr_audit(max_k: int) -> dict[str, object]:
    results: dict[str, object] = {}
    for k in range(1, max_k + 1):
        n = 1 << k
        tables = 1 << (n - 1)
        worst_ratio = 0.0
        witness = None
        profiles: Counter[tuple[int, int]] = Counter()
        for table in range(tables):
            rejects = sum(
                bit_at(table, a) ^ bit_at(table, b) ^ bit_at(table, a ^ b)
                for a in range(n)
                for b in range(n)
            )
            alpha, correlation = best_fourier_character(table, k)
            distance = (n - correlation) // 2
            if distance * n > rejects:
                raise AssertionError(
                    f"scalar BLR failed: k={k}, table={table}, "
                    f"distance={distance}, rejects={rejects}"
                )
            profiles[(distance, rejects)] += 1
            if rejects:
                ratio = distance * n / rejects
                if ratio > worst_ratio:
                    worst_ratio = ratio
                    witness = {
                        "table": table,
                        "nearest_linear_alpha": alpha,
                        "distance_count": distance,
                        "rejection_count_ordered_pairs": rejects,
                    }
        results[str(k)] = {
            "truth_tables_checked": tables,
            "worst_distance_over_rejection_ratio": worst_ratio,
            "worst_witness": witness,
            "distance_rejection_profiles": {
                f"{distance},{rejects}": count
                for (distance, rejects), count in sorted(profiles.items())
            },
        }
    return results


def scalar_coset_representatives(k: int) -> list[int]:
    n = 1 << k
    linears = linear_tables(k)
    seen: set[int] = set()
    representatives: list[int] = []
    for table in range(1 << (n - 1)):
        if table in seen:
            continue
        coset = {table ^ linear for linear in linears}
        seen.update(coset)
        representatives.append(min(coset))
    return representatives


def joint_average_sharpness(max_even_k: int = 8) -> dict[str, object]:
    """Check the paired bent construction making the average constant sharp."""

    results: dict[str, object] = {}
    for k in range(2, max_even_k + 1, 2):
        n = 1 << k

        def phi(x: int) -> int:
            return sum(
                ((x >> (2 * i)) & 1) * ((x >> (2 * i + 1)) & 1)
                for i in range(k // 2)
            ) & 1

        first = sum(phi(x) << (x - 1) for x in range(1, n))
        second = sum(((1 ^ phi(x)) << (x - 1)) for x in range(1, n))
        maximum_defect = max(
            (
                bit_at(first, a) ^ bit_at(first, b) ^ bit_at(first, a ^ b)
            )
            + (
                bit_at(second, a)
                ^ bit_at(second, b)
                ^ bit_at(second, a ^ b)
            )
            for a in range(n)
            for b in range(n)
        )
        distances = []
        for table in (first, second):
            _, correlation = best_fourier_character(table, k)
            distances.append((n - correlation) // 2)
        expected = [
            (n - (1 << (k // 2))) // 2,
            (n - (1 << (k // 2))) // 2 - 1,
        ]
        if maximum_defect != 1 or distances != expected:
            raise AssertionError("paired bent sharpness formula failed")
        results[str(k)] = {
            "pointwise_defect": maximum_defect,
            "coordinate_distances_to_linear": distances,
            "minimum_average_vector_error": sum(distances) / n,
            "formula": 1 - 2 ** (-k / 2) - 2 ** (-k),
        }
    return results


def vector_map_audit(k: int, output_dimension: int) -> dict[str, object]:
    n = 1 << k
    representatives = scalar_coset_representatives(k)
    linears = linear_tables(k)
    checked = 0
    worst_average_ratio = 0.0
    worst_uniform_ratio = 0.0
    worst_optimal_uniform_ratio = 0.0
    witnesses: dict[str, object] = {}

    defect_tables = {
        table: [
            bit_at(table, a) ^ bit_at(table, b) ^ bit_at(table, a ^ b)
            for a in range(n)
            for b in range(n)
        ]
        for table in representatives
    }
    error_tables = {
        table: [
            [bit_at(table ^ linear, x) for x in range(n)]
            for linear in linears
        ]
        for table in representatives
    }

    for coordinate_tables in itertools.product(
        representatives, repeat=output_dimension
    ):
        checked += 1
        delta = max(
            sum(defect_tables[table][pair] for table in coordinate_tables)
            for pair in range(n * n)
        )
        if delta == 0:
            continue

        chosen = [best_fourier_character(table, k)[0] for table in coordinate_tables]
        errors = [
            sum(
                error_tables[table][chosen[j]][x]
                for j, table in enumerate(coordinate_tables)
            )
            for x in range(n)
        ]
        average_numerator = sum(errors)
        uniform_error = max(errors)
        if average_numerator > delta * n or uniform_error > 3 * delta:
            raise AssertionError(
                f"vector BLR failed: k={k}, D={output_dimension}, "
                f"tables={coordinate_tables}, delta={delta}, errors={errors}"
            )

        optimal_uniform = min(
            max(
                sum(
                    error_tables[table][alphas[j]][x]
                    for j, table in enumerate(coordinate_tables)
                )
                for x in range(n)
            )
            for alphas in itertools.product(range(n), repeat=output_dimension)
        )

        average_ratio = average_numerator / (n * delta)
        uniform_ratio = uniform_error / delta
        optimal_ratio = optimal_uniform / delta
        for name, ratio, payload in (
            (
                "average",
                average_ratio,
                {"tables": coordinate_tables, "delta": delta, "errors": errors},
            ),
            (
                "uniform_for_fourier_choice",
                uniform_ratio,
                {"tables": coordinate_tables, "delta": delta, "errors": errors},
            ),
            (
                "optimal_uniform",
                optimal_ratio,
                {
                    "tables": coordinate_tables,
                    "delta": delta,
                    "optimal_uniform_error": optimal_uniform,
                },
            ),
        ):
            current = {
                "average": worst_average_ratio,
                "uniform_for_fourier_choice": worst_uniform_ratio,
                "optimal_uniform": worst_optimal_uniform_ratio,
            }[name]
            if ratio > current:
                witnesses[name] = payload
                if name == "average":
                    worst_average_ratio = ratio
                elif name == "uniform_for_fourier_choice":
                    worst_uniform_ratio = ratio
                else:
                    worst_optimal_uniform_ratio = ratio

    return {
        "k": k,
        "D": output_dimension,
        "coordinate_coset_classes": len(representatives),
        "vector_maps_modulo_coordinatewise_linear_addition": checked,
        "worst_average_ratio": worst_average_ratio,
        "worst_uniform_ratio_for_fourier_choice": worst_uniform_ratio,
        "worst_optimal_uniform_ratio": worst_optimal_uniform_ratio,
        "witnesses": witnesses,
    }


def cycle_contracting_k3_audit() -> dict[str, object]:
    """Exhaust all nonlinear coordinate multisets allowed by all Q=F_2^3 cycles."""

    k = 3
    n = 1 << k
    nonzero = tuple(range(1, n))
    linears = linear_tables(k)
    representatives = [x for x in scalar_coset_representatives(k) if x != 0]
    cycles: list[tuple[int, tuple[int, ...]]] = []
    for mask in range(1 << len(nonzero)):
        selected = tuple(
            nonzero[i] for i in range(len(nonzero)) if (mask >> i) & 1
        )
        if selected and xor_sum(selected) == 0:
            cycles.append((len(selected), selected))

    cycle_defects = [
        [sum(bit_at(table, q) for q in selected) & 1 for _, selected in cycles]
        for table in representatives
    ]
    triple_indices = [i for i, (size, _) in enumerate(cycles) if size == 3]
    minimum_nonzero_triple_weight = min(
        sum(defects[i] for i in triple_indices) for defects in cycle_defects
    )
    if minimum_nonzero_triple_weight != 3:
        raise AssertionError("unexpected scalar quotient-code minimum weight")

    # Seven triple constraints of capacity three imply at most seven nonlinear
    # coordinates.  Coordinates in the zero class are linear and can be fitted
    # exactly, so omitting them loses no maps relevant to the optimum.
    maximum_nonlinear_coordinates = 7 * 3 // minimum_nonzero_triple_weight
    error_tables = [
        [[bit_at(table ^ linear, x) for x in range(n)] for linear in linears]
        for table in representatives
    ]

    admissible = 0
    largest_optimal_uniform = 0
    witness = None
    by_size: Counter[int] = Counter()

    def feasible_with_radius(indices: tuple[int, ...], radius: int) -> bool:
        loads = [0] * n
        order = sorted(indices, key=lambda i: -sum(cycle_defects[i]))

        def recurse(position: int) -> bool:
            if position == len(order):
                return True
            coordinate = order[position]
            seen_rows: set[tuple[int, ...]] = set()
            for row in error_tables[coordinate]:
                row_tuple = tuple(row)
                if row_tuple in seen_rows:
                    continue
                seen_rows.add(row_tuple)
                if all(loads[x] + row[x] <= radius for x in range(n)):
                    for x in range(n):
                        loads[x] += row[x]
                    if recurse(position + 1):
                        return True
                    for x in range(n):
                        loads[x] -= row[x]
            return False

        return recurse(0)

    for coordinate_count in range(1, maximum_nonlinear_coordinates + 1):
        for indices in itertools.combinations_with_replacement(
            range(len(representatives)), coordinate_count
        ):
            if any(
                sum(cycle_defects[index][cycle] for index in indices)
                > cycles[cycle][0]
                for cycle in range(len(cycles))
            ):
                continue
            admissible += 1
            by_size[coordinate_count] += 1
            optimal_uniform = next(
                radius
                for radius in range(coordinate_count + 1)
                if feasible_with_radius(indices, radius)
            )
            if optimal_uniform > largest_optimal_uniform:
                largest_optimal_uniform = optimal_uniform
                witness = {
                    "coordinate_class_indices": indices,
                    "coordinate_tables": [representatives[i] for i in indices],
                }

    return {
        "nonlinear_scalar_classes": len(representatives),
        "minimum_nonzero_triple_syndrome_weight": minimum_nonzero_triple_weight,
        "maximum_nonlinear_coordinates_for_a_cycle_contractive_map": (
            maximum_nonlinear_coordinates
        ),
        "admissible_multisets": admissible,
        "admissible_by_coordinate_count": dict(sorted(by_size.items())),
        "largest_optimal_uniform_distance": largest_optimal_uniform,
        "witness": witness,
    }


def rank_binary(vectors, dimension: int) -> int:
    pivots = [0] * dimension
    rank = 0
    for vector in vectors:
        x = vector
        while x:
            pivot = x.bit_length() - 1
            if pivots[pivot]:
                x ^= pivots[pivot]
            else:
                pivots[pivot] = x
                rank += 1
                break
    return rank


def span(vectors) -> set[int]:
    result = {0}
    for vector in vectors:
        result |= {x ^ vector for x in tuple(result)}
    return result


def coordinate_lookup(basis) -> dict[int, int]:
    return {
        xor_sum(b for i, b in enumerate(basis) if (mask >> i) & 1): mask
        for mask in range(1 << len(basis))
    }


def word_lengths(support, dimension: int) -> list[int]:
    distances = [-1] * (1 << dimension)
    distances[0] = 0
    queue = deque([0])
    while queue:
        x = queue.popleft()
        for generator in support:
            y = x ^ generator
            if distances[y] < 0:
                distances[y] = distances[x] + 1
                queue.append(y)
    return distances


def blr_graph_section(
    basis: tuple[int, ...],
    quotient_representatives: tuple[int, ...],
    transversal: dict[int, int],
) -> tuple[dict[int, int], int, int]:
    d = len(basis)
    k = len(quotient_representatives)
    n = 1 << k
    full_basis = basis + quotient_representatives
    coordinates = coordinate_lookup(full_basis)

    coordinate_tables: list[int] = []
    for output_coordinate in range(d):
        table = 0
        for q in range(1, n):
            value = (coordinates[transversal[q]] >> output_coordinate) & 1
            table |= value << (q - 1)
        coordinate_tables.append(table)

    alphas = [best_fourier_character(table, k)[0] for table in coordinate_tables]
    section: dict[int, int] = {0: 0}
    selected_errors: list[int] = []
    for q in range(1, n):
        lift = xor_sum(
            quotient_representatives[i]
            for i in range(k)
            if (q >> i) & 1
        )
        correction_mask = sum(
            ((popcount(alpha & q) & 1) << j) for j, alpha in enumerate(alphas)
        )
        correction = xor_sum(
            basis[j] for j in range(d) if (correction_mask >> j) & 1
        )
        section[q] = lift ^ correction
        selected_errors.append(
            popcount(coordinates[transversal[q] ^ section[q]])
        )
    return section, sum(selected_errors), max(selected_errors, default=0)


def quotient_shortest_words(labels: tuple[int, ...], k: int) -> tuple[list[int], list[tuple[int, ...]]]:
    """Cayley distances and one deterministic shortest label word."""

    distances = [-1] * (1 << k)
    words: list[tuple[int, ...]] = [tuple() for _ in range(1 << k)]
    distances[0] = 0
    queue = deque([0])
    while queue:
        q = queue.popleft()
        for label in labels:
            target = q ^ label
            if distances[target] < 0:
                distances[target] = distances[q] + 1
                words[target] = words[q] + (label,)
                queue.append(target)
    if any(distance < 0 for distance in distances):
        raise AssertionError("quotient labels do not span")
    return distances, words


def partial_graph_section(
    basis: tuple[int, ...],
    quotient_representatives: tuple[int, ...],
    selected_lifts: dict[int, int],
    shortest_words: list[tuple[int, ...]],
) -> tuple[dict[int, int], int]:
    """Construct the partial-projection BLR section from the words C_q."""

    d = len(basis)
    k = len(quotient_representatives)
    full_coordinates = coordinate_lookup(basis + quotient_representatives)
    g_values = [
        xor_sum(selected_lifts[label] for label in shortest_words[q])
        for q in range(1 << k)
    ]
    coordinate_tables: list[int] = []
    for output_coordinate in range(d):
        table = 0
        for q in range(1, 1 << k):
            lift_0 = xor_sum(
                quotient_representatives[i]
                for i in range(k)
                if (q >> i) & 1
            )
            value = (full_coordinates[g_values[q] ^ lift_0] >> output_coordinate) & 1
            table |= value << (q - 1)
        coordinate_tables.append(table)

    alphas = [best_fourier_character(table, k)[0] for table in coordinate_tables]
    section = {0: 0}
    maximum_error = 0
    for q in range(1, 1 << k):
        lift_0 = xor_sum(
            quotient_representatives[i]
            for i in range(k)
            if (q >> i) & 1
        )
        correction_mask = sum(
            ((popcount(alpha & q) & 1) << j) for j, alpha in enumerate(alphas)
        )
        correction = xor_sum(
            basis[j] for j in range(d) if (correction_mask >> j) & 1
        )
        section[q] = lift_0 ^ correction
        maximum_error = max(
            maximum_error,
            popcount(full_coordinates[g_values[q] ^ section[q]]),
        )
    return section, maximum_error


def exhaustive_partial_projection_audit(max_w: int) -> dict[str, object]:
    """Check fibre stripping and partial synchronization in every binary case."""

    results: dict[str, object] = {}
    for w in range(1, max_w + 1):
        universe = tuple(range(1, 1 << w))
        selected_transversals = 0
        future_profiles = 0
        largest_h = 0
        largest_stripping_gap = 0
        largest_partial_graph_minus_source = 0
        largest_source_minus_partial_graph = 0
        largest_partial_section_error = 0
        worst_witness = None

        for support_mask in range(1 << len(universe)):
            support = tuple(
                universe[i]
                for i in range(len(universe))
                if (support_mask >> i) & 1
            )
            if rank_binary(support, w) != w:
                continue
            source_lengths = word_lengths(support, w)
            for d in range(1, w + 1):
                for basis in itertools.combinations(support, d):
                    if rank_binary(basis, w) != d:
                        continue
                    if source_lengths[xor_sum(basis)] != d:
                        continue
                    k = w - d
                    if k == 0:
                        # Cycle contraction forces S=B; this is the exact endpoint.
                        if set(support) != set(basis):
                            raise AssertionError("nontrivial zero-quotient geodesic support")
                        continue

                    quotient_representatives: list[int] = []
                    current_rank = d
                    for generator in support:
                        if generator in basis:
                            continue
                        new_rank = rank_binary(
                            list(basis) + quotient_representatives + [generator], w
                        )
                        if new_rank > current_rank:
                            quotient_representatives.append(generator)
                            current_rank = new_rank
                    if len(quotient_representatives) != k:
                        raise AssertionError("failed to choose quotient basis")
                    quotient_representatives_tuple = tuple(quotient_representatives)
                    full_coordinates = coordinate_lookup(
                        basis + quotient_representatives_tuple
                    )
                    fibres: dict[int, list[int]] = {}
                    for generator in support:
                        if generator in basis:
                            continue
                        label = full_coordinates[generator] >> d
                        fibres.setdefault(label, []).append(generator)
                    labels = tuple(sorted(fibres))
                    quotient_distances, shortest_words = quotient_shortest_words(
                        labels, k
                    )
                    h = max(quotient_distances)
                    largest_h = max(largest_h, h)

                    for chosen in itertools.product(*(fibres[label] for label in labels)):
                        selected_lifts = dict(zip(labels, chosen))
                        stripped = tuple(basis) + tuple(chosen)
                        section, section_error = partial_graph_section(
                            basis,
                            quotient_representatives_tuple,
                            selected_lifts,
                            shortest_words,
                        )
                        if section_error > 9 * h:
                            raise AssertionError("partial BLR section bound failed")
                        graph = tuple(basis) + tuple(
                            section[q] for q in range(1, 1 << k)
                        )
                        selected_transversals += 1
                        largest_partial_section_error = max(
                            largest_partial_section_error, section_error
                        )

                        for appended_mask in range(1 << len(universe)):
                            appended = tuple(
                                universe[i]
                                for i in range(len(universe))
                                if (appended_mask >> i) & 1
                            )
                            source_future = word_lengths(
                                tuple(set(support) | set(appended)), w
                            )
                            stripped_future = word_lengths(
                                tuple(set(stripped) | set(appended)), w
                            )
                            graph_future = word_lengths(
                                tuple(set(graph) | set(appended)), w
                            )
                            future_profiles += 1
                            stripping_gap = max(
                                b - a for a, b in zip(source_future, stripped_future)
                            )
                            stripping_reverse = max(
                                a - b for a, b in zip(source_future, stripped_future)
                            )
                            graph_minus_source = max(
                                b - a for a, b in zip(source_future, graph_future)
                            )
                            source_minus_graph = max(
                                a - b for a, b in zip(source_future, graph_future)
                            )
                            if stripping_reverse != 0 or stripping_gap > 2 * h:
                                raise AssertionError("binary fibre-stripping bound failed")
                            if (
                                graph_minus_source > 10 * h + 1
                                or source_minus_graph > 10 * h - 1
                            ):
                                raise AssertionError(
                                    "partial-projection response bound failed"
                                )
                            if (
                                stripping_gap > largest_stripping_gap
                                or graph_minus_source > largest_partial_graph_minus_source
                                or source_minus_graph > largest_source_minus_partial_graph
                            ):
                                worst_witness = {
                                    "support": support,
                                    "basis": basis,
                                    "h": h,
                                    "selected_lifts": selected_lifts,
                                    "appended": appended,
                                    "stripping_gap": stripping_gap,
                                    "graph_minus_source": graph_minus_source,
                                    "source_minus_graph": source_minus_graph,
                                }
                            largest_stripping_gap = max(
                                largest_stripping_gap, stripping_gap
                            )
                            largest_partial_graph_minus_source = max(
                                largest_partial_graph_minus_source,
                                graph_minus_source,
                            )
                            largest_source_minus_partial_graph = max(
                                largest_source_minus_partial_graph,
                                source_minus_graph,
                            )

        results[str(w)] = {
            "selected_transversals_checked": selected_transversals,
            "future_profiles_checked": future_profiles,
            "largest_quotient_diameter": largest_h,
            "largest_fibre_stripping_gap": largest_stripping_gap,
            "largest_partial_section_error": largest_partial_section_error,
            "largest_partial_graph_minus_source": (
                largest_partial_graph_minus_source
            ),
            "largest_source_minus_partial_graph": (
                largest_source_minus_partial_graph
            ),
            "worst_witness": worst_witness,
        }
    return results


def sharp_fibre_stripping_audit(max_h: int = 5) -> dict[str, object]:
    """Check the factor-two construction and its selector response cube."""

    results: dict[str, object] = {}
    for h in range(1, max_h + 1):
        kernel_dimension = 2 * h
        ambient_dimension = 3 * h
        basis = tuple(1 << i for i in range(kernel_dimension))
        plain = tuple(1 << (kernel_dimension + i) for i in range(h))
        perturbed = tuple(
            plain[i] ^ (1 << (2 * i)) ^ (1 << (2 * i + 1))
            for i in range(h)
        )
        source = basis + plain + perturbed
        stripped = basis + plain

        cycle_words = 0
        lifts = plain + perturbed
        quotient_mask = ((1 << h) - 1) << kernel_dimension
        kernel_mask = (1 << kernel_dimension) - 1
        for mask in range(1 << len(lifts)):
            selected = tuple(
                lifts[i] for i in range(len(lifts)) if (mask >> i) & 1
            )
            total = xor_sum(selected)
            if total & quotient_mask:
                continue
            cycle_words += 1
            if popcount(total & kernel_mask) > len(selected):
                raise AssertionError("sharp construction violates cycle contraction")

        target = xor_sum(perturbed)
        source_lengths = word_lengths(source, ambient_dimension)
        stripped_lengths = word_lengths(stripped, ambient_dimension)
        if source_lengths[target] != h or stripped_lengths[target] != 3 * h:
            raise AssertionError("factor-two sharp target has wrong length")

        selector_checks = 0
        largest_selector_distance = 0
        selector_profiles: dict[int, list[int]] = {}
        for j_mask in range(1 << h):
            support_j = basis + plain + tuple(
                perturbed[i] for i in range(h) if (j_mask >> i) & 1
            )
            lengths_j = word_lengths(support_j, ambient_dimension)
            profile: list[int] = []
            for p_mask in range(1 << h):
                x_p = xor_sum(
                    perturbed[i] for i in range(h) if (p_mask >> i) & 1
                )
                predicted = popcount(p_mask) + 2 * popcount(
                    p_mask & (~j_mask & ((1 << h) - 1))
                )
                if lengths_j[x_p] != predicted:
                    raise AssertionError("selector root formula failed")
                profile.append(lengths_j[x_p])
                selector_checks += 1
            selector_profiles[j_mask] = profile
        for j_mask in range(1 << h):
            for k_mask in range(1 << h):
                actual = max(
                    abs(a - b)
                    for a, b in zip(
                        selector_profiles[j_mask], selector_profiles[k_mask]
                    )
                )
                predicted = 2 * max(
                    popcount(j_mask & ~k_mask), popcount(k_mask & ~j_mask)
                )
                if actual != predicted:
                    raise AssertionError("selector response metric formula failed")
                largest_selector_distance = max(largest_selector_distance, actual)

        results[str(h)] = {
            "ambient_dimension": ambient_dimension,
            "projected_zero_subsets_checked": cycle_words,
            "source_target_length": source_lengths[target],
            "stripped_target_length": stripped_lengths[target],
            "sharp_gap": stripped_lengths[target] - source_lengths[target],
            "selector_root_values_checked": selector_checks,
            "largest_selector_response_distance": largest_selector_distance,
        }
    return results


def exhaustive_response_audit(max_w: int) -> dict[str, object]:
    results: dict[str, object] = {}
    for w in range(1, max_w + 1):
        universe = tuple(range(1, 1 << w))
        cases = 0
        appended_profiles = 0
        maximum_graph_minus_source = 0
        maximum_source_minus_graph = 0
        maximum_selected_average = 0.0
        maximum_selected_uniform = 0
        worst_witness = None

        for support_mask in range(1 << len(universe)):
            support = tuple(
                universe[i]
                for i in range(len(universe))
                if (support_mask >> i) & 1
            )
            if rank_binary(support, w) != w:
                continue
            source_lengths = word_lengths(support, w)
            for d in range(1, w + 1):
                for basis in itertools.combinations(support, d):
                    if rank_binary(basis, w) != d:
                        continue
                    target = xor_sum(basis)
                    if source_lengths[target] != d:
                        continue

                    quotient_representatives: list[int] = []
                    current_rank = d
                    for generator in support:
                        if generator in basis:
                            continue
                        new_rank = rank_binary(
                            list(basis) + quotient_representatives + [generator], w
                        )
                        if new_rank > current_rank:
                            quotient_representatives.append(generator)
                            current_rank = new_rank
                    k = w - d
                    if len(quotient_representatives) != k:
                        raise AssertionError("failed to choose a quotient basis")

                    full_coordinates = coordinate_lookup(
                        basis + tuple(quotient_representatives)
                    )
                    fibres = {q: [] for q in range(1, 1 << k)}
                    for generator in support:
                        if generator in basis:
                            continue
                        q = full_coordinates[generator] >> d
                        if q:
                            fibres[q].append(generator)
                    if any(not fibre for fibre in fibres.values()):
                        continue

                    fibre_order = tuple(range(1, 1 << k))
                    for selected in itertools.product(
                        *(fibres[q] for q in fibre_order)
                    ):
                        transversal = dict(zip(fibre_order, selected))
                        section, average_numerator, uniform_error = blr_graph_section(
                            basis,
                            tuple(quotient_representatives),
                            transversal,
                        )
                        quotient_size = 1 << k
                        if average_numerator > 3 * quotient_size or uniform_error > 9:
                            raise AssertionError("BLR section violated its certified bounds")

                        graph = tuple(basis) + tuple(
                            section[q] for q in range(1, 1 << k)
                        )
                        cases += 1
                        maximum_selected_average = max(
                            maximum_selected_average,
                            average_numerator / quotient_size,
                        )
                        maximum_selected_uniform = max(
                            maximum_selected_uniform, uniform_error
                        )

                        for appended_mask in range(1 << len(universe)):
                            appended = tuple(
                                universe[i]
                                for i in range(len(universe))
                                if (appended_mask >> i) & 1
                            )
                            source_future = word_lengths(
                                tuple(set(support) | set(appended)), w
                            )
                            graph_future = word_lengths(
                                tuple(set(graph) | set(appended)), w
                            )
                            appended_profiles += 1
                            graph_minus_source = max(
                                b - a for a, b in zip(source_future, graph_future)
                            )
                            source_minus_graph = max(
                                a - b for a, b in zip(source_future, graph_future)
                            )
                            if graph_minus_source > 11 or source_minus_graph > 9:
                                raise AssertionError(
                                    "future-response constants failed: "
                                    f"w={w}, S={support}, B={basis}, U={appended}"
                                )
                            if (
                                graph_minus_source > maximum_graph_minus_source
                                or source_minus_graph > maximum_source_minus_graph
                            ):
                                worst_witness = {
                                    "support": support,
                                    "basis": basis,
                                    "transversal": transversal,
                                    "graph": graph,
                                    "appended": appended,
                                    "graph_minus_source": graph_minus_source,
                                    "source_minus_graph": source_minus_graph,
                                }
                            maximum_graph_minus_source = max(
                                maximum_graph_minus_source, graph_minus_source
                            )
                            maximum_source_minus_graph = max(
                                maximum_source_minus_graph, source_minus_graph
                            )

        results[str(w)] = {
            "dense_geodesic_transversal_cases": cases,
            "appended_support_profiles_checked": appended_profiles,
            "maximum_BLR_selected_average_error": maximum_selected_average,
            "maximum_BLR_selected_uniform_error": maximum_selected_uniform,
            "maximum_graph_minus_source_profile": maximum_graph_minus_source,
            "maximum_source_minus_graph_profile": maximum_source_minus_graph,
            "worst_witness": worst_witness,
        }
    return results


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-scalar-k", type=int, default=4)
    parser.add_argument("--max-geodesic-w", type=int, default=3)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name(
            "phase3_vector_blr_response_audit_results.json"
        ),
    )
    args = parser.parse_args()

    vector_cases = [
        vector_map_audit(2, output_dimension) for output_dimension in range(1, 5)
    ]
    vector_cases.extend(
        vector_map_audit(3, output_dimension) for output_dimension in range(1, 4)
    )
    result = {
        "claim": (
            "Pointwise vector-Hamming additive defect delta admits one linear "
            "map with average error at most delta and uniform error at most "
            "3 delta; a dense geodesic transversal therefore has an all-future "
            "linear-graph word-profile approximation with directional errors "
            "11 and 9."
        ),
        "scalar_blr": scalar_blr_audit(args.max_scalar_k),
        "joint_average_sharpness": joint_average_sharpness(),
        "vector_blr_modulo_linear_coordinates": vector_cases,
        "all_cycle_contractive_maps_on_F2_3": cycle_contracting_k3_audit(),
        "dense_geodesic_future_responses": exhaustive_response_audit(
            args.max_geodesic_w
        ),
        "partial_projection_and_fibre_stripping": (
            exhaustive_partial_projection_audit(args.max_geodesic_w)
        ),
        "sharp_fibre_stripping_selector_family": sharp_fibre_stripping_audit(),
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
