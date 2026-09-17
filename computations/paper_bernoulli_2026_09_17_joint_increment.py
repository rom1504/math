"""Exact finite replay of exchangeable full-sign increments.

All graph probabilities, conditional increments, widths, and VC checks
are rational/integer. Information and final logarithmic inequalities
are evaluated in double precision; these are diagnostics, not proofs.
Run from the repository with .venv/bin/python.
"""

from collections import defaultdict
from fractions import Fraction
from itertools import combinations, product
import json
import math


def words(n):
    return list(product((-1, 1), repeat=n))


def entropy(probabilities):
    return -sum(float(p) * math.log(float(p)) for p in probabilities if p)


def graph_law(base_n, base_signs, added):
    size = base_n + added
    spin_words = words(size)
    all_edges = list(combinations(range(size), 2))
    fixed_edges = list(combinations(range(base_n), 2))
    random_edges = [edge for edge in all_edges if edge[1] >= base_n]
    fixed = dict(zip(fixed_edges, base_signs))
    for random_signs in words(len(random_edges)):
        edges = dict(fixed)
        edges.update(zip(random_edges, random_signs))
        energies = [sum(a * x[i] * x[j] for (i, j), a in edges.items())
                    for x in spin_words]
        cap = max(map(abs, energies))
        yield edges, energies, cap, spin_words


def information_audit(base_n, base_signs, added):
    """Uniform absolute optimizer, including genuinely randomized ties."""
    size = base_n + added
    random_edge_count = base_n * added + added * (added - 1) // 2
    sample_count = 2 ** random_edge_count
    pz = defaultdict(Fraction)
    joints = [defaultdict(Fraction) for _ in range(added)]
    row_marginals = [defaultdict(Fraction) for _ in range(added)]
    mean_cap = Fraction(0)
    mean_row_response = [Fraction(0) for _ in range(added)]
    for edges, energies, cap, spin_words in graph_law(base_n, base_signs, added):
        mean_cap += Fraction(cap, sample_count)
        ground_indices = [i for i, val in enumerate(energies) if abs(val) == cap]
        atom_mass = Fraction(1, sample_count * len(ground_indices))
        for index in ground_indices:
            z = spin_words[index]
            pz[z] += atom_mass
            for label, vertex in enumerate(range(base_n, size)):
                remaining = [i for i in range(size) if i != vertex]
                row = tuple(edges[tuple(sorted((vertex, i)))] for i in remaining)
                joints[label][z, row] += atom_mass
                row_marginals[label][row] += atom_mass
                response = abs(sum(a * z[i] for a, i in zip(row, remaining)))
                mean_row_response[label] += atom_mass * response
    hz = entropy(pz.values())
    mis = []
    for joint, marginal in zip(joints, row_marginals):
        assert len(marginal) == 2 ** (size - 1)
        assert set(marginal.values()) == {Fraction(1, 2 ** (size - 1))}
        mi = sum(float(p) * math.log(float(p / (pz[z] * marginal[row])))
                 for (z, row), p in joint.items())
        mis.append(mi)
    assert max(mis) - min(mis) < 2e-12
    assert sum(mis) <= 2 * hz + 2e-12
    assert len(set(mean_row_response)) == 1
    return mean_cap, {
        "parent_order": size,
        "random_graphs": sample_count,
        "optimizer_entropy": hz,
        "row_mutual_information": mis,
        "information_sum_over_entropy": sum(mis) / hz if hz else 0,
        "mean_cap": str(mean_cap),
        "mean_row_response": str(mean_row_response[0]),
    }


def conditional_windows_audit(base_n, base_signs, added):
    size = base_n + added
    random_edge_count = base_n * added + added * (added - 1) // 2
    sample_count = 2 ** random_edge_count
    mean_cap, mean_delta = Fraction(0), Fraction(0)
    checked_windows = 0
    minimum_slack = None
    for _, energies, cap, spin_words in graph_law(base_n, base_signs, added):
        mean_cap += Fraction(cap, sample_count)
        dot_products = [[sum(a * b for a, b in zip(h, x)) for x in spin_words]
                        for h in words(size)]
        added_caps = [max(abs(e) + abs(dot) for e, dot in zip(energies, dots))
                      for dots in dot_products]
        delta = Fraction(sum(added_caps), 2 ** size) - cap
        assert delta >= 0
        mean_delta += delta / sample_count
        # Near-codes only change at these exact integer deficits.
        for threshold in sorted(set(cap - abs(e) for e in energies)):
            indices = [i for i, e in enumerate(energies) if cap - abs(e) <= threshold]
            width = Fraction(sum(max(dots[i] for i in indices)
                                 for dots in dot_products), 2 ** size)
            slack = threshold + delta - width
            assert slack >= 0
            minimum_slack = slack if minimum_slack is None else min(minimum_slack, slack)
            checked_windows += 1
    return mean_cap, mean_delta, checked_windows, minimum_slack


def vc_dimension(code, n):
    for size in range(n, -1, -1):
        for subset in combinations(range(n), size):
            patterns = {tuple(x[i] for i in subset) for x in code}
            if len(patterns) == 2 ** size:
                return size
    raise AssertionError("Nonempty code should shatter the empty set")


def code_audit():
    cases = 0
    max_ratio = 0.0
    for n in range(1, 4):
        ambient = words(n)
        for mask in range(1, 2 ** len(ambient)):
            code = [x for i, x in enumerate(ambient) if mask & (1 << i)]
            width = Fraction(sum(max(sum(a * b for a, b in zip(h, x)) for x in code)
                                 for h in ambient), 2 ** n)
            dimension = vc_dimension(code, n)
            assert dimension <= width
            sauer = sum(math.comb(n, i) for i in range(dimension + 1))
            assert len(code) <= sauer
            rhs = float(width) * math.log(math.e * n / float(width)) if width else 0.0
            assert math.log(len(code)) <= rhs + 1e-13
            if rhs:
                max_ratio = max(max_ratio, math.log(len(code)) / rhs)
            cases += 1
    return {"all_boolean_codes_orders_1_to_3": cases,
            "max_entropy_over_width_bound": max_ratio}


def main():
    specifications = [(0, (), 2), (1, (), 2), (2, (1,), 1),
                      (2, (-1,), 2), (3, (1, 1, 1), 1),
                      (3, (1, 1, -1), 2)]
    rows = []
    for base_n, signs, added in specifications:
        parent_mean, info = information_audit(base_n, signs, added + 1)
        old_mean, delta, count, slack = conditional_windows_audit(base_n, signs, added)
        assert parent_mean - old_mean == delta
        size = base_n + added
        mu = Fraction(sum(abs(sum(h)) for h in words(size)), 2 ** size)
        bound = float(mu) + 2 * math.sqrt(size * (size + 1) * math.log(2) / (added + 1))
        assert float(delta) <= bound + 1e-12
        rows.append({"base_order": base_n, "added_vertices": added,
                     "base_signs": signs, "mean_conditional_increment": str(delta),
                     "increment_bound": bound, "windows_checked": count,
                     "minimum_exact_width_slack": str(slack), "information": info})
    print(json.dumps({"status": "PASS", "graph_cases": rows,
                      "code_checks": code_audit()}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
