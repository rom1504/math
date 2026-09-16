"""Exact n=16 planted-clique tail check; no solver or optimality claim."""

from collections import Counter
from fractions import Fraction
import json

from bh_2026_09_16_power_checks import energies, sylvester_signing


def main():
    n, k = 16, 8
    n_three_halves = 64
    base = sylvester_signing(n)
    planted = base.copy()
    original_clique_sum = sum(int(base[i, j]) for i in range(k)
                              for j in range(i + 1, k))
    assert original_clique_sum == k // 2
    for i in range(k):
        for j in range(i + 1, k):
            planted[i, j] = planted[j, i] = 1
    delta_l1 = sum(int(planted[i, j] - base[i, j])
                   for i in range(n) for j in range(i + 1, n))
    assert delta_l1 == k * (k - 2) // 2
    values = energies(planted)
    cap = max(abs(v) for v in values)
    assert cap == 46
    cap_upper = n_three_halves - k
    assert cap <= cap_upper
    clique_mask = (1 << k) - 1
    conditioned = [v for mask, v in enumerate(values)
                   if mask & clique_mask == 0]
    mean = Fraction(sum(conditioned), len(conditioned))
    assert mean == Fraction(k * (k - 1), 2)
    threshold = n_three_halves // 4
    probability = Fraction(sum(v >= threshold for v in conditioned),
                           len(conditioned))
    assert probability == Fraction(115, 128)
    assert probability >= Fraction(mean - threshold, cap_upper - threshold)
    print(json.dumps({
        "status": "PASS", "n": n, "k": k,
        "actual_cap": cap, "proved_cap_upper": cap_upper,
        "clique_original_sign_sum": original_clique_sum,
        "modification_l1": delta_l1,
        "conditional_mean": str(mean),
        "conditional_histogram": dict(sorted(Counter(conditioned).items())),
        "threshold": threshold,
        "conditional_tail_probability": str(probability),
        "unconditional_tail_lower": str(probability / (1 << k)),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
