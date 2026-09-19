#!/usr/bin/env python3
"""Independent scalar replay of the order3/4 stable first-moment averages.

No conditional permanents or NumPy field updates are used in the scalar
replay.  Equal bridge matrices are aggregated with their exact group weights.
The stored direct witnesses also disprove a projective tail cutoff above1.
"""
from collections import Counter
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def parent_matrix(a, p, s, d):
    n = len(a)
    b = [[s[i] * s[j] * a[p[i]][p[j]] for j in range(n)] for i in range(n)]
    c = [[b[i][j] + (d[i] if i == j else 0) for j in range(n)] for i in range(n)]
    return [a[i] + c[i] for i in range(n)] + [c[i] + [-v for v in a[i]] for i in range(n)]


def histograms(matrix):
    size = len(matrix)
    raw, stable = Counter(), Counter()
    stable_words = []
    for tail in itertools.product([-1, 1], repeat=size - 1):
        x = (1,) + tail
        signed_fields = [x[i] * sum(matrix[i][j] * x[j] for j in range(size)) for i in range(size)]
        assert all(v % 2 for v in signed_fields)
        energy = sum(signed_fields) // 2
        raw[energy] += 1
        if min(signed_fields) > 0:
            stable[energy] += 1
            stable_words.append({"spin": x, "energy": energy, "signed_fields": signed_fields})
    return raw, stable, stable_words


def fraction_data(numerator, denominator):
    value = Fraction(numerator, denominator)
    return {"numerator": value.numerator, "denominator": value.denominator}


def main():
    source = ROOT / "computations/results/twisted_chiral_stable_average_2026_09_19.json"
    saved = json.loads(source.read_text())
    results = []
    for rec in saved["records"]:
        n, a = rec["n"], rec["matrix"]
        signs = list(itertools.product([-1, 1], repeat=n))
        orbit = Counter()
        for p in itertools.permutations(range(n)):
            for s in signs:
                b = tuple(tuple(s[i] * s[j] * a[p[i]][p[j]] for j in range(n)) for i in range(n))
                orbit[b] += 1
        assert sum(orbit.values()) == math.factorial(n) * 2**n
        raw, stable = Counter(), Counter()
        for b, weight in orbit.items():
            for d in signs:
                c = [[b[i][j] + (d[i] if i == j else 0) for j in range(n)] for i in range(n)]
                matrix = [a[i] + c[i] for i in range(n)] + [c[i] + [-v for v in a[i]] for i in range(n)]
                one_raw, one_stable, _ = histograms(matrix)
                raw.update({e: weight * count for e, count in one_raw.items()})
                stable.update({e: weight * count for e, count in one_stable.items()})
        assert dict(sorted(raw.items())) == {int(e): v for e, v in rec["positive_raw_energy_histogram_numerators"].items()}
        assert dict(sorted(stable.items())) == {int(e): v for e, v in rec["positive_stable_energy_histogram_numerators"].items()}
        denominator = 4**n * math.factorial(n)
        assert sum(raw.values()) == denominator * 2**(2*n - 1)
        scale = rec["factorial_root_upper_scale"]
        for degree, ceiling in enumerate(rec["factorial_root_upper_numerators"]):
            if degree == 0:
                assert ceiling == 0
            else:
                target = math.factorial(degree) * scale**degree
                assert (ceiling - 1)**degree < target <= ceiling**degree
        tails = []
        for row in rec["thresholds"]:
            threshold = row["L"]
            value = fraction_data(sum(count for e, count in stable.items() if e > threshold), denominator)
            assert value == {k: row["exact_positive_stable_above"][k] for k in value}
            tails.append({"L": threshold, "stable_tail": value})
        first = next(row["L"] for row in tails if row["stable_tail"]["numerator"] < row["stable_tail"]["denominator"])
        assert first == {3: 7, 4: 12}[n]

        # For n3 use the analytic clique construction.  For the actual n4
        # minimizing child use a tail-specific unique-violation counterexample.
        p = list(range(n))
        s = [1] * n if n == 3 else [1, -1, -1, -1]
        d = [1] * n if n == 3 else [-1] * n
        matrix = parent_matrix(a, p, s, d)
        one_raw, one_stable, words = histograms(matrix)
        threshold = 5 if n == 3 else 10
        assert sum(count for e, count in one_stable.items() if e > threshold) == 1
        assert sorted(one_stable.elements()) == ({3: [9], 4: [10, 10, 14]}[n])
        results.append({"n": n, "scalar_full_group_replay_pass": True,
                        "group_size": denominator, "distinct_bridges": len(orbit),
                        "first_exact_stable_certificate_L": first,
                        "stable_tail_at_family_minimum": next(row for row in tails if row["L"] == threshold),
                        "cutoff1_witness": {"p": p, "s": s, "d": d, "L": threshold,
                            "parent_matrix": matrix, "stable_states": words,
                            "stable_histogram": dict(sorted(one_stable.items()))}})
    output = {"scope": "Complete scalar replay for the two stated children and exact counterexamples to any larger universal projective tail multiplicity.",
              "all_pass": True, "records": results}
    path = ROOT / "computations/results/twisted_chiral_stable_average_adversary_2026_09_19.json"
    path.write_text(json.dumps(output, indent=2) + "\n")
    for rec in results:
        print(json.dumps({k: v for k, v in rec.items() if k != "cutoff1_witness"}))


if __name__ == "__main__":
    main()
