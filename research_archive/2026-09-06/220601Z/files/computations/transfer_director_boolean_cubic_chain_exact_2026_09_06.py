"""Exact finite regressions of cubic conditional cancellation and L4 estimate.

No asymptotic or Gaussian comparison is inferred from these checks.
Python standard library only; seed and all test parameters are recorded.
"""
import argparse
from fractions import Fraction
import itertools
import json
import random
from pathlib import Path


def value(poly, spin):
    return sum(c * prod(spin[a] for a in support)
               for support, c in poly.items())


def prod(items):
    result = 1
    for item in items:
        result *= item
    return result


def derivative(poly, coordinate, spin):
    return sum(c * prod(spin[a] for a in support if a != coordinate)
               for support, c in poly.items() if coordinate in support)


def polynomial(rng, n, degree, pure=False):
    result = {}
    degrees = [degree] if pure else range(degree + 1)
    for k in degrees:
        for support in itertools.combinations(range(n), k):
            coefficient = rng.randint(-2, 2)
            if coefficient:
                result[support] = coefficient
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    seed = 202609062141
    rng = random.Random(seed)
    cases = []
    fourth_moment_checks = 0
    for n in range(2, 9):
        for degree in range(1, min(n, 3) + 1):
            for probe_degree in range(1, min(n, 3) + 1):
                v1 = polynomial(rng, n, degree)
                v2 = polynomial(rng, n, degree)
                z = polynomial(rng, n, probe_degree, pure=True)
                lhs_sum = bracket_sum = cubic_sum = 0
                fourth = [0, 0, 0]
                for spin in itertools.product((-1, 1), repeat=n):
                    a, b, c = (value(p, spin) for p in (v1, v2, z))
                    lhs_sum += c * a * a * b
                    for coordinate in range(n):
                        da, db, dc = (derivative(p, coordinate, spin)
                                      for p in (v1, v2, z))
                        bracket_sum += dc * (2 * a * b * da + a * a * db)
                        cubic_sum += dc * da * da * db
                    for j, val in enumerate((a, b, c)):
                        fourth[j] += val ** 4
                denominator = (2 ** n) * probe_degree
                defect = Fraction(probe_degree * lhs_sum - bracket_sum,
                                  denominator)
                exact_remainder = Fraction(-2 * cubic_sum, denominator)
                assert defect == exact_remainder
                for p, moment_sum in zip((v1, v2, z), fourth):
                    weighted = sum((3 ** len(s)) * c * c for s, c in p.items())
                    assert moment_sum <= (2 ** n) * weighted * weighted
                    fourth_moment_checks += 1
                cases.append({'cube_dimension': n, 'query_degree': degree,
                              'probe_degree': probe_degree,
                              'chain_defect': str(defect),
                              'predicted_cubic_remainder': str(exact_remainder)})
    report = {'status': 'PASS', 'evidence': 'exact integer/Fraction finite regression',
              'seed': seed, 'test_function': 'f(v1,v2)=v1^2*v2',
              'coefficient_range': [-2, 2], 'case_count': len(cases),
              'nonzero_defects': sum(c['chain_defect'] != '0' for c in cases),
              'fourth_moment_checks': fourth_moment_checks, 'cases': cases,
              'scope': 'No assertion about arbitrary rich signing feedback closure.'}
    rendered = json.dumps(report, indent=2) + '\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
