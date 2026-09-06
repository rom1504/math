"""Exact full-cube checks of the coherent twin-response counterexample.

The growing-order claim is proved separately. These are finite integer
identities, not a claim that these signings are near-minimizers.
"""
import argparse
from fractions import Fraction
import json
from pathlib import Path
import numpy as np


def hadamard(k):
    h = np.ones((1, 1), dtype=np.int64)
    while len(h) < k:
        h = np.block([[h, h], [h, -h]])
    return h


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    rows = []
    for k in (2, 4, 8):
        n = 2 * k
        h = hadamard(k)
        big = np.kron(np.ones((2, 2), dtype=np.int64), h)
        a = big - np.diag(np.diag(big))
        t = np.kron(np.ones((2, 2), dtype=np.int64), np.eye(k, dtype=np.int64))
        q_num = a @ a
        defect_num = q_num - (n - 1) * t
        states = 2 ** n
        ids = np.arange(states, dtype=np.uint64)
        spins = (2 * ((ids[:, None] >> np.arange(n, dtype=np.uint64)) & 1)
                 .astype(np.int64) - 1)
        field_num = spins @ q_num
        f = np.sign(field_num) * (np.abs(field_num) > n - 1)
        u = np.column_stack(((spins[:, :k] + spins[:, k:]) // 2,
                             (spins[:, :k] + spins[:, k:]) // 2))
        mse = Fraction(int(np.sum((f - u) ** 2)), n * states)
        rhs = Fraction(4 * int(np.sum(defect_num ** 2)), n * (n - 1) ** 2)
        assert mse <= rhs
        f2 = Fraction(int(np.sum(f * f)), n * states)
        coefficient_num = f.T @ spins
        degree_one_mass = Fraction(int(np.sum(coefficient_num ** 2)),
                                   n * states ** 2)
        assert np.all(f.sum(axis=0) == 0)
        higher_mass = f2 - degree_one_mass
        assert 0 <= higher_mass <= mse
        quadratic = np.sum((spins @ a) * spins, axis=1)
        assert np.all(quadratic % 2 == 0)
        rows.append({'k': k, 'n': n, 'full_cube_states': states,
                     'actual_cap': int(np.max(np.abs(quadratic))) // 2,
                     'marked_density': str(f2), 'hole_density': str(1 - f2),
                     'degree_ge_2_mass': str(higher_mass),
                     'distance_to_linear_twin_response': str(mse),
                     'exact_frobenius_error_bound': str(rhs)})
    report = {'status': 'PASS', 'evidence': 'exact integer and Fraction finite checks',
              'construction': 'Sylvester H; A=J2 tensor H-diag; F=sign(QS) 1{|QS|>1}',
              'scope': 'Bounded-operator signing family, not near-minimizers.', 'cases': rows}
    rendered = json.dumps(report, indent=2) + '\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
