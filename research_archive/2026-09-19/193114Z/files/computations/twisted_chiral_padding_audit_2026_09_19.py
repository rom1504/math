#!/usr/bin/env python3
"""Exact finite checks of the universal low-cap padding proof.

Exhausts all principal m-subsets and all ordered disjoint pairs for small
inputs. These checks support the algebra; the asymptotic proof is separate.
"""

from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def spins(n):
    ids = np.arange(1 << n, dtype=np.uint64)
    return 1 - 2 * ((ids[:, None] >> np.arange(n, dtype=np.uint64)) & 1).astype(np.int64)


def norms(a):
    xs = spins(len(a))
    ax = xs @ a
    q = int(np.max(np.abs(np.sum(ax * xs, axis=1)))) // 2
    values = np.sum(np.abs(ax), axis=1)
    j = int(np.argmax(values))
    y = xs[j]
    x = np.where(a @ y >= 0, 1, -1)
    beta = int(values[j])
    assert int(x @ a @ y) == beta
    return q, beta, x, y


def tr4(a):
    squared = a @ a
    return int(np.sum(squared * squared))


def falling(n, k):
    value = 1
    for j in range(k):
        value *= n-j
    return value


def repeated(n):
    return n*(n-1)*(2*n-3)


def audit(label, a, m):
    n = len(a)
    q, beta, x, y = norms(a)
    assert beta <= 4*q
    a4 = tr4(a)
    assert a4 <= n*(n-1)*beta
    sets = list(combinations(range(n), m))
    traces = {s: tr4(a[np.ix_(s, s)]) for s in sets}
    average = Fraction(sum(traces.values()), len(sets))
    predicted = Fraction(falling(m, 4), falling(n, 4))*(a4-repeated(n))+repeated(m)
    assert average == predicted
    bound_t = Fraction(m**4 * beta, n*n) + 2*m**3
    assert average <= bound_t
    expected_deficit = (1-Fraction(falling(n-2*m, 2), falling(n, 2)))*beta
    pairs = 0
    good = 0
    sum_deficit = 0
    selected = None
    for s1 in sets:
        leftover = tuple(i for i in range(n) if i not in s1)
        for s2 in combinations(leftover, m):
            comp = tuple(i for i in leftover if i not in s2)
            f = int(x[list(comp)] @ a[np.ix_(comp, comp)] @ y[list(comp)])
            deficit = beta-f
            assert deficit >= 0
            pairs += 1
            sum_deficit += deficit
            if traces[s1] <= 4*bound_t and traces[s2] <= 4*bound_t and deficit <= 4*expected_deficit:
                good += 1
                if selected is None:
                    selected = (s1, s2, comp, f)
    assert Fraction(sum_deficit, pairs) == expected_deficit
    assert 4*good >= pairs
    s1, s2, comp, f = selected
    patched = a.copy()
    for s, sign in [(s1, 1), (s2, -1)]:
        for i in s:
            for j in s:
                patched[i, j] = sign if i != j else 0
    qp, bp, _, _ = norms(patched)
    qs = [norms(a[np.ix_(s, s)])[0] for s in [s1, s2]]
    assert qp <= q + sum(qs) + m*m//2
    assert bp >= f + 2*m*(m-1)
    assert bp >= beta + 2*m*(m-1) - 4*expected_deficit
    for small_q in qs:
        assert 16*small_q**4 <= m**4 * 4*bound_t
    return dict(label=label, n=n, m=m, Q=q, beta=beta, tr4=a4,
                subsets=len(sets), exact_expected_tr4=str(average),
                trace_upper=str(bound_t), disjoint_pairs=pairs,
                simultaneous_good_pairs=good, exact_expected_deficit=str(expected_deficit),
                selected_blocks=[s1, s2], patched_Q=qp, patched_beta=bp,
                all_checks_pass=True)


def main():
    rng = np.random.default_rng(20260919)
    records = []
    for n in [6, 7, 8, 9, 10]:
        a = rng.choice([-1, 1], size=(n, n)).astype(np.int64)
        a = np.triu(a, 1)
        a += a.T
        records.append(audit(f"random{n}", a, min(4, n//2)))
    census = json.loads((ROOT/'computations/results/twisted_chiral_2026_09_18_classify10.json').read_text())
    for cls in census['classes']:
        records.append(audit(f"min10class{cls['class']}", np.array(cls['matrix'], dtype=np.int64), 4))
    output = dict(method='exact integer/fraction exhaustive small principal-block audit', records=records)
    target = ROOT/'computations/results/twisted_chiral_padding_audit_2026_09_19.json'
    target.write_text(json.dumps(output, indent=2)+'\n')
    print(json.dumps(output, indent=2))


if __name__ == '__main__':
    main()
