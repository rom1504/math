#!/usr/bin/env python3
"""Exact endpoint/beta profiles and a finite single-clique proof audit."""

from fractions import Fraction
from itertools import combinations
import hashlib
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def profile(a):
    n = len(a)
    ids = np.arange(1 << max(0, n-1), dtype=np.int64)
    x = np.column_stack((np.ones(len(ids), dtype=np.int64),
                         1-2*((ids[:, None] >> np.arange(n-1)) & 1)))
    ax = x @ a
    h = np.sum(ax*x, axis=1)//2
    up, down = int(h.argmax()), int(h.argmin())
    bv = np.sum(abs(ax), axis=1)
    j = int(bv.argmax())
    beta_y = x[j]
    beta_x = np.where(a @ beta_y >= 0, 1, -1)
    u, v, beta = int(h[up]), int(-h[down]), int(bv[j])
    assert beta <= 2*(u+v)
    return u, v, beta, x[up], x[down], beta_x, beta_y


def single_clique_audit(a, m):
    n = len(a)
    u, v, beta, xu, xv, xb, yb = profile(a)
    theta = 1-Fraction((n-m)*(n-m-1), n*(n-1))
    # The analytic bound contains sqrt(m); choose m=4 so this audit stays rational.
    assert m == 4
    beta_bound = Fraction(32*m*m*beta, n*n) + 64*8
    sets = list(combinations(range(n), m))
    deficits = [0, 0, 0]
    betas = 0
    good = []
    for s in sets:
        c = tuple(i for i in range(n) if i not in s)
        sub = a[np.ix_(c, c)]
        du = u-int(xu[list(c)] @ sub @ xu[list(c)])//2
        dv = v+int(xv[list(c)] @ sub @ xv[list(c)])//2
        db = beta-int(xb[list(c)] @ sub @ yb[list(c)])
        assert min(du, dv, db) >= 0
        small_beta = profile(a[np.ix_(s, s)])[2]
        betas += small_beta
        deficits = [old+inc for old, inc in zip(deficits, [du, dv, db])]
        if (small_beta <= 5*beta_bound and du <= 5*theta*u
                and dv <= 5*theta*v and db <= 5*theta*beta):
            good.append(s)
    assert [Fraction(t, len(sets)) for t in deficits] == [theta*t for t in [u, v, beta]]
    assert Fraction(betas, len(sets)) <= beta_bound
    assert 5*len(good) >= len(sets)
    s = good[0]
    patched = a.copy()
    for i in s:
        for j in s:
            patched[i, j] = -1 if i != j else 0
    up, vp, bp, *_ = profile(patched)
    assert u+m//2-5*theta*u <= up <= u+m//2+Fraction(5, 2)*beta_bound
    assert v+m*(m-1)//2-5*theta*v <= vp <= v+m*(m-1)//2+Fraction(5, 2)*beta_bound
    assert beta+m*(m-1)-5*theta*beta <= bp <= beta+m*(m-1)+5*beta_bound
    return dict(n=n, m=m, subsets=len(sets), simultaneous_good_subsets=len(good),
                exact_mean_deficits=[str(Fraction(t, len(sets))) for t in deficits],
                exact_mean_block_beta=str(Fraction(betas, len(sets))),
                beta_expectation_bound=str(beta_bound), selected_block=s,
                old_U=u, old_V=v, old_beta=beta, new_U=up, new_V=vp, new_beta=bp,
                finite_inequalities_26_verified=True)


def main():
    inputs = []
    old = json.loads((ROOT/'computations/results/twisted_chiral_bilinear_audit_2026_09_18.json').read_text())
    inputs.extend((row['matrix'], row['sources']) for row in old['profiles'])
    for order in [9, 10]:
        source = f'computations/results/twisted_chiral_2026_09_18_classify{order}.json'
        census = json.loads((ROOT/source).read_text())
        inputs.extend((row['matrix'], [dict(path=source, class_index=row['class'])]) for row in census['classes'])
    records, audits, seen = [], [], set()
    for matrix, sources in inputs:
        a = np.asarray(matrix, dtype=np.int64)
        digest = hashlib.sha256(a.astype(np.int8).tobytes()).hexdigest()
        if digest in seen:
            continue
        seen.add(digest)
        u, v, beta, *_ = profile(a)
        q = max(u, v)
        joint = beta+2*abs(u-v)
        record = dict(n=len(a), Q=q, U=u, V=v, beta=beta,
                      balanced_padding_joint_invariant=joint,
                      joint_exceeds_2sqrt2_Q=joint*joint > 8*q*q,
                      matrix_sha256_int8=digest, sources=sources)
        records.append(record)
        if len(a) == 10:
            audits.append(single_clique_audit(a, 4))
    output = dict(profiles=records, single_clique_audits=audits,
                  scope='Finite exact profiles only; no asymptotic seed family or exact-minimizer preservation inferred.')
    target = ROOT/'computations/results/twisted_chiral_endpoint_balance_2026_09_19.json'
    target.write_text(json.dumps(output, indent=2)+'\n')
    for row in records:
        print(json.dumps({key: row[key] for key in ['n', 'Q', 'U', 'V', 'beta', 'balanced_padding_joint_invariant', 'joint_exceeds_2sqrt2_Q']}))
    print('single-clique audits:', len(audits), 'PASS')


if __name__ == '__main__':
    main()
