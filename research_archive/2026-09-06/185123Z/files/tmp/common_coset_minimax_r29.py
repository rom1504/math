#!/usr/bin/env python3
"""Finite audit of the row-good block-coset minimax game from (10.838).

This is evidence only.  It enumerates all translates of every block subgroup
of a prescribed shape, discards cosets whose *whole* spin support violates a
row-square cap, forms their exact selector hit sets (up to the one irrational
p^(3/2) factor), and solves the dual finite zero-sum game.
"""

from __future__ import annotations

import argparse
import math
from fractions import Fraction

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix

from envelope_block_cover_r27 import (
    A6,
    A8,
    A9,
    all_data,
    projective_spins,
    set_partitions_with_sizes,
)


MATRICES = {"A6": A6, "A8": A8, "A9": A9}


def normalize_projective(x: tuple[int, ...]) -> tuple[int, ...]:
    if x[0] == 1:
        return x
    return tuple(-v for v in x)


def subgroup(partition: tuple[tuple[int, ...], ...], n: int):
    """Projective block-constant subgroup as normalized sign tuples."""
    words = set()
    k = len(partition)
    for mask in range(1 << k):
        x = [1] * n
        for j, block in enumerate(partition):
            sign = -1 if (mask >> j) & 1 else 1
            for i in block:
                x[i] = sign
        words.add(normalize_projective(tuple(x)))
    assert len(words) == 1 << (k - 1)
    return tuple(sorted(words))


def all_distinct_cosets(partition, spins, spin_to_id):
    group = subgroup(partition, len(spins[0]))
    seen = set()
    for base in spins:
        ids = frozenset(
            spin_to_id[normalize_projective(tuple(a * b for a, b in zip(base, h)))]
            for h in group
        )
        if ids not in seen:
            seen.add(ids)
            yield ids


def solve(name: str, sizes: tuple[int, ...], m: int, cap: int, tolerance: float):
    A = MATRICES[name]
    n = len(A)
    spins_np = list(projective_spins(n))
    spins = [tuple(int(v) for v in x) for x in spins_np]
    spin_to_id = {x: i for i, x in enumerate(spins)}

    q, selectors, records, qS, den = all_data(A, m)
    # records contains both orientations.  Spin quantities below use one
    # projective representative and an absolute value, which is equivalent.
    energies = np.array([int(x @ A @ x) for x in spins_np], dtype=np.int64)
    row = np.array([int((A @ x) @ (A @ x)) for x in spins_np], dtype=np.int64)
    child = np.array([
        [int(x[list(S)] @ A[np.ix_(S, S)] @ x[list(S)]) for S in selectors]
        for x in spins_np
    ], dtype=np.int64)
    p2 = Fraction(m * (m - 1), n * (n - 1))
    Xnum = np.abs(
        p2.denominator * child - p2.numerator * energies[:, None]
    ).astype(np.int64)
    X = Xnum.astype(float) / p2.denominator
    Y = qS.astype(float) - q * (m / n) ** 1.5

    candidates = set()
    partitions = list(set_partitions_with_sizes(n, sizes))
    for partition in partitions:
        for ids in all_distinct_cosets(partition, spins, spin_to_id):
            if int(row[list(ids)].max()) <= cap:
                candidates.add(ids)

    hit_rows = []
    cover_counts = np.zeros(len(selectors), dtype=np.int64)
    best_cover = 0
    full = 0
    min_margin = math.inf
    for ids in candidates:
        xmax_num = Xnum[list(ids), :].max(axis=0)
        margin = xmax_num.astype(float) / p2.denominator - (Y - tolerance)
        if tolerance == 0.0:
            # Compare xmax >= qS-q(m/n)^(3/2) without floating point.
            # When xmax-qS<0 both sides of the rearranged inequality are
            # nonnegative, so squaring is equivalent.
            hit = np.empty(len(selectors), dtype=bool)
            alpha2 = Fraction(q * q * m**3, n**3)
            for j, num in enumerate(xmax_num):
                diff = Fraction(int(num), p2.denominator) - int(qS[j])
                hit[j] = diff >= 0 or diff * diff <= alpha2
        else:
            hit = margin >= -1e-10
        cols = np.flatnonzero(hit)
        hit_rows.append(cols)
        cover_counts[cols] += 1
        best_cover = max(best_cover, len(cols))
        full += int(len(cols) == len(selectors))
        min_margin = min(min_margin, float(np.min(np.abs(margin))))

    # Dual game: min_{w in Delta(selectors)} max_a sum_S w_S 1{a hits S}.
    # Variables are (w_1,...,w_N,v), constraints H_a w-v <= 0.
    nr, nc = len(hit_rows), len(selectors)
    hrows, hcols, hvals = [], [], []
    rows, cols, vals = [], [], []
    for r, js in enumerate(hit_rows):
        hrows.extend([r] * len(js))
        hcols.extend(int(j) for j in js)
        hvals.extend([1.0] * len(js))
        rows.extend([r] * len(js))
        cols.extend(int(j) for j in js)
        vals.extend([1.0] * len(js))
        rows.append(r)
        cols.append(nc)
        vals.append(-1.0)
    Aub = csr_matrix((vals, (rows, cols)), shape=(nr, nc + 1))
    H = csr_matrix((hvals, (hrows, hcols)), shape=(nr, nc))
    c = np.zeros(nc + 1)
    c[-1] = 1.0
    Aeq = np.zeros((1, nc + 1))
    Aeq[0, :nc] = 1.0
    result = linprog(
        c,
        A_ub=Aub,
        b_ub=np.zeros(nr),
        A_eq=Aeq,
        b_eq=np.array([1.0]),
        bounds=[(0.0, None)] * nc + [(0.0, 1.0)],
        method="highs",
    )
    assert result.success, result.message
    w = result.x[:nc]
    support = [(selectors[i], float(a)) for i, a in enumerate(w) if a > 1e-8]

    # Primal game: max delta subject to H^T mu >= delta and mu in Delta(A).
    pAub = csr_matrix(np.hstack((-H.T.toarray(), np.ones((nc, 1)))))
    pc = np.zeros(nr + 1)
    pc[-1] = -1.0
    pAeq = np.zeros((1, nr + 1))
    pAeq[0, :nr] = 1.0
    primal = linprog(
        pc,
        A_ub=pAub,
        b_ub=np.zeros(nc),
        A_eq=pAeq,
        b_eq=np.array([1.0]),
        bounds=[(0.0, None)] * nr + [(0.0, 1.0)],
        method="highs",
    )
    assert primal.success, primal.message
    assert abs(result.fun + primal.fun) < 1e-7

    # Try to turn both floating LP optima into exact rational certificates.
    # H is a zero-one matrix, so the following verification is exact once
    # componentwise rational reconstruction succeeds.
    def rationalize(v):
        return [Fraction(float(a)).limit_denominator(1_000_000) for a in v]

    wr = rationalize(w)
    mur = rationalize(primal.x[:nr])
    rational_value = None
    dual_pair_probability = None
    dual_triple_probability = None
    missing_pairs = None
    missing_triples = None
    if sum(wr) == 1 and sum(mur) == 1:
        dual_max = max(sum(wr[j] for j in js) for js in hit_rows)
        primal_min = min(
            sum(mur[a] for a in range(nr) if H[a, j])
            for j in range(nc)
        )
        if dual_max == primal_min:
            rational_value = str(dual_max)
            dual_pair_probability = sum(
                wr[i] * wr[j]
                for i in range(nc)
                for j in range(nc)
                if wr[i] and wr[j] and H[:, i].multiply(H[:, j]).nnz
            )
            col_bits = [0] * nc
            for a, js in enumerate(hit_rows):
                bit = 1 << a
                for j in js:
                    col_bits[int(j)] |= bit
            missing_pairs = sum(
                not (col_bits[i] & col_bits[j])
                for i in range(nc) for j in range(i, nc)
            )
            missing_triples = sum(
                not (col_bits[i] & col_bits[j] & col_bits[k])
                for i in range(nc)
                for j in range(i, nc)
                for k in range(j, nc)
            )
            dual_triple_probability = sum(
                wr[i] * wr[j] * wr[k]
                for i in range(nc)
                for j in range(nc)
                for k in range(nc)
                if wr[i] and wr[j] and wr[k]
                and (col_bits[i] & col_bits[j] & col_bits[k])
            )
    uniform_value = float(cover_counts.min() / max(1, len(candidates)))
    out = {
        "matrix": name,
        "n": n,
        "m": m,
        "shape": sizes,
        "cap": cap,
        "tolerance": tolerance,
        "partitions": len(partitions),
        "row_good_distinct_cosets": len(candidates),
        "selectors": len(selectors),
        "deterministic_full_cosets": full,
        "best_single_coverage": f"{best_cover}/{len(selectors)}",
        "uniform_coset_min_hit": uniform_value,
        "game_value": float(result.fun),
        "exact_rational_value": rational_value,
        "dual_pair_common_probability": (
            str(dual_pair_probability) if dual_pair_probability is not None else None
        ),
        "dual_triple_common_probability": (
            str(dual_triple_probability) if dual_triple_probability is not None else None
        ),
        "missing_unordered_pairs": missing_pairs,
        "missing_unordered_triples": missing_triples,
        "primal_support_size": int(np.sum(primal.x[:nr] > 1e-8)),
        "dual_support_size": len(support),
        "dual_support": support,
        "closest_hit_margin": min_margin,
    }
    print(out)
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("matrix", choices=sorted(MATRICES))
    parser.add_argument("m", type=int)
    parser.add_argument("cap", type=int)
    parser.add_argument("sizes", nargs="+", type=int)
    parser.add_argument("--tolerance", type=float, default=0.0)
    parser.add_argument("--expect-fraction")
    parser.add_argument("--expect-cosets", type=int)
    parser.add_argument("--expect-full", type=int)
    args = parser.parse_args()
    assert sum(args.sizes) == len(MATRICES[args.matrix])
    out = solve(args.matrix, tuple(args.sizes), args.m, args.cap, args.tolerance)
    if args.expect_fraction is not None:
        assert out["exact_rational_value"] == args.expect_fraction
    if args.expect_cosets is not None:
        assert out["row_good_distinct_cosets"] == args.expect_cosets
    if args.expect_full is not None:
        assert out["deterministic_full_cosets"] == args.expect_full
    if any(v is not None for v in (args.expect_fraction, args.expect_cosets, args.expect_full)):
        print("PASS: expected exact finite minimax data verified")


if __name__ == "__main__":
    main()
