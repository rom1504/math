#!/usr/bin/env python3
"""Exact finite audit of cylinder certificates for the bare tail (10.795).

The chart consists of one orientation bit and n-1 projective spin bits
(x_0=1).  For every selector, we find a minimum-codimension coordinate
subcube all of whose oriented cuts satisfy widehat ell <= t.  All energies
and row squares are exact integers; only p^(3/2) is evaluated in floating
point, and the minimum distance from a decision boundary is reported.
"""

from __future__ import annotations

import itertools
import math
from collections import Counter
from fractions import Fraction

import numpy as np

from envelope_block_cover_r27 import A6, A8, A9, projective_spins


def qnorm(A: np.ndarray) -> int:
    return max(abs(int(x @ A @ x)) for x in projective_spins(len(A)))


def leq_p32_q(lhs: Fraction, q: int, m: int, n: int) -> bool:
    """Decide lhs <= q*(m/n)^(3/2) using integer arithmetic."""
    if lhs <= 0:
        return True
    return (lhs.numerator ** 2) * (n ** 3) <= (
        (q ** 2) * (m ** 3) * (lhs.denominator ** 2)
    )


def chart_words(n: int):
    """Return (bits, sigma, x) in lexicographic Boolean-chart order."""
    out = []
    for obit in (0, 1):
        sigma = -1 if obit == 0 else 1
        for tailbits in itertools.product((0, 1), repeat=n - 1):
            x = np.array((1,) + tuple(-1 if b == 0 else 1 for b in tailbits),
                         dtype=np.int64)
            out.append(((obit,) + tailbits, sigma, x))
    return out


def all_cylinders(words, rows: np.ndarray):
    """Precompute coordinate cylinders as integer bit masks by codimension."""
    bits = np.asarray([w[0] for w in words], dtype=np.int8)
    nbits = bits.shape[1]
    groups = []
    for k in range(nbits + 1):
        group = []
        for coords in itertools.combinations(range(nbits), k):
            for vals in itertools.product((0, 1), repeat=k):
                keep = np.ones(len(words), dtype=bool)
                for j, v in zip(coords, vals):
                    keep &= bits[:, j] == v
                ids = np.flatnonzero(keep)
                mask = sum(1 << int(i) for i in ids)
                group.append((
                    mask,
                    float(np.mean(rows[ids])),
                    int(np.max(rows[ids])),
                    coords,
                    vals,
                    len(ids),
                ))
        groups.append(group)
    return groups


def min_cylinder(cylinders, favorable: np.ndarray):
    """Minimum fixed-coordinate cylinder contained in favorable."""
    favmask = sum(1 << int(i) for i in np.flatnonzero(favorable))
    for k, group in enumerate(cylinders):
        candidates = [entry for entry in group if entry[0] & ~favmask == 0]
        if candidates:
            mask, meanrow, maxrow, coords, vals, size = min(
                candidates, key=lambda z: z[1:]
            )
            return k, meanrow, maxrow, coords, vals, size
    raise AssertionError("the favorable set is nonempty for every selector")


def audit(A: np.ndarray, name: str, m: int, t: float = 0.0):
    n = len(A)
    q = qnorm(A)
    p = m / n
    p2 = m * (m - 1) / (n * (n - 1))
    p2_exact = Fraction(m * (m - 1), n * (n - 1))
    t_exact = Fraction(t)
    p32 = p ** 1.5
    words = chart_words(n)
    raw = np.asarray([int(x @ A @ x) for _, _, x in words], dtype=np.int64)
    rows = np.asarray([int((A @ x) @ (A @ x)) for _, _, x in words],
                      dtype=np.int64)
    cylinders = all_cylinders(words, rows)
    codims = []
    means = []
    maxima = []
    sizes = []
    fav_counts = []
    strong_counts = []
    fav_row_means = []
    output = np.zeros(len(words), dtype=float)
    cover_counts = np.zeros(len(words), dtype=np.int64)
    boundary_margin = math.inf
    for S in itertools.combinations(range(n), m):
        AS = A[np.ix_(S, S)]
        qs = qnorm(AS)
        child = np.asarray([int(x[list(S)] @ AS @ x[list(S)])
                            for _, _, x in words], dtype=np.int64)
        h = np.asarray([
            qs - sigma * c + p2 * sigma * e - p32 * q
            for (_, sigma, _), c, e in zip(words, child, raw)
        ])
        boundary_margin = min(boundary_margin, float(np.min(np.abs(h - t))))
        favorable = np.asarray([
            leq_p32_q(
                Fraction(qs - sigma * int(c))
                + p2_exact * sigma * int(e) - t_exact,
                q, m, n,
            )
            for (_, sigma, _), c, e in zip(words, child, raw)
        ], dtype=bool)
        # The older local-deficit incidence uses only E_d <= q and is
        # therefore strong when delta <= B+t.  Favorable cuts outside it use
        # the retained full deficit genuinely.
        delta = np.asarray([qs - sigma * c
                            for (_, sigma, _), c in zip(words, child)])
        Ballow = (p32 - p2) * q
        strong = np.asarray([
            leq_p32_q(
                Fraction(int(dlt)) + p2_exact * q - t_exact,
                q, m, n,
            )
            for dlt in delta
        ], dtype=bool)
        full_energy = np.asarray([sigma * e for (_, sigma, _), e
                                  in zip(words, raw)])
        profile_favorable = (
            full_energy <= q + (Ballow + t - delta) / p2 + 1e-12
        )
        assert np.array_equal(favorable, profile_favorable)
        assert np.all(favorable[strong])
        cert = min_cylinder(cylinders, favorable)
        k, meanrow, maxrow, _, _, size = cert
        codims.append(k)
        means.append(meanrow)
        maxima.append(maxrow)
        sizes.append(size)
        nfav = int(np.sum(favorable))
        fav_counts.append(nfav)
        strong_counts.append(int(np.sum(strong)))
        fav_row_means.append(float(np.mean(rows[favorable])))
        output[favorable] += 1.0 / nfav
        cover_counts[favorable] += 1

    output /= len(codims)
    uniform_prob = 1.0 / len(words)
    output_kl = float(np.sum(output[output > 0]
                             * np.log(output[output > 0] / uniform_prob)))
    total_kl = float(np.mean([math.log(len(words) / f) for f in fav_counts]))
    output_row = float(output @ rows)
    annealed_output = cover_counts.astype(float) / sum(fav_counts)
    annealed_output_kl = float(np.sum(
        annealed_output[annealed_output > 0]
        * np.log(annealed_output[annealed_output > 0] / uniform_prob)
    ))
    annealed_output_row = float(annealed_output @ rows)
    annealed_mass = float(np.mean(fav_counts) / len(words))
    annealed_L = float(-math.log(annealed_mass))

    # For a k-physical-bit cylinder, the general conditional row bound is
    # n(n-1)+2qk.  Counting the orientation bit only overestimates k.
    generic = n * (n - 1) + 2 * q * max(codims)
    result = {
        "name": name,
        "n": n,
        "m": m,
        "q": q,
        "codim_hist": dict(sorted(Counter(codims).items())),
        "worst_codim": max(codims),
        "smallest_cylinder": min(sizes),
        "worst_best_mean_row": max(means),
        "worst_best_max_row": max(maxima),
        "generic_row_bound": generic,
        "decision_margin": boundary_margin,
        "favorable_count_range": (min(fav_counts), max(fav_counts)),
        "annealed_favorable_mass": annealed_mass,
        "annealed_L": annealed_L,
        "annealed_strong_local_mass": float(np.mean(strong_counts) / len(words)),
        "retained_only_share": float(
            (sum(fav_counts) - sum(strong_counts)) / sum(fav_counts)
        ),
        "total_incidence_K": total_kl,
        "fiber_output_KL": output_kl,
        "fiber_output_mean_row": output_row,
        "annealed_output_KL": annealed_output_kl,
        "annealed_output_mean_row": annealed_output_row,
        "worst_conditional_favorable_row": max(fav_row_means),
        "best_column_coverage": float(np.max(cover_counts) / len(codims)),
    }
    print(result)
    assert boundary_margin > 1e-9
    # In the x_0=1 chart, fixing an oriented local ground uses at most the
    # orientation bit and all m local bits when 0 is not in the selector.
    assert max(codims) <= m + 1
    assert min(sizes) == 2 ** (n - max(codims))
    assert output_kl <= total_kl + 1e-10
    assert annealed_output_kl <= annealed_L + 1e-10
    extracted = [
        i for i in range(len(words))
        if cover_counts[i] / len(codims) >= annealed_mass / 2 - 1e-12
        and rows[i] <= 2 * annealed_output_row + 1e-12
    ]
    assert extracted
    return result


def main():
    results = []
    for A, name, ms in (
        (A6, "A6", (3, 4, 5)),
        (A8, "A8", (4, 5, 6, 7)),
        (A9, "A9", (5, 6, 7, 8)),
    ):
        for m in ms:
            results.append(audit(A, name, m))
    # Nontriviality: retained tolerance enlarges at least one exact-ground
    # completion cylinder in the audited fixed-density cases.
    assert any(r["worst_codim"] < r["m"] for r in results)
    print("PASS cylinder-certificate enumeration")


if __name__ == "__main__":
    main()
