#!/usr/bin/env python3
"""Wave 49 exact audits for project-row coarea and completion mass.

All finite calculations are exhaustive for each matrix passed to ``audit``.
The order-ten matrices themselves are deterministic MILP samples, not an
exhaustive list of exact minimizers.
"""

from __future__ import annotations

import itertools
import math
import sys
from collections import Counter
from fractions import Fraction

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from coarea_alignment_r47_check import exact_boundary_data
from cross_block_dual_r48_check import sample_order_ten
from envelope_block_cover_r27 import A6, A8, A9, projective_spins
from hard_center_branching_r45_check import down_up_kernel


def quadratic_energy(a: np.ndarray, x: np.ndarray) -> int:
    return int(x @ a @ x)


def row_square(a: np.ndarray, x: np.ndarray) -> int:
    return int(np.sum((a @ x) ** 2))


def exact_center_terms(f: np.ndarray, kernel, ell: int, n: int, m: int):
    """Return (a,B,collision,Xi) for every center, using exact rationals."""
    nsel = f.shape[1]
    lam = Fraction(ell * (n - m), m * (n - ell)) if ell else Fraction(0)
    delta = 1 - lam
    out = []
    for row in f:
        active = np.flatnonzero(row)
        a = Fraction(len(active), nsel)
        collision = sum(
            (kernel[int(i), int(j)] for i in active for j in active), Fraction(0)
        ) / nsel
        b = a - collision
        # B = delta*(a-a^2) + Xi exactly.
        xi = b - delta * (a - a * a)
        assert xi >= 0
        out.append((a, b, collision, xi))
    return lam, delta, out


def aggregate_metrics(terms, keep, lam: Fraction, delta: Fraction):
    """Exact coarea, cubic-harmonic, and triple-retention normalizations."""
    selected = [terms[i] for i in range(len(terms)) if keep[i]]
    d2 = sum((a * a for a, _, _, _ in selected), Fraction(0))
    if not d2:
        return None
    lhs = sum((a * b for a, b, _, _ in selected), Fraction(0))
    cubic_lhs = sum((a * xi for a, _, _, xi in selected), Fraction(0))
    d3 = sum((a**3 for a, _, _, _ in selected), Fraction(0))
    triple = sum((a * collision for a, _, collision, _ in selected), Fraction(0))
    coarea = lhs / (delta * d2)
    cubic = cubic_lhs / (delta * d3) if d3 else None
    retention = triple / (lam * d2) if lam else None
    # Three exactly equivalent forms of the non-strict target.
    assert (coarea <= 1) == (cubic <= 1)
    if lam:
        assert (coarea <= 1) == (retention >= 1)
        assert lhs + triple == d2
    assert lhs - delta * d2 == cubic_lhs - delta * d3
    return coarea, cubic, retention, d2, d3


def completion_pareto(a: np.ndarray, m: int):
    """Exhaust all selectors and child grounds and audit completion formulas."""
    n = len(a)
    q = max(abs(quadratic_energy(a, np.asarray(x))) for x in projective_spins(n))
    records = []
    identities = 0
    for s in itertools.combinations(range(n), m):
        s = tuple(s)
        t = tuple(i for i in range(n) if i not in s)
        child = a[np.ix_(s, s)]
        ys = np.asarray(list(projective_spins(m)), dtype=np.int64)
        energies = np.einsum("bi,ij,bj->b", ys, child, ys, optimize=True).astype(int)
        qs = int(np.max(np.abs(energies)))
        delta_q = q - qs
        for y, e in zip(ys, energies):
            if abs(int(e)) != qs:
                continue
            internal_vec = child @ y
            external_vec = a[np.ix_(t, s)] @ y
            gi = int(internal_vec @ internal_vec)
            ge = int(external_vec @ external_vec)
            k = len(t)
            mean_row = Fraction(gi + ge + k * (n - 1))
            mean_energy_sq = Fraction(qs * qs + 4 * ge + 2 * k * (k - 1))
            assert mean_energy_sq <= q * q
            certificate = Fraction(gi) + Fraction(delta_q * (2 * q - delta_q), 4)
            baseline = Fraction(k * (n - 1)) - Fraction(k * (k - 1), 2)
            upper_mean = certificate + baseline
            assert mean_row <= upper_mean

            rows = []
            for w in itertools.product((-1, 1), repeat=k):
                z = np.empty(n, dtype=np.int64)
                z[list(s)] = y
                z[list(t)] = w
                rows.append(row_square(a, z))
            assert Fraction(sum(rows), len(rows)) == mean_row
            min_row = min(rows)
            records.append(
                {
                    "S": s,
                    "Q_S": qs,
                    "Delta": delta_q,
                    "G_int": gi,
                    "G_ext": ge,
                    "P": certificate,
                    "mean_row": mean_row,
                    "upper_mean": upper_mean,
                    "min_row": min_row,
                }
            )
            identities += 1

    def key(rec, fields):
        return tuple(rec[x] for x in fields)

    # These minima are over all child-ground incidences, exhaustively.
    best_p = min(records, key=lambda r: key(r, ("P", "mean_row", "min_row")))
    best_mean = min(records, key=lambda r: key(r, ("mean_row", "P", "min_row")))
    best_exact = min(records, key=lambda r: key(r, ("min_row", "P", "mean_row")))
    best_short = min(records, key=lambda r: key(r, ("Delta", "G_int", "mean_row")))
    return q, identities, best_p, best_mean, best_exact, best_short


def compact_record(r):
    return {
        x: r[x]
        for x in ("S", "Q_S", "Delta", "G_int", "G_ext", "P", "mean_row", "upper_mean", "min_row")
    }


def audit(name: str, a: np.ndarray, m: int) -> None:
    n = len(a)
    z, selectors, q, parent_e, rows, deficits, h0, f, degree = exact_boundary_data(a, m)
    assert np.array_equal(np.asarray([row_square(a, x) for x in z]), rows)
    active_rows = sorted({int(rows[i]) for i, az in enumerate(degree) if az > 0})
    old_cap = 2 * n * (n - 1)
    print(f"\n{name}: n={n},m={m},q={q}, centers={len(z)}, selectors={len(selectors)}")
    print(
        f"  active center rows: min={min(active_rows)}, max={max(active_rows)}, "
        f"old_cap={old_cap}, active_under_old_cap={sum(degree[i]>0 and rows[i]<=old_cap for i in range(len(z)))}"
    )
    success_pattern = []
    key_metrics = {}
    for ell in range(1, m):
        kernel = down_up_kernel(selectors, n, m, ell)
        lam, delta, terms = exact_center_terms(f, kernel, ell, n, m)
        all_metrics = aggregate_metrics(terms, np.ones(len(z), dtype=bool), lam, delta)
        old_metrics = aggregate_metrics(terms, rows <= old_cap, lam, delta)
        # Full row-threshold phase diagram, summarized by success intervals.
        phases = []
        last = None
        for cap in active_rows:
            metrics = aggregate_metrics(terms, rows <= cap, lam, delta)
            success = metrics is not None and metrics[0] <= 1
            if success != last:
                phases.append((cap, success, None if metrics is None else metrics[0]))
                last = success
        ca, cc, cr, _, _ = all_metrics
        success_pattern.append(ca <= 1)
        key_metrics[ell] = (ca, cc, cr)
        oa, oc, orr, _, _ = old_metrics
        print(
            f"  ell={ell}: lambda={lam}, all(coarea,cubic,triple)=({ca},{cc},{cr}); "
            f"oldcap=({oa},{oc},{orr}); phase_changes={phases}"
        )

    pq, count, best_p, best_mean, best_exact, best_short = completion_pareto(a, m)
    assert pq == q
    print(f"  completion-ground incidences exhausted={count}")
    print(f"    min joint P: {compact_record(best_p)}")
    print(f"    min mean row: {compact_record(best_mean)}")
    print(f"    min exact completed row: {compact_record(best_exact)}")
    print(f"    min shortfall then internal Gram: {compact_record(best_short)}")
    expected_success = {
        "A6": (True, True, True, True),
        "A8": (True, True, True, True, True),
        "A9": (True, True, True, True, True, False),
        "A10_sample_seed_0": (True, True, True, True, False, False, False),
    }
    assert tuple(success_pattern) == expected_success[name]
    expected_p = {"A6": 33, "A8": 73, "A9": 107, "A10_sample_seed_0": 105}
    assert best_p["P"] == expected_p[name]
    if name == "A9":
        assert key_metrics[6] == (
            Fraction(10154, 9891), Fraction(10541, 9489), Fraction(12925, 13188)
        )
    if name == "A10_sample_seed_0":
        assert key_metrics[5] == (
            Fraction(11369, 11200), Fraction(16025, 14504), Fraction(10693, 11200)
        )


def order_ten_sample_summary(seeds: int = 12) -> None:
    """Summarize exhaustive-within-matrix audits of sampled exact minimizers."""
    patterns = Counter()
    pareto = Counter()
    codes = set()
    for seed in range(seeds):
        a = sample_order_ten(seed)
        codes.add(tuple(int(a[i, j]) for i in range(1, 10) for j in range(i + 1, 10)))
        z, selectors, q, _, rows, _, _, f, degree = exact_boundary_data(a, 8)
        assert q == 26
        ratios = []
        for ell in range(1, 8):
            kernel = down_up_kernel(selectors, 10, 8, ell)
            lam, delta, terms = exact_center_terms(f, kernel, ell, 10, 8)
            metrics = aggregate_metrics(terms, np.ones(len(z), dtype=bool), lam, delta)
            ratios.append(metrics[0])
        patterns[tuple(r <= 1 for r in ratios)] += 1
        _, _, best_p, _, best_exact, best_short = completion_pareto(a, 8)
        pareto[(best_p["P"], best_p["Delta"], best_p["G_int"], best_p["G_ext"], best_exact["min_row"])] += 1
        assert max(int(rows[i]) for i, az in enumerate(degree) if az > 0) <= 180
    print("\nA10 deterministic MILP sample summary (matrix sample is nonexhaustive):")
    print(f"  sampled={seeds}, distinct_normalized_matrices={len(codes)}")
    print(f"  coarea success patterns ell=1..7: {dict(patterns)}")
    print(f"  joint-P/min-row profiles: {dict(pareto)}")
    assert patterns == Counter({(True, True, True, True, False, False, False): 12})
    assert pareto == Counter({(Fraction(105), 2, 80, 0, 74): 7,
                              (Fraction(113), 2, 88, 0, 74): 5})


def fixed_density_mass_diagnostics() -> None:
    """Audit lower-ratio mass, where near-diagonal shortfall evidence is irrelevant."""
    cases = [
        ("A6", A6, 4),
        ("A8", A8, 4),
        ("A8", A8, 5),
        ("A9", A9, 5),
        ("A9", A9, 6),
        ("A10_seed_0", sample_order_ten(0), 6),
    ]
    reports = {}
    print("\nLower-ratio exhaustive-within-matrix mass diagnostics:")
    for name, a, m in cases:
        q, count, best_p, best_mean, best_exact, _ = completion_pareto(a, m)
        report = (
            q,
            count,
            best_p["P"],
            best_p["Delta"],
            best_mean["G_int"] + best_mean["G_ext"],
            best_mean["Delta"],
            best_exact["min_row"],
            best_exact["Delta"],
        )
        reports[(name, m)] = report
        print(
            f"  {name},m={m}: q={q}, incidences={count}, "
            f"best_P={best_p['P']} (Delta={best_p['Delta']}), "
            f"best_Gtot={best_mean['G_int'] + best_mean['G_ext']} "
            f"(Delta={best_mean['Delta']}), exact_min_row={best_exact['min_row']} "
            f"(Delta={best_exact['Delta']})"
        )
    assert reports[("A9", 5)] == (24, 595, 136, 8, 28, 16, 16, 16)
    assert reports[("A10_seed_0", 6)] == (26, 520, 134, 4, 38, 16, 10, 16)


def main() -> None:
    audit("A6", A6, 5)
    audit("A8", A8, 6)
    audit("A9", A9, 7)
    audit("A10_sample_seed_0", sample_order_ten(0), 8)
    order_ten_sample_summary()
    fixed_density_mass_diagnostics()
    print("\nPASS: cubic/triple equivalence, row-cap phase diagrams, and completion Pareto formulas")


if __name__ == "__main__":
    main()
