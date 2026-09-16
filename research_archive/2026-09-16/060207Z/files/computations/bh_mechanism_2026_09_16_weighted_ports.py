#!/usr/bin/env python3
"""Replay weighted reciprocal-port bounds on actual repaired sign frames."""

import argparse
import json
from fractions import Fraction
from itertools import combinations
from math import isqrt, sqrt
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, minimize_scalar

from principle_synthesis_2026_09_07_sparse_active_check import walsh


def root_bounds(value, digits=40):
    scale = 10**digits
    floor = isqrt(value.numerator * scale * scale // value.denominator)
    return Fraction(floor, scale), Fraction(floor + 1, scale)


def rational_phase_certificate():
    p = Fraction(31, 32)
    theta = Fraction(999, 1000)
    a_lo, _ = root_bounds((1 - p) * theta)
    _, b_hi = root_bounds(p * (1 - theta))
    p_lo, _ = root_bounds(p)
    assert a_lo > b_hi
    deficit_lo = (a_lo - b_hi)**2
    coefficient_hi = (1 - deficit_lo) / (2 * p_lo)
    assert coefficient_hi < Fraction(497237, 1000000)
    assert coefficient_hi < Fraction(498, 1000)
    p_float = float(p)
    target = 0.499
    delta = 1 - 2 * target * sqrt(p_float)
    critical = brentq(lambda t: (sqrt((1 - p_float) * t)
                                - sqrt(p_float * (1 - t)))**2 - delta,
                      p_float, 1)
    return {"p": str(p), "core_variance_fraction": str(theta),
            "exact_rational_coefficient_upper": str(coefficient_hi),
            "decimal_upper": float(coefficient_hi),
            "certified_less_than": "497237/1000000",
            "eventual_actual_coefficient": "498/1000",
            "numeric_threshold_for_point_499": critical,
            "threshold_qualification": "The threshold is diagnostic; the 999/1000 certificate is rational."}


def make_frames(m, rng):
    q = 3 * m // 4
    had = walsh(m)
    frames = []
    for i in range(m):
        selected = []
        for group in range(m // 4):
            selected.extend(group * 4 + rng.choice(4, 3, replace=False))
        frame = had[np.asarray(selected)][:, np.arange(m) ^ i].copy()
        frame[:, i] = 0
        for j in range(m):
            excess = int(frame[:, j].sum())
            if excess:
                candidates = np.flatnonzero(frame[:, j] == np.sign(excess))
                frame[rng.choice(candidates, abs(excess) // 2, replace=False), j] *= -1
        assert np.all(frame.sum(axis=0) == 0)
        frames.append(frame)
    return q, frames


def exact_model_checks(seed=2026091605):
    rng = np.random.default_rng(seed)
    reports = []
    for m in (8, 16, 32, 64):
        q, frames = make_frames(m, rng)
        global_cap = max(float(np.linalg.eigvalsh(f.T @ f).max()) for f in frames) + 1e-9
        signs = np.triu(rng.choice((-1, 1), (m, m)), 1)
        signs = signs + signs.T
        checks = 0
        min_slack = float("inf")
        improved_cases = 0
        for _ in range(120):
            words = np.ones((m, q), dtype=np.int64)
            for i in range(m):
                # A broad collection of actual Boolean minority densities.
                r = rng.choice((0, 1 / q, 2 / q, 0.125, 0.25, 0.5))
                count = int(round(r * q))
                if count:
                    words[i, rng.choice(q, count, replace=False)] = -1
            variance = q - words.sum(axis=1).astype(float)**2 / q
            if not variance.sum():
                continue
            core_size = int(rng.integers(1, min(6, m) + 1))
            core = np.argsort(-variance)[:core_size]
            restricted_cap = max(float(np.linalg.eigvalsh(f[:, core].T @ f[:, core]).max())
                                 for f in frames) + 1e-9
            restricted_cap = min(restricted_cap, global_cap)
            v_h = float(variance[core].sum())
            v_l = float(variance.sum() - v_h)
            ports = np.asarray([frames[i].T @ words[i] for i in range(m)])
            energy = float(np.sum(np.triu(signs * ports * ports.T, 1)))
            gap = max(0.0, sqrt(max(0.0, (global_cap - restricted_cap) * v_h))
                      - sqrt(restricted_cap * v_l))
            bound = (global_cap * (v_h + v_l) - gap * gap) / 2
            assert abs(energy) <= bound + 1e-7
            min_slack = min(min_slack, bound - abs(energy))
            improved_cases += gap > 0
            if v_l > 0 and restricted_cap > 0:
                optimum = max(1.0, sqrt((global_cap - restricted_cap) * v_h
                                        / (restricted_cap * v_l)))
                weighted_bound = ((global_cap + (optimum - 1) * restricted_cap)
                                  * (v_h / optimum + v_l) / 2)
                assert abs(weighted_bound - bound) <= 1e-7
            for t in (1.0, 1.1, 2.0, 5.0):
                lam = np.ones(m)
                lam[core] = t
                weighted_square = float(np.sum(ports**2 * lam[None, :] / lam[:, None]))
                row_cap = (global_cap + (t - 1) * restricted_cap) * (v_h / t + v_l)
                assert 2 * abs(energy) <= weighted_square + 1e-7
                assert weighted_square <= row_cap + 1e-7
                checks += 1
        reports.append({"m": m, "q": q, "computed_global_column_cap": global_cap,
                        "weighted_checks": checks, "positive_deficit_cases": int(improved_cases),
                        "minimum_energy_slack": min_slack})
    return reports


def potential_optimum_checks(seed=2026091606):
    rng = np.random.default_rng(seed)
    largest_error = 0.0
    for _ in range(500):
        full = rng.uniform(0.2, 5)
        restricted = full * rng.uniform(0.02, 0.98)
        vh, vl = rng.uniform(0.01, 4, size=2)
        bound = (full * (vh + vl)
                 - max(0, sqrt((full - restricted) * vh) - sqrt(restricted * vl))**2) / 2
        result = minimize_scalar(lambda z: (full + (np.exp(z) - 1) * restricted)
                                 * (vh / np.exp(z) + vl) / 2,
                                 bounds=(0, 15), method="bounded",
                                 options={"xatol": 1e-12})
        error = abs(result.fun - bound)
        largest_error = max(largest_error, error)
        assert error <= 1e-7
    return {"numeric_checks": 500, "largest_error": largest_error}


def all_temperature_packet_checks():
    checked = 0
    for order in (4, 8):
        had = walsh(order)
        temperatures = (Fraction(1, 10), Fraction(1, 2), Fraction(1),
                        Fraction(order, 2), Fraction(order - 1), Fraction(order),
                        Fraction(2 * order), Fraction(10 * order))
        for c in temperatures:
            expected = min(1 / c, Fraction(2, order), order / c - 1)
            minimum = None
            for mask in range(1, (1 << order) - 1):
                subset = [j for j in range(order) if (mask >> j) & 1]
                coefficients = had[:, subset].sum(axis=1)
                observed = (sum(min(Fraction(int(v * v), 1) / c, 1)
                                for v in coefficients) - 1) / len(subset)
                assert observed >= expected
                minimum = observed if minimum is None else min(minimum, observed)
                checked += 1
            assert minimum == expected
    rng = np.random.default_rng(2026091607)
    for _ in range(1000):
        order = 32
        retained = order - 1
        weights = rng.dirichlet(np.ones(5))
        tau = np.exp(rng.uniform(-5, 5, size=5))

        def gamma(z):
            return max(-4 * z / order, -2 * retained / (order * order),
                       retained / order - 4 * z)

        assert sum(w * gamma(t) for w, t in zip(weights, tau)) >= gamma(weights @ tau) - 1e-12
    return {"exact_packet_temperature_tests": checked, "temperature_jensen_tests": 1000,
            "c_scope": "All positive c; replay includes c below 1 and above L."}


def actual_column_obstructions(seed=2026091608):
    rng = np.random.default_rng(seed)
    rank_reports = []
    for m in (8, 16):
        q, frames = make_frames(m, rng)
        frame = frames[0][:, 1:]
        d = frame.shape[1]
        gram = frame.T @ frame
        for s in (2, 3):
            lower = q + (s - 1) * q * (d - q + 1) / ((d - 1) * (q - 1))
            maximum = max(float(np.linalg.eigvalsh(gram[np.ix_(subset, subset)]).max())
                          for subset in combinations(range(d), s))
            assert maximum >= lower - 1e-10
            rank_reports.append({"m": m, "q": q, "s": s,
                                 "rank_trace_lower": lower, "actual_maximum": maximum})
    coset_reports = []
    for child_order, coset_size in ((64, 4), (256, 8)):
        node = 4
        ambient = node * child_order
        q = (node - 1) * child_order
        had = walsh(ambient)
        for attempt in range(20):
            omitted = rng.integers(node, size=child_order)
            block = next((b for b in range(child_order // coset_size)
                          if np.all(omitted[b * coset_size:(b + 1) * coset_size] != 0)), None)
            if block is not None:
                break
        assert block is not None
        selected = np.asarray([a * child_order + g for a in range(node)
                               for g in range(child_order) if omitted[g] != a])
        assert len(selected) == q
        physical = np.zeros(ambient)
        support = block * coset_size + np.arange(coset_size)
        physical[support] = np.where(np.arange(coset_size) % 2, -1, 1) / sqrt(coset_size)
        assert abs(physical.sum()) < 1e-12
        coefficient = had.T @ physical / sqrt(ambient)
        columns = np.flatnonzero(np.abs(coefficient) > 1e-10)
        assert len(columns) == ambient // coset_size
        assert 0 not in columns
        assert abs(float(coefficient @ coefficient) - 1) < 1e-12
        raw = had[selected].copy()
        image = raw @ coefficient
        assert abs(float(image @ image) - ambient) < 1e-8
        assert abs(float(image.sum())) < 1e-8
        repaired = raw.copy()
        repaired[:, 0] = 0
        for j in range(1, ambient):
            excess = int(repaired[:, j].sum())
            if excess:
                candidates = np.flatnonzero(repaired[:, j] == np.sign(excess))
                repaired[rng.choice(candidates, abs(excess) // 2, replace=False), j] *= -1
        assert np.all(repaired.sum(axis=0) == 0)
        repaired_image = repaired[:, columns] @ coefficient[columns]
        repaired_rayleigh = float(repaired_image @ repaired_image)
        coset_reports.append({"ambient": ambient, "physical": q, "coset_size": coset_size,
                              "column_support": len(columns), "raw_centered_squared_norm": ambient,
                              "repaired_sign_frame_witness_squared_norm": repaired_rayleigh,
                              "sampling_attempts": attempt + 1,
                              "repair_qualification": "Finite replay uses exact random balancing; asymptotic survival uses the separately proved operator-controlled repair."})
    return {"rank_trace": rank_reports, "retained_affine_cosets": coset_reports}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = {"status": "PASS", "rational_phase": rational_phase_certificate(),
              "actual_sign_frame_checks": exact_model_checks(),
              "potential_optimization": potential_optimum_checks(),
              "all_temperature_packet_curve": all_temperature_packet_checks(),
              "actual_column_obstructions": actual_column_obstructions()}
    rendered = json.dumps(report, indent=2, sort_keys=True)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n")
    print(rendered)


if __name__ == "__main__":
    main()
