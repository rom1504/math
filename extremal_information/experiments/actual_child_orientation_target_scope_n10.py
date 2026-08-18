#!/usr/bin/env python3
"""Reproduce the finite N=10 failure of exact balanced-orientation optimality.

This is a complete finite enumeration, not an asymptotic experiment.  The
child energies, histograms, signing cube, and bridge cube are enumerated
exactly.  Transcendental child pressures are compared with mpmath at 100
decimal digits; bridge pressures and soft values use the long-double Walsh
implementation audited in ``actual_child_bridge_law_exact.py``.

Run from the repository root with

    .venv/bin/python \
      extremal_information/experiments/actual_child_orientation_target_scope_n10.py

The assertions deliberately use margins much wider than floating error.
"""

from __future__ import annotations

import json
import math

import mpmath as mp
import numpy as np

import actual_child_bridge_law_exact as exact


TOTAL_ORDER = 10
BETA_TEXT = "4"
BETA = float(BETA_TEXT)
CHILD_ORDERS = (3, 7)
LAMBDAS = (2.0, 4.0)


def sector_bias(energies: np.ndarray, raw_t: mp.mpf) -> mp.mpf:
    count = int(len(energies))
    positive = mp.fsum(mp.exp(raw_t * int(value)) for value in energies) / count
    negative = mp.fsum(mp.exp(-raw_t * int(value)) for value in energies) / count
    return mp.log(positive / negative) / 2


def soft_value(pressure: np.ndarray, inverse_disorder: float) -> float:
    exponent = -inverse_disorder * np.asarray(pressure, dtype=np.float64)
    anchor = float(np.max(exponent))
    return -(
        anchor + math.log(float(np.mean(np.exp(exponent - anchor))))
    ) / inverse_disorder


def json_safe(value):
    """Replace the selector's ``inf`` no-competitor marker by JSON null."""
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if math.isfinite(float(value)) else None
    if isinstance(value, np.integer):
        return int(value)
    return value


def orientation_record(
    left: np.ndarray,
    right: np.ndarray,
    orientation: int,
) -> dict:
    pressure, audit = exact.bridge_pressures(
        left,
        right,
        BETA,
        TOTAL_ORDER,
        orientation,
    )
    order = np.argsort(pressure, kind="stable")
    minimum = float(pressure[order[0]])
    distinct = pressure[pressure > minimum + 1e-12]
    return {
        "orientation": orientation,
        "bridge_pressure_minimum": minimum,
        "minimizing_bridge_mask": int(order[0]),
        "minimum_multiplicity_at_1e-12": int(
            np.count_nonzero(np.abs(pressure - minimum) <= 1e-12)
        ),
        "next_distinct_bridge_pressure_at_1e-12": float(np.min(distinct)),
        "bridge_pressure_maximum": float(np.max(pressure)),
        "soft_values": {
            format(value, ".12g"): soft_value(pressure, value)
            for value in LAMBDAS
        },
        "bridge_audit": audit,
    }


def main() -> None:
    mp.mp.dps = 100
    raw_t = mp.mpf(BETA_TEXT) / mp.sqrt(TOTAL_ORDER)

    selected: dict[int, dict] = {}
    certificates: dict[int, dict] = {}
    matrices: dict[int, np.ndarray] = {}
    biases: dict[int, mp.mpf] = {}

    for order in CHILD_ORDERS:
        space = exact.build_signing_space(order)
        classes, certificate = exact.thermal_minimizer_classes(
            space,
            BETA_TEXT,
            TOTAL_ORDER,
        )
        assert len(classes) == 1
        representative = classes[0]
        matrix = np.asarray(representative["representative_matrix"], dtype=np.int8)
        energies = exact.energies_for_matrix(matrix, exact.projective_spins(order))
        selected[order] = {
            "representative_mask": representative["representative_mask"],
            "representative_sha256": representative["representative_sha256"],
            "root_gauged_member_count": representative["root_gauged_member_count"],
        }
        certificates[order] = certificate
        matrices[order] = matrix
        biases[order] = sector_bias(energies, raw_t)

    assert mp.mpf(certificates[7]["mp_gap_to_next_histogram"]) > mp.mpf("0.17")
    bias_product = biases[3] * biases[7]
    assert bias_product != 0
    balanced_orientation = -1 if bias_product > 0 else 1
    assert balanced_orientation == 1

    records = {
        orientation: orientation_record(
            matrices[3],
            matrices[7],
            orientation,
        )
        for orientation in (-1, 1)
    }
    balanced = records[balanced_orientation]
    opposite = records[-balanced_orientation]

    target = float(
        mp.mpf(certificates[3]["mp_optimum"])
        + mp.mpf(certificates[7]["mp_optimum"])
    )
    minimum_loss = (
        balanced["bridge_pressure_minimum"]
        - opposite["bridge_pressure_minimum"]
    )
    soft_losses = {
        key: balanced["soft_values"][key] - opposite["soft_values"][key]
        for key in balanced["soft_values"]
    }

    # Wide-margin regression assertions for the advertised scope check.
    assert minimum_loss > 0.52
    assert soft_losses["2"] > 0.02
    assert soft_losses["4"] > 0.30
    assert opposite["bridge_pressure_minimum"] - target > 3.0
    assert max(
        records[orientation]["bridge_audit"]["maximum_direct_log_pressure_error"]
        for orientation in records
    ) < 1e-12

    output = {
        "schema": "actual-child-orientation-target-scope-n10-v1",
        "classification": (
            "complete finite signing and bridge enumeration; integer-exact "
            "energies/histograms and numerical transcendental evaluation"
        ),
        "scope": {
            "N": TOTAL_ORDER,
            "split": list(CHILD_ORDERS),
            "beta": BETA,
            "raw_t": mp.nstr(raw_t, 100),
            "inverse_disorder_values": list(LAMBDAS),
        },
        "children": {
            str(order): {
                **selected[order],
                "sector_bias": mp.nstr(biases[order], 100),
                "minimizer_certificate": certificates[order],
            }
            for order in CHILD_ORDERS
        },
        "balanced_orientation": balanced_orientation,
        "same_raw_temperature_child_target": target,
        "orientation_records": {
            str(orientation): records[orientation]
            for orientation in (-1, 1)
        },
        "balanced_minus_opposite": {
            "minimum_bridge_pressure": minimum_loss,
            "soft_values": soft_losses,
        },
        "scope_limit": (
            "finite counterexample to exact orientation optimality only; "
            "neither orientation reaches the child target and no asymptotic "
            "orientation-loss conclusion is inferred"
        ),
    }
    print(json.dumps(json_safe(output), indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
