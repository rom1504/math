"""Replay finite isotropic-ground certificates and test Paley cap values.

No asymptotic inference is licensed by the finite cap values. Stored
optimality labels are imported, while each cap and covariance law is
checked independently here with integer/rational arithmetic.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import itertools
import json
from pathlib import Path
import subprocess

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def replay_stored():
    data = json.loads((ROOT / "computations/results/transfer_adversary_minimizer_isotropy_2026_09_06.json").read_text())
    records = []
    for case in data["cases"]:
        if case["minimum_uniform_slack_for_isotropy"] != 0:
            continue
        source = case["source"]
        matrix = json.loads((ROOT / source["path"]).read_text())
        for key in source["json_pointer"].strip("/").split("/"):
            matrix = matrix[int(key)] if isinstance(matrix, list) else matrix[key]
        matrix = np.asarray(matrix, dtype=np.int64)
        n = len(matrix)
        words = np.asarray([(1,) + row for row in itertools.product((-1, 1), repeat=n - 1)])
        energies = np.einsum("bi,ij,bj->b", words, matrix, words) // 2
        cap = int(max(abs(energies)))
        law = case["slack_screen"][-1]["law"]
        weights = [Fraction(atom["weight"]) for atom in law]
        assert sum(weights) == 1 and min(weights) >= 0
        for atom in law:
            x = np.asarray(atom["spin"], dtype=np.int64)
            assert abs(int(x @ matrix @ x)) == 2 * cap
        for i in range(n):
            for j in range(n):
                value = sum(w * atom["spin"][i] * atom["spin"][j]
                            for atom, w in zip(law, weights))
                assert value == int(i == j)
        records.append({"case": case["case"], "n": n, "cap": cap,
                        "normalized_cap": cap / n ** 1.5,
                        "isotropic_law_atoms": len(law),
                        "factorization_norm_squared": n,
                        "status": "exact cap and rational isotropic-ground law replayed"})
    return records


def paley(prime: int):
    assert prime % 4 == 1
    squares = {i * i % prime for i in range(1, prime)}
    n = prime + 1
    matrix = np.ones((n, n), dtype=np.int64) - np.eye(n, dtype=np.int64)
    for i in range(prime):
        for j in range(i):
            matrix[i, j] = matrix[j, i] = 1 if (i - j) % prime in squares else -1
    assert np.array_equal(matrix @ matrix, prime * np.eye(n, dtype=np.int64))
    return matrix


def paley_records(max_prime: int):
    directory = ROOT / "tmp/paper_portfolio_2026_09_17/discrepancy"
    directory.mkdir(parents=True, exist_ok=True)
    binary = directory / "exact_fixed_signing_gray"
    subprocess.run(["g++", "-O3", "-std=c++17",
                    str(ROOT / "computations/exact_fixed_signing_gray.cpp"),
                    "-o", str(binary)], check=True)
    records = []
    for prime in [5, 13, 17, 29]:
        if prime > max_prime:
            continue
        matrix = paley(prime)
        n = len(matrix)
        payload = str(n) + "\n" + "\n".join(" ".join(map(str, row)) for row in matrix) + "\n"
        result = subprocess.run([str(binary)], input=payload, capture_output=True,
                                text=True, check=True)
        record = json.loads(result.stdout)
        assert record["max_energy"] == -record["min_energy"]
        record["normalized_cap"] = record["cap"] / n ** 1.5
        record["factorization_norm_squared"] = n
        record["isotropy_source"] = "exact signed pair-transitive and anti-symmetry argument"
        record["scope"] = "finite example, not an asymptotic subhalf family"
        print(json.dumps({"paley_progress": record}), flush=True)
        records.append(record)
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-prime", type=int, default=17)
    options = parser.parse_args()
    report = {"stored_replay": replay_stored(),
              "paley": paley_records(options.max_prime),
              "asymptotic_isotropic_subhalf_claim": "open; not decided by this script"}
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
