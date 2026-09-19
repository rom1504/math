#!/usr/bin/env python3
"""Independent integer checks of the nonlinear diamond identity and clones.

The proof is algebraic. These checks use exact integer arithmetic for the
identity, all Boolean spins for the small seed/parent caps, and no optimizer.
"""

import itertools
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def cube(n):
    return 1 - 2 * ((np.arange(1 << n)[:, None] >> np.arange(n)) & 1)


def cap(a):
    x = cube(len(a))
    values = np.einsum("bi,ij,bj->b", x, a, x) // 2
    return int(np.max(np.abs(values)))


def parent_cap(a, b, diagonal):
    x = cube(len(a))
    energy = np.einsum("bi,ij,bj->b", x, a, x) // 2
    bridge = x @ (b + np.diag(diagonal)) @ x.T
    return int(np.max(np.abs(energy[:, None] - energy[None, :]) + np.abs(bridge)))


def main():
    rng = np.random.default_rng(20260919)
    denominator = 16
    identity_checks = 0
    bound_checks = 0
    for n in range(2, 13):
        for trial in range(400):
            upper = np.triu(rng.integers(-5, 6, (n, n)), 1)
            matrix = upper + upper.T
            if trial % 2:
                matrix += np.diag(rng.integers(-4, 5, n))
            b = rng.integers(0, denominator + 1, n)
            d = rng.integers(-denominator, denominator + 1, n)
            a = np.array([rng.integers(-(denominator - v), denominator - v + 1) for v in b])
            c = np.array([rng.integers(-(denominator - abs(v)), denominator - abs(v) + 1) for v in d])
            h = np.maximum(b + d - denominator, 0)
            u, v, w = b + d - h, b - h, d - h
            arguments = [u, v, a + h, a - h, c + w, c - w]
            assert max(int(np.max(np.abs(z))) for z in arguments) <= denominator

            def q(z):
                return int(z @ matrix @ z)

            twice_energy = q(a) - q(c) + 2 * int(b @ matrix @ d)
            twice_rhs = 2 * q(u) - 2 * q(v) + q(a + h) + q(a - h) - q(c + w) - q(c - w)
            assert 2 * twice_energy == twice_rhs
            identity_checks += 1
            if n <= 8 and trial % 2 == 0:
                assert abs(twice_energy) <= 8 * cap(matrix) * denominator**2
                bound_checks += 1

    clones = []
    for n in range(2, 5):
        for trial in range(8):
            upper = np.triu(rng.choice([-1, 1], (n, n)), 1)
            a = upper + upper.T
            twin = rng.choice([-1, 1], n)
            child = np.kron(a, np.ones((2, 2), dtype=int))
            for i, sign in enumerate(twin):
                child[2 * i, 2 * i + 1] = child[2 * i + 1, 2 * i] = sign
            switch = np.tile([1, -1], n)
            bridge = child * switch[:, None] * switch[None, :]
            diagonal = rng.choice([-1, 1], 2 * n)
            qa, qc, qp = cap(a), cap(child), parent_cap(child, bridge, diagonal)
            assert abs(qc - 4 * qa) <= n
            assert qp <= 16 * qa + 6 * n
            assert qp <= 4 * qc + 10 * n
            clones.append({
                "base_order": n, "trial": trial, "base_matrix": a.tolist(),
                "twin_signs": twin.tolist(), "parent_matching": diagonal.tolist(),
                "base_cap": qa, "child_cap": qc, "parent_cap": qp,
                "base_bound": 16 * qa + 6 * n, "child_bound": 4 * qc + 10 * n,
            })

    sharpness = []
    for m in [2, 4, 8, 16, 64]:
        qa, attained = m * m // 2, 2 * m * (m - 1)
        if m <= 4:
            block = np.ones((m, m), dtype=int) - np.eye(m, dtype=int)
            a = np.block([[block, np.zeros_like(block)], [np.zeros_like(block), -block]])
            assert cap(a) == qa
            b = np.ones(2 * m, dtype=int)
            d = np.r_[np.ones(m, dtype=int), -np.ones(m, dtype=int)]
            assert int(b @ a @ d) == attained
        sharpness.append({"block_size": m, "seed_cap": qa, "diamond_value": attained,
                          "ratio": attained / qa})

    result = {
        "normalization": "Q(A)=max_x |sum_(i<j) A_ij x_i x_j|; hollow matrices",
        "identity_integer_checks": identity_checks,
        "hollow_diamond_bound_checks": bound_checks,
        "all_checks_pass": True,
        "clone_records": clones,
        "sharpness_family": sharpness,
    }
    output = ROOT / "computations/results/twisted_chiral_diamond_2026_09_19.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: result[key] for key in ["identity_integer_checks", "hollow_diamond_bound_checks", "all_checks_pass"]}))
    print(f"Verified {len(clones)} completed clone witnesses; output {output.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
