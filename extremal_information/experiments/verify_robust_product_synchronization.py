#!/usr/bin/env python3
"""Finite adversarial checks for robust product synchronization.

The theorem is analytic; this script searches small exact contraction/port
instances for violations of the positivity, response, Schur, and tensor
claims.  It deliberately includes indefinite Rayleigh kernels.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product
from math import prod
import random

import numpy as np

from verify_boolean_port_product_algebra_closure import h16_seed


TOL = 2.0e-9


def all_involutions(n: int):
    """All involutions as tuples."""
    def rec(unused, pairs):
        if not unused:
            pi = list(range(n))
            for i, j in pairs:
                pi[i] = j
                pi[j] = i
            yield tuple(pi)
            return
        i = min(unused)
        rest = set(unused)
        rest.remove(i)
        yield from rec(rest, pairs + [(i, i)])
        for j in sorted(rest):
            tail = set(rest)
            tail.remove(j)
            yield from rec(tail, pairs + [(i, j)])
    yield from rec(set(range(n)), [])


def signed_symmetric_permutations(n: int):
    """Exact symmetric orthogonal matrices with entries in {0,+-1}."""
    out = []
    for pi in all_involutions(n):
        orbits = []
        seen = set()
        for i in range(n):
            if i not in seen:
                orbit = tuple(sorted({i, pi[i]}))
                seen.update(orbit)
                orbits.append(orbit)
        for signs in product((-1, 1), repeat=len(orbits)):
            a = np.zeros((n, n), dtype=float)
            for orbit, sign in zip(orbits, signs):
                if len(orbit) == 1:
                    a[orbit[0], orbit[0]] = sign
                else:
                    i, j = orbit
                    a[i, j] = a[j, i] = sign
            assert np.array_equal(a, a.T)
            assert np.array_equal(a @ a, np.eye(n))
            out.append(a)
    return out


def selector_table(p: int):
    """Odd majority with first-coordinate tie rule."""
    table = {}
    for a in product((-1, 1), repeat=p):
        total = sum(a)
        table[a] = (1 if total > 0 else -1 if total < 0 else a[0])
    return table


def fourier_support(p: int):
    table = selector_table(p)
    support = []
    for mask in range(1 << p):
        val = sum(
            Fraction(table[a] * prod(a[i] for i in range(p) if mask >> i & 1), 1 << p)
            for a in table
        )
        if val:
            support.append((mask, val))
            assert bin(mask).count("1") % 2 == 1
    assert sum(v * v for _, v in support) == 1
    return support


def active_matrix(base_ports: np.ndarray, support):
    n, p = base_ports.shape
    cols = []
    for mask, _ in support:
        col = np.ones(n)
        for i in range(p):
            if mask >> i & 1:
                col *= base_ports[:, i]
        cols.append(col)
    return np.column_stack(cols)


def state(t: np.ndarray, v: np.ndarray):
    n = v.shape[0]
    g = v.T @ v / n
    r = v.T @ t @ v / n
    return g, r, g - r


def endpoint_coeff(eps, support):
    return np.array([
        float(hat) * prod(eps[i] for i in range(len(eps)) if mask >> i & 1)
        for mask, hat in support
    ])


def exact_response(t, base_ports, eps, m):
    n = len(t)
    z = base_ports @ np.asarray(eps, dtype=float)
    best = -float("inf")
    for x_tuple in product((-1, 1), repeat=n):
        x = np.asarray(x_tuple, dtype=float)
        q = float(x @ t @ x) / 2
        field = m * float(z @ x)
        best = max(best, q + field, -q + field)
    return best, n / 2 + m * np.sum(np.abs(z))


def assert_psd(a, label):
    least = np.linalg.eigvalsh((a + a.T) / 2)[0]
    assert least >= -TOL, (label, least, a)


def check_growing_affine_coset_boundary():
    """Raw defect grows in PC.3; Gram-relative defect does not."""
    r, h, ports, _ = h16_seed()
    t0 = h / r
    support3 = fourier_support(3)
    v0 = active_matrix(ports.T, support3)
    g0, ray0, d0 = state(t0, v0)
    expected_g0 = np.asarray([
        [1, 0.5, 0, -0.5],
        [0.5, 1, -0.5, 0],
        [0, -0.5, 1, 0.5],
        [-0.5, 0, 0.5, 1],
    ])
    assert np.max(np.abs(g0 - expected_g0)) < TOL
    assert np.max(np.abs(ray0 - g0)) < TOL
    assert np.max(np.abs(d0)) < TOL
    assert abs(np.linalg.norm(g0, 2) - 2) < TOL

    # Check that the level-two PC.3 products are exactly the Cartesian
    # products of the seed active poles.
    a, b, c = ports
    one = np.ones(len(a))
    c1, c2 = a * b, a * c
    a2 = np.kron(a, a)
    generators = [
        np.kron(c1, one), np.kron(c2, one),
        np.kron(one, c1), np.kron(one, c2),
    ]
    ports2 = np.asarray([a2] + [a2 * g for g in generators]).T
    v2_pc = active_matrix(ports2, fourier_support(5))
    v2_cart = np.column_stack([
        np.kron(v0[:, i], v0[:, j])
        for i in range(v0.shape[1]) for j in range(v0.shape[1])
    ])
    assert {tuple(v2_pc[:, i]) for i in range(v2_pc.shape[1])} == {
        tuple(v2_cart[:, i]) for i in range(v2_cart.shape[1])
    }

    # Lower the active seed span by eta while leaving a contraction.
    projection = v0 @ np.linalg.pinv(v0)
    assert np.linalg.norm(projection @ projection - projection, 2) < TOL
    eta = 0.125
    t_eta = t0 - eta * projection
    assert np.linalg.norm(t_eta, 2) <= 1 + TOL
    g_eta, r_eta, d_eta = state(t_eta, v0)
    assert np.max(np.abs(r_eta - (1 - eta) * g0)) < TOL
    assert np.max(np.abs(d_eta - eta * g0)) < TOL

    # Actual level-two state, not merely a formal matrix calculation.
    t2 = np.kron(t_eta, t0)
    g2, r2, d2 = state(t2, v2_cart)
    assert np.max(np.abs(g2 - np.kron(g0, g0))) < TOL
    assert np.max(np.abs(d2 - eta * g2)) < 2e-8
    assert abs(np.linalg.norm(d_eta, 2) - 2 * eta) < TOL
    assert abs(np.linalg.norm(d2, 2) - 4 * eta) < 2e-8

    # Algebraic continuation through growing arities p_j=2j+1.
    raw = []
    for j in range(1, 5):
        gj = g0.copy()
        for _ in range(j - 1):
            gj = np.kron(gj, g0)
        dj = eta * gj
        raw.append(np.linalg.norm(dj, 2))
        assert abs(raw[-1] - eta * (2 ** j)) < 2e-8
        # Generalized D <= delta G is exact with delta=eta.
        assert_psd(eta * gj - dj, "relative defect equality")

    # Random contractions verify Cartesian relative subadditivity using the
    # compressed operators directly.  Here the active spans may be singular.
    relative_checks = 0
    rng = np.random.default_rng(431)
    for _ in range(100):
        dims = (rng.integers(1, 5), rng.integers(1, 5))
        compressed = []
        deltas = []
        for dim in dims:
            q, _ = np.linalg.qr(rng.normal(size=(dim, dim)))
            eig = rng.uniform(-1, 1, size=dim)
            aa = q @ np.diag(eig) @ q.T
            compressed.append(aa)
            deltas.append(1 - np.min(eig))
        a1, a2c = compressed
        delta12 = 1 - np.linalg.eigvalsh(np.kron(a1, a2c))[0]
        assert delta12 <= sum(deltas) + 2e-8
        relative_checks += 1
    return raw, relative_checks


def check_endpoint_table_noncongruence():
    """Equal p=3 query defects split under corresponding-port self tensor."""
    # Four uniform projective row types; the active monomials are the full
    # four-character Walsh matrix.
    w1 = np.ones(4)
    w2 = np.asarray([1, 1, -1, -1], dtype=float)
    w3 = np.asarray([1, -1, 1, -1], dtype=float)
    base = np.column_stack([w1, w2, w3])
    support = fourier_support(3)
    v = active_matrix(base, support)
    assert np.max(np.abs(v.T @ v / 4 - np.eye(4))) < TOL

    coeffs = []
    for eps in product((-1, 1), repeat=3):
        c = endpoint_coeff(eps, support)
        if not any(np.max(np.abs(c - old)) < TOL or
                   np.max(np.abs(c + old)) < TOL for old in coeffs):
            coeffs.append(c)
    cmat = np.column_stack(coeffs)
    assert len(coeffs) == 4
    assert np.max(np.abs(cmat.T @ cmat - np.eye(4))) < TOL

    a, t = 0.5, 0.25
    b = a * np.eye(4)
    b[0, 1] = b[1, 0] = t
    d0 = a * np.eye(4)
    d1 = cmat @ b @ cmat.T
    assert_psd(d0, "table D0")
    assert_psd(2 * np.eye(4) - d0, "table 2I-D0")
    assert_psd(d1, "table D1")
    assert_psd(2 * np.eye(4) - d1, "table 2I-D1")
    initial0 = np.diag(cmat.T @ d0 @ cmat)
    initial1 = np.diag(cmat.T @ d1 @ cmat)
    assert np.max(np.abs(initial0 - 0.5)) < TOL
    assert np.max(np.abs(initial1 - 0.5)) < TOL

    # Realize the two kernels by honest symmetric contractions on R^4.
    u = v / 2
    t0 = u @ (np.eye(4) - d0) @ u.T
    t1 = u @ (np.eye(4) - d1) @ u.T
    assert np.linalg.norm(t0, 2) <= 1 + TOL
    assert np.linalg.norm(t1, 2) <= 1 + TOL
    _, _, realized0 = state(t0, v)
    _, _, realized1 = state(t1, v)
    assert np.max(np.abs(realized0 - d0)) < TOL
    assert np.max(np.abs(realized1 - d1)) < TOL

    # Formula and actual tensor states agree.
    out0 = np.eye(4) - (np.eye(4) - d0) * (np.eye(4) - d0)
    out1 = np.eye(4) - (np.eye(4) - d1) * (np.eye(4) - d1)
    base2 = np.column_stack([np.kron(base[:, i], base[:, i]) for i in range(3)])
    v2 = active_matrix(base2, support)
    _, _, actual0 = state(np.kron(t0, t0), v2)
    _, _, actual1 = state(np.kron(t1, t1), v2)
    assert np.max(np.abs(actual0 - out0)) < TOL
    assert np.max(np.abs(actual1 - out1)) < TOL
    output0 = np.diag(cmat.T @ out0 @ cmat)
    output1 = np.diag(cmat.T @ out1 @ cmat)
    assert np.max(np.abs(output0 - 48 / 64)) < TOL
    assert np.max(np.abs(output1 - 47 / 64)) < TOL
    return initial0, output0, output1


def main():
    rng = random.Random(20260817)
    np_rng = np.random.default_rng(20260817)

    # Exact rational contractions: averages of two signed involutions.
    pools = {n: signed_symmetric_permutations(n) for n in (3, 4)}
    response_checks = 0
    schur_checks = 0
    tensor_checks = 0
    indefinite_cases = 0

    instances = []
    for p in (2, 3, 4):
        support = fourier_support(p)
        for n in (3, 4):
            pool = pools[n]
            for _ in range(30):
                t = (rng.choice(pool) + rng.choice(pool)) / 2
                assert np.linalg.norm(t, 2) <= 1 + TOL
                base = np.asarray(
                    [[rng.choice((-1, 1)) for _ in range(p)] for _ in range(n)],
                    dtype=float,
                )
                v = active_matrix(base, support)
                g, ray, defect = state(t, v)
                assert_psd(g + ray, "G+R")
                assert_psd(g - ray, "G-R")
                assert_psd(defect, "D")

                if np.linalg.eigvalsh(ray)[0] < -1e-7:
                    indefinite_cases += 1

                delta = np.linalg.norm(defect, 2)
                m = 0.7
                for eps in product((-1, 1), repeat=p):
                    c = endpoint_coeff(eps, support)
                    x = v @ c
                    expected = np.array([
                        selector_table(p)[tuple(eps[i] * int(base[j, i]) for i in range(p))]
                        for j in range(n)
                    ])
                    assert np.max(np.abs(x - expected)) < TOL
                    assert abs(c @ c - 1) < TOL
                    assert abs(c @ g @ c - 1) < TOL
                    actual, roof = exact_response(t, base, eps, m)
                    query_loss = n * float(c @ defect @ c) / 2
                    assert actual <= roof + TOL
                    assert roof - actual <= query_loss + TOL
                    assert query_loss <= n * delta / 2 + TOL
                    response_checks += 1

                # Test the indefinite contraction-kernel Schur lemma.
                for _ in range(12):
                    x = np_rng.normal(size=ray.shape)
                    x = (x + x.T) / 2
                    x /= max(np.linalg.norm(x, 2), 1e-15)
                    assert np.linalg.norm(ray * x, 2) <= 1 + 2e-8
                    schur_checks += 1
                instances.append((p, support, t, base, v, g, ray, defect))

    # Tensor only matching arity/support labels; use small factors.
    by_p = {}
    for item in instances:
        by_p.setdefault(item[0], []).append(item)
    for p, items in by_p.items():
        for _ in range(80):
            one = rng.choice(items)
            two = rng.choice(items)
            _, support, t1, b1, v1, g1, r1, d1 = one
            _, _, t2, b2, v2, g2, r2, d2 = two
            # Corresponding tensor base ports and active products.
            base12 = np.column_stack([
                np.kron(b1[:, i], b2[:, i]) for i in range(p)
            ])
            v12 = active_matrix(base12, support)
            t12 = np.kron(t1, t2)
            g12, r12, d12 = state(t12, v12)
            assert np.max(np.abs(g12 - g1 * g2)) < TOL
            assert np.max(np.abs(r12 - r1 * r2)) < TOL
            assert np.max(np.abs(d12 - (d1 * g2 + r1 * d2))) < TOL
            assert np.max(np.abs(d12 - (d1 * g2 + g1 * d2 - d1 * d2))) < TOL
            assert_psd(d12, "tensor D")
            rhs = np.linalg.norm(d1, 2) + np.linalg.norm(d2, 2)
            assert np.linalg.norm(d12, 2) <= rhs + 3e-8
            tensor_checks += 1

    # The one-port family gives exact response loss and first-order tensor
    # sharpness.  Use a concrete Boolean pole rather than a scalar surrogate.
    n = 6
    w = np.asarray([1, -1, 1, 1, -1, -1], dtype=float)
    for d in (0.001, 0.1, 0.7):
        t = np.eye(n) - d * np.outer(w, w) / n
        base = w[:, None]
        actual, roof = exact_response(t, base, (1,), max(d, 0.2))
        assert abs((roof - actual) - d * n / 2) < TOL
        _, _, defect = state(t, base)
        assert abs(defect[0, 0] - d) < TOL
    for d1 in (1e-5, 0.03, 0.4):
        for d2 in (2e-5, 0.07, 0.6):
            d12 = d1 + d2 - d1 * d2
            assert d12 <= d1 + d2 + TOL
    sharp_ratio = (1e-5 + 2e-5 - 2e-10) / 3e-5
    assert sharp_ratio > 0.99999

    raw_affine_defects, relative_checks = check_growing_affine_coset_boundary()
    initial_table, output_table0, output_table1 = check_endpoint_table_noncongruence()

    assert indefinite_cases > 0
    print(f"response checks: {response_checks}")
    print(f"indefinite Rayleigh cases: {indefinite_cases}")
    print(f"Hermitian Schur checks: {schur_checks}")
    print(f"tensor checks: {tensor_checks}")
    print(f"first-order sharpness ratio: {sharp_ratio:.9f}")
    print("PC.3 raw defects:", ", ".join(f"{x:.3f}" for x in raw_affine_defects))
    print(f"Cartesian relative-defect checks: {relative_checks}")
    print(
        "endpoint-table noncongruence: "
        f"initial={initial_table[0]:.6f}, "
        f"outputs={output_table0[0]:.6f}/{output_table1[0]:.6f}"
    )
    print("robust product synchronization verification: PASS")


if __name__ == "__main__":
    main()
