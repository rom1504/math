#!/usr/bin/env python3
"""Finite wind tunnel for the contracting-fibre/cocycle decomposition."""

from __future__ import annotations

import math
import random

import numpy as np
from scipy.optimize import linprog


def stationary(P: np.ndarray) -> np.ndarray:
    vals, vecs = np.linalg.eig(P.T)
    i = int(np.argmin(np.abs(vals - 1.0)))
    v = np.real(vecs[:, i])
    if np.sum(v) < 0:
        v = -v
    v /= np.sum(v)
    return v


def l2(v: np.ndarray, pi: np.ndarray) -> float:
    return float(math.sqrt(np.dot(pi, v * v)))


def centred_rho(P: np.ndarray, pi: np.ndarray) -> float:
    # Conjugate to ordinary Euclidean L2 and remove the constant projection.
    root = np.sqrt(pi)
    S = root[:, None] * P / root[None, :]
    Pi = np.outer(root, root)
    return float(np.linalg.svd(S - Pi, compute_uv=False)[0])


def random_kernel(n: int, rng: random.Random) -> np.ndarray:
    A = np.array([[0.2 + rng.random() for _ in range(n)] for _ in range(n)])
    return A / A.sum(axis=1, keepdims=True)


def rectangular_rho(P: np.ndarray, pi_s: np.ndarray, pi_t: np.ndarray) -> float:
    rs = np.sqrt(pi_s)
    rt = np.sqrt(pi_t)
    S = rs[:, None] * P / rt[None, :]
    centred = S - np.outer(rs, rt)
    return float(np.linalg.svd(centred, compute_uv=False)[0])


def check_response_decomposition() -> int:
    rng = random.Random(20260817)
    checks = 0
    for _ in range(200):
        n = rng.choice([2, 3, 4])
        # Use one common kernel so its stationary law is transported exactly.
        P = random_kernel(n, rng)
        pi = stationary(P)
        rho = centred_rho(P, pi)
        assert rho < 1 + 1e-10
        rewards = [
            np.array([rng.uniform(-2, 2) for _ in range(n)])
            for _ in range(3)
        ]
        terminal = np.array([rng.uniform(-2, 2) for _ in range(n)])
        word = [rng.randrange(3) for _ in range(rng.randrange(1, 10))]

        value = terminal.copy()
        for e in reversed(word):
            value = rewards[e] + P @ value

        means = [float(np.dot(pi, a)) for a in rewards]
        bars = [a - m for a, m in zip(rewards, means)]
        ubar = float(np.dot(pi, terminal))
        v = terminal - ubar
        B = max(l2(b, pi) for b in bars)
        R = l2(v, pi)
        t = len(word)
        scalar = ubar + sum(means[e] for e in word)
        bound = rho**t * R + B * (1 - rho**t) / (1 - rho)
        centred = value - scalar
        assert abs(float(np.dot(pi, centred))) < 2e-9
        assert abs(l2(value, pi) ** 2 - scalar**2 - l2(centred, pi) ** 2) < 2e-8
        assert l2(centred, pi) <= bound + 2e-9
        checks += 1
    return checks


def check_rectangular_transport() -> int:
    rng = random.Random(27182818)
    checks = 0
    for _ in range(300):
        dims = [rng.choice([2, 3, 4, 5]) for _ in range(rng.randrange(2, 7))]
        pi = np.array([0.2 + rng.random() for _ in range(dims[0])])
        pi /= pi.sum()
        laws = [pi]
        kernels = []
        rhos = []
        rewards = []
        for ns, nt in zip(dims[:-1], dims[1:]):
            P = np.array([[0.2 + rng.random() for _ in range(nt)] for _ in range(ns)])
            P /= P.sum(axis=1, keepdims=True)
            pit = laws[-1] @ P
            kernels.append(P)
            laws.append(pit)
            rhos.append(rectangular_rho(P, laws[-2], pit))
            rewards.append(np.array([rng.uniform(-2, 2) for _ in range(ns)]))
        terminal = np.array([rng.uniform(-2, 2) for _ in range(dims[-1])])
        value = terminal.copy()
        for P, a in reversed(list(zip(kernels, rewards))):
            value = a + P @ value
        means = [float(np.dot(laws[i], a)) for i, a in enumerate(rewards)]
        bars = [a - means[i] for i, a in enumerate(rewards)]
        ubar = float(np.dot(laws[-1], terminal))
        v = terminal - ubar
        B = max(l2(b, laws[i]) for i, b in enumerate(bars))
        R = l2(v, laws[-1])
        rho = max(rhos)
        t = len(kernels)
        scalar = ubar + sum(means)
        bound = rho**t * R + B * (1 - rho**t) / (1 - rho)
        assert abs(float(np.dot(laws[0], value - scalar))) < 3e-9
        assert l2(value - scalar, laws[0]) <= bound + 3e-9
        checks += 1
    return checks


def simple_cycles(vertices: range, edges: list[tuple[int, int, float]]):
    adj: dict[int, list[tuple[int, float]]] = {v: [] for v in vertices}
    for a, b, w in edges:
        adj[a].append((b, w))
    out = []
    # Canonical enumeration of simple cycles by requiring the start be least.
    for s in vertices:
        stack = [(s, [s], [])]
        while stack:
            x, path, ws = stack.pop()
            for y, w in adj[x]:
                if y == s:
                    out.append(tuple(ws + [w]))
                elif y not in path and y >= s:
                    stack.append((y, path + [y], ws + [w]))
    return out


def check_cycle_law() -> int:
    rng = random.Random(17082026)
    checks = 0
    for _ in range(200):
        n = rng.choice([2, 3, 4])
        edges = []
        # Complete digraph guarantees strong connectivity and arbitrary paths.
        for a in range(n):
            for b in range(n):
                edges.append((a, b, rng.uniform(-2, 2)))
        cycles = simple_cycles(range(n), edges)
        chi = max(abs(sum(c) / len(c)) for c in cycles)
        M = max(abs(w) for _, _, w in edges)
        table = {(a, b): w for a, b, w in edges}
        max_end = np.zeros(n)
        min_end = np.zeros(n)
        for t in range(1, 25):
            max_end = np.array([
                max(max_end[a] + table[a, b] for a in range(n))
                for b in range(n)
            ])
            min_end = np.array([
                min(min_end[a] + table[a, b] for a in range(n))
                for b in range(n)
            ])
            worst = max(float(np.max(max_end)), float(-np.min(min_end)))
            assert worst <= t * chi + (n - 1) * M + 1e-9
            checks += 1

    delta = 1.7
    d = np.array([[-delta / 4, 3 * delta / 4],
                  [-delta / 4, -delta / 4]])
    edges = [(i, j, float(d[i, j])) for i in range(2) for j in range(2)]
    cycles = simple_cycles(range(2), edges)
    assert abs(max(abs(sum(c) / len(c)) for c in cycles) - delta / 4) < 1e-12
    checks += 1
    return checks


def check_stationary_flow_dual() -> int:
    rng = random.Random(31415926)
    checks = 0
    for _ in range(100):
        n = rng.choice([2, 3, 4, 5])
        edges = [(a, b, rng.uniform(-3, 3)) for a in range(n) for b in range(n)]
        cycles = simple_cycles(range(n), edges)
        chi = max(abs(sum(c) / len(c)) for c in cycles)
        Aeq = np.zeros((n + 1, len(edges)))
        beq = np.zeros(n + 1)
        beq[-1] = 1
        for j, (a, b, _) in enumerate(edges):
            Aeq[a, j] += 1
            Aeq[b, j] -= 1
            Aeq[-1, j] = 1
        weights = np.array([w for _, _, w in edges])
        hi = -linprog(-weights, A_eq=Aeq, b_eq=beq, bounds=(0, None), method="highs").fun
        lo = linprog(weights, A_eq=Aeq, b_eq=beq, bounds=(0, None), method="highs").fun
        assert abs(max(abs(hi), abs(lo)) - chi) < 2e-8
        checks += 1
    return checks


def check_potentials_and_diamond() -> int:
    rng = random.Random(16180339)
    checks = 0
    for n in range(2, 8):
        psi = [rng.uniform(-2, 2) for _ in range(n)]
        edges = [(a, b, psi[b] - psi[a]) for a in range(n) for b in range(n)]
        cycles = simple_cycles(range(n), edges)
        assert max(abs(sum(c)) for c in cycles) < 2e-10
        checks += 1
    # Acyclic diamond: unequal coterminal path sums, no global potential,
    # but every path has bounded length and hence no pumpable response.
    diamond = [(0, 1, 0.0), (0, 2, 0.0), (1, 3, 0.0), (2, 3, 1.0)]
    A = []
    b = []
    for s, t, w in diamond:
        row = np.zeros(4)
        row[t], row[s] = 1, -1
        A.append(row)
        b.append(w)
    sol = np.linalg.lstsq(np.array(A), np.array(b), rcond=None)[0]
    assert np.linalg.norm(np.array(A) @ sol - np.array(b)) > 0.1
    assert not simple_cycles(range(4), diamond)
    checks += 1
    return checks


def check_variance_tax() -> int:
    rng = random.Random(8172026)
    checks = 0
    for _ in range(2000):
        n = rng.choice([2, 3, 4, 5])
        P = random_kernel(n, rng)
        pi = stationary(P)
        rho = centred_rho(P, pi)
        f = np.array([rng.uniform(-3, 3) for _ in range(n)])
        noise = np.array([rng.uniform(-0.3, 0.3) for _ in range(n)])
        g = f + noise
        omega = float(np.max(np.abs(noise)))
        epsilon = max(0.0, float(np.max(f - P @ f)))
        sigma = l2(g - np.dot(pi, g), pi)
        denom = float(np.max(g) - np.min(g) + 2 * omega)
        rhs = 0.0 if denom == 0 else (1 - rho) ** 2 * max(0, sigma - omega) ** 2 / denom
        assert epsilon + 2e-9 >= rhs
        checks += 1
    return checks


def check_sharp_channels() -> int:
    checks = 0
    for rho in (0.0, 0.2, 0.7, 0.95):
        P = np.array([[(1 + rho) / 2, (1 - rho) / 2],
                      [(1 - rho) / 2, (1 + rho) / 2]])
        v = np.array([1.0, -1.0])
        for t in range(1, 15):
            assert np.max(np.abs(np.linalg.matrix_power(P, t) @ v - rho**t * v)) < 1e-11
            total = np.zeros(2)
            power = np.eye(2)
            for _ in range(t):
                total += power @ v
                power = power @ P
            target = (1 - rho**t) / (1 - rho) * v
            assert np.max(np.abs(total - target)) < 1e-11
            checks += 2
    return checks


def check_nonlinear_secant_and_discount() -> int:
    checks = 0
    for rho in (0.0, 0.3, 0.8):
        P = np.array([[(1 + rho) / 2, (1 - rho) / 2],
                      [(1 - rho) / 2, (1 + rho) / 2]])
        b = np.array([0.7, -0.7])
        m = 0.4
        eta = m + b
        # F(x)=Px and Fhat(x)=Px-eta.  Their secant is exactly P and their
        # same-input defect is eta, independently of the nonlinear path.
        x = np.array([1.2, -0.2])
        y = np.array([-0.4, 0.8])
        z0 = x - y
        c0 = float(np.mean(z0))
        R = l2(z0 - c0, np.array([0.5, 0.5]))
        for t in range(1, 20):
            x = P @ x
            y = P @ y - eta
            centred = x - y - (c0 + t * m)
            bound = rho**t * R + 0.7 * (1 - rho**t) / (1 - rho)
            assert l2(centred, np.array([0.5, 0.5])) <= bound + 2e-10
            checks += 1

        for lam in (0.2, 0.6, 0.95):
            scalar = 0.0
            centred = np.zeros(2)
            for t in range(1, 30):
                scalar = 1.3 + lam * scalar
                centred = b + lam * (P @ centred)
                assert abs(scalar) <= 1.3 / (1 - lam) + 1e-10
                assert l2(centred, np.array([0.5, 0.5])) <= 0.7 / (1 - lam * rho) + 1e-10
                checks += 2
    return checks


def main() -> None:
    count = 0
    count += check_response_decomposition()
    count += check_rectangular_transport()
    count += check_cycle_law()
    count += check_stationary_flow_dual()
    count += check_potentials_and_diamond()
    count += check_variance_tax()
    count += check_sharp_channels()
    count += check_nonlinear_secant_and_discount()
    print(f"contracting-fibre/cocycle checks passed: {count}")


if __name__ == "__main__":
    main()
