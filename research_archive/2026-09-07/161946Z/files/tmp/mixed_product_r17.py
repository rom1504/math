"""Numerical/exact diagnostics for A9 x Paley-conference products.

Scratch only.  All displayed large-order norms are heuristic lower bounds from
coordinate ascent, not exact maxima.
"""

from __future__ import annotations

import itertools
import math
import numpy as np


A9 = np.array([
    [0,1,-1,1,1,1,1,-1,-1],
    [1,0,1,-1,-1,-1,1,-1,-1],
    [-1,1,0,1,1,1,1,1,-1],
    [1,-1,1,0,1,1,1,-1,1],
    [1,-1,1,1,0,1,-1,1,-1],
    [1,-1,1,1,1,0,-1,-1,-1],
    [1,1,1,1,-1,-1,0,1,1],
    [-1,-1,1,-1,1,-1,1,0,-1],
    [-1,-1,-1,1,-1,-1,1,-1,0],
], dtype=np.int64)


def paley_conference(q: int) -> np.ndarray:
    """Symmetric Paley conference matrix of order q+1, q prime = 1 mod 4."""
    assert q % 4 == 1
    residues = {a * a % q for a in range(1, q)}
    chi = lambda a: 1 if a % q in residues else -1
    C = np.zeros((q + 1, q + 1), dtype=np.int64)
    C[0, 1:] = 1
    C[1:, 0] = 1
    for i in range(q):
        for j in range(q):
            if i != j:
                C[i + 1, j + 1] = chi(i - j)
    assert np.array_equal(C, C.T)
    assert np.array_equal(C @ C, q * np.eye(q + 1, dtype=np.int64))
    return C


def mixed(A: np.ndarray, d: np.ndarray, C: np.ndarray, e=None) -> np.ndarray:
    if e is None:
        e = np.ones(len(C), dtype=np.int64)
    P = A + np.diag(d)
    H = C + np.diag(e)
    T = np.kron(H, P)
    np.fill_diagonal(T, 0)
    assert set(np.unique(T)) <= {-1, 0, 1}
    return T


def ascent_one(M: np.ndarray, x: np.ndarray) -> tuple[int, np.ndarray]:
    # Maximize x'Mx. A flip changes energy by -4 x_i (Mx)_i.
    h = M @ x
    E = int(x @ h)
    n = len(x)
    while True:
        gains = -4 * x * h
        i = int(np.argmax(gains))
        if gains[i] <= 0:
            break
        old = int(x[i])
        x[i] = -old
        E += int(gains[i])
        # diagonal is zero
        h += -2 * old * M[:, i]
    assert E == int(x @ M @ x)
    return E, x


def heuristic_q(M: np.ndarray, restarts=80, seed=1) -> tuple[int, str, np.ndarray]:
    rng = np.random.default_rng(seed)
    n = len(M)
    best = -10**30
    kind = ""
    bestx = None
    for sign, label in ((1, "+"), (-1, "-")):
        N = sign * M
        for _ in range(restarts):
            x = rng.choice(np.array([-1, 1], dtype=np.int64), n)
            E, _ = ascent_one(N, x)
            if E > best:
                best = E
                kind = label
                bestx = x.copy()
    return best, kind, bestx


def q_exact(A: np.ndarray) -> tuple[int, int]:
    n = len(A)
    vals = []
    for bits in itertools.product((-1, 1), repeat=n - 1):
        x = np.array((1,) + bits, dtype=np.int64)
        vals.append(int(x @ A @ x))
    return min(vals), max(vals)


def main():
    print("A9 exact", q_exact(A9), "eig", np.linalg.eigvalsh(A9))
    choices = {
        "+all": np.ones(9, dtype=np.int64),
        "-all": -np.ones(9, dtype=np.int64),
        "trace1": np.array([1] * 5 + [-1] * 4, dtype=np.int64),
        "trace-1": np.array([1] * 4 + [-1] * 5, dtype=np.int64),
        "alt1": np.array([1,-1,1,-1,1,-1,1,-1,1], dtype=np.int64),
    }
    for q in (197,):
        C = paley_conference(q)
        m = q + 1
        print("\nouter", m)
        for name, d in choices.items():
            T = mixed(A9, d, C)
            val, kind, bestx = heuristic_q(T, restarts=8, seed=1000 + q + int(d.sum()))
            norm = val / (len(T) * math.sqrt(len(T) - 1))
            print(name, "delta", int(d.sum()), "lb", val, kind, "ratio", f"{norm:.6f}")
            if q == 197 and name == "trace1":
                X = bestx.reshape(m, 9)
                pats, cnt = np.unique(X, axis=0, return_counts=True)
                order = np.argsort(-cnt)
                print("patterns", [(int(cnt[i]), ''.join('+' if z > 0 else '-' for z in pats[i])) for i in order])


if __name__ == "__main__":
    main()
