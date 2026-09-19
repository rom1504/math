import itertools
import random

import numpy as np


def extrema(A):
    vals = []
    for bits in itertools.product((-1, 1), repeat=len(A)):
        x = np.array(bits, dtype=np.int64)
        vals.append((int(x @ A @ x), x))
    P = max(e for e, _ in vals)
    N = -min(e for e, _ in vals)
    return P, N, max(P, N), vals


def make_path(A, y, rho, rng):
    R = list(range(len(A)))
    nodes = []
    while len(R) > 1:
        C0 = A[np.ix_(R, R)]
        P, N, q, vals = extrema(C0)
        orientations = []
        if P == q:
            orientations.append(1)
        if N == q:
            orientations.append(-1)
        # Deliberately encourage resets whenever possible.
        tau = -rho if -rho in orientations and rng.random() < 0.75 else rng.choice(orientations)
        grounds = [x for e, x in vals if tau * e == q]
        rng.shuffle(grounds)
        choice = None
        for x in grounds:
            for z in (x, -x):
                Di = np.flatnonzero(z != y)
                Ei = np.flatnonzero(z == y)
                if len(Di) and len(Ei):
                    choice = z.copy(), Di, Ei
                    break
            if choice:
                break
        if choice is None:
            break
        x, Di, Ei = choice
        B = C0[np.ix_(Di, Ei)]
        D = C0[np.ix_(Di, Di)]
        E = C0[np.ix_(Ei, Ei)]
        HD = int(tau * (x[Di] @ D @ x[Di]))
        HE = int(tau * (x[Ei] @ E @ x[Ei]))
        PD, ND, qD, _ = extrema(D)
        PE, NE, qE, _ = extrema(E)
        inherited_oriented = int(rho * (y @ C0 @ y))
        g = q - inherited_oriented
        sigma = rho * tau
        cut_a = q - int(tau * (y @ C0 @ y))
        assert (sigma == 1 and cut_a == g) or (sigma == -1 and cut_a == 2 * q - g)
        c = int(tau * (x[Di] @ B @ x[Ei]))
        assert cut_a == 4 * c
        h = None
        if sigma == 1:
            gammaD = qD - HD
            gammaE = qE - HE
            LF = int(np.abs(B @ x[Ei]).sum())
            zeta = qD + qE - q
            f = max(2 * LF - gammaE, 0)
            h = max(g // 2 - zeta - f, 0)
            assert h <= 2 * c
        nodes.append({
            "R": R.copy(), "tau": tau, "rho": rho, "sigma": sigma,
            "q": q, "P": P, "N": N, "g": g, "a": cut_a,
            "c": c, "h": h,
        })
        R = [R[i] for i in Ei]
        y = x[Ei].copy()
        rho = tau
    return nodes


def audit(A, y, rho, rng):
    nodes = make_path(A, y, rho, rng)
    if not nodes:
        return 0
    P0, N0, _, _ = extrema(A)
    lhs = sum(z["h"] for z in nodes if z["sigma"] == 1)
    reset_a = sum(z["a"] for z in nodes if z["sigma"] == -1)
    assert lhs <= P0 + N0 + reset_a, (lhs, P0, N0, reset_a, nodes)
    return len(nodes)


def main():
    rng = random.Random(90210)
    paths = steps = 0
    for n in range(3, 10):
        for _ in range(500):
            A = np.zeros((n, n), dtype=np.int64)
            for i in range(n):
                for j in range(i + 1, n):
                    A[i, j] = A[j, i] = rng.choice((-1, 1))
            y = np.array([rng.choice((-1, 1)) for _ in range(n)], dtype=np.int64)
            rho = rng.choice((-1, 1))
            steps += audit(A, y, rho, rng)
            paths += 1
    print(f"audited {paths} paths and {steps} compatible steps with resets")


if __name__ == "__main__":
    main()
