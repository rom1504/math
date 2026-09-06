import itertools
import numpy as np


def build(m, d):
    half = m // 2
    s = np.r_[np.ones(half, dtype=int), -np.ones(half, dtype=int)]
    t = s.copy()
    G = np.zeros((m, m), dtype=int)
    for i in range(half):
        for a in range(d):
            j = half + (i + a) % half
            G[i, j] = G[j, i] = 1
    H = G.copy()
    K = np.zeros((m, m), dtype=int)
    for i in range(m):
        for a in range(d):
            K[i, (i + a) % m] = 1
    B = np.ones((m, m), dtype=int) - np.eye(m, dtype=int) - 2 * G
    D = -np.ones((m, m), dtype=int) + np.eye(m, dtype=int) + 2 * H
    C = np.ones((m, m), dtype=int) - 2 * K
    A = np.block([[B, s[:, None], C],
                  [s[None, :], np.zeros((1, 1), dtype=int), t[None, :]],
                  [C.T, t[:, None], D]])
    return A, s, t, G, H, K


def energy(A, x):
    return int(x @ A @ x // 2)


for m, d in [(4, 1), (6, 1), (8, 2)]:
    A, s, t, G, H, K = build(m, d)
    delta = m - 2 * d
    M0 = m * delta
    e = int(G.sum() // 2)
    plane = [
        (1, np.r_[np.ones(m, int), 1, np.ones(m, int)]),
        (1, np.r_[np.ones(m, int), -1, np.ones(m, int)]),
        (-1, np.r_[np.ones(m, int), 1, -np.ones(m, int)]),
        (-1, np.r_[np.ones(m, int), -1, -np.ones(m, int)]),
    ]
    scores = [rho * energy(A, x) for rho, x in plane]
    fields = [rho * x * (A @ x) for rho, x in plane]
    assert scores == [M0] * 4
    assert all(np.all(r >= 0) for r in fields)
    assert np.all(np.min(np.stack(fields), axis=0) == 0)
    vals = []
    n = 2 * m + 1
    for bits in itertools.product((-1, 1), repeat=n - 1):
        x = np.r_[1, bits]
        vals.append(energy(A, x))
    bound = (m * m + m - M0) // 2
    print({"m": m, "d": d, "M0": M0, "e": e,
           "min": min(vals), "max": max(vals), "bound": bound,
           "fields": [sorted(set(map(int, r))) for r in fields]})
