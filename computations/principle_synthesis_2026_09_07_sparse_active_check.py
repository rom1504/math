"""Exact checks of the balanced-column sparse-active inequality."""
import json
from fractions import Fraction

import numpy as np


def walsh(n):
    return np.array([[(-1) ** bin(i & j).count("1") for j in range(n)]
                     for i in range(n)], dtype=np.int64)


def main():
    rng = np.random.default_rng(202609072110)
    results = []
    for m in (8, 16, 32):
        block = 4
        q = 3 * m // 4
        had = walsh(m)
        frames = []
        for i in range(m):
            selected = []
            for group in range(m // block):
                selected.extend(group * block + rng.choice(block, block - 1, replace=False))
            frame = had[np.array(selected)][:, np.arange(m) ^ i].copy()
            frame[:, i] = 0
            for j in range(m):
                excess = int(frame[:, j].sum())
                if excess:
                    locations = np.flatnonzero(frame[:, j] == np.sign(excess))
                    flips = rng.choice(locations, abs(excess) // 2, replace=False)
                    frame[flips, j] *= -1
            assert np.all(frame.sum(axis=0) == 0)
            frames.append(frame)
        gram = [frame.T @ frame for frame in frames]
        mu = max(int(np.max(np.abs(g - np.diag(np.diag(g))))) for g in gram)
        signs = rng.choice((-1, 1), (m, m))
        signs = np.triu(signs, 1)
        signs += signs.T
        checked = 0
        for ell in (1, 2, 3, 4):
            for unused in range(100):
                active = rng.choice(m, ell, replace=False)
                words = np.ones((m, q), dtype=np.int64)
                for i in active:
                    while True:
                        word = rng.choice((-1, 1), q)
                        if abs(int(word.sum())) < q:
                            break
                    words[i] = word
                ports = np.array([frames[i].T @ words[i] for i in range(m)])
                energy = sum(int(signs[i, j] * ports[i, j] * ports[j, i])
                             for i in range(m) for j in range(i + 1, m))
                variance = sum(Fraction(q * q - int(word.sum()) ** 2, q) for word in words)
                if ell == 1:
                    assert energy == 0
                else:
                    assert 2 * abs(energy) <= (q + (ell - 2) * mu) * variance
                checked += 1
        sharp = np.ones((m, q), dtype=np.int64)
        sharp[0] = frames[0][:, 1]
        sharp[1] = signs[0, 1] * frames[1][:, 0]
        ports = np.array([frames[i].T @ sharp[i] for i in range(m)])
        sharp_energy = sum(int(signs[i, j] * ports[i, j] * ports[j, i])
                           for i in range(m) for j in range(i + 1, m))
        sharp_variance = sum(Fraction(q * q - int(word.sum()) ** 2, q) for word in sharp)
        assert sharp_energy == q * q
        assert sharp_variance == 2 * q
        results.append({"m": m, "q": q, "mu": mu, "inequalities_checked": checked,
                        "sharp_energy": sharp_energy, "sharp_variance": str(sharp_variance),
                        "status": "PASS"})
    print(json.dumps({"status": "PASS", "results": results}, indent=2))


if __name__ == "__main__":
    main()
