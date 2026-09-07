"""Finite algebra checks for the actual-sign integrated payment obstruction."""
import json
import math
import numpy as np


def hadamard(n):
    h = np.ones((1, 1), dtype=np.int64)
    while len(h) < n:
        h = np.block([[h, h], [h, -h]])
    return h


def check(k):
    ell = 64 * k
    n = k * ell
    h = hadamard(ell // 2)
    b = np.block([[h, -h], [-h, h]])
    s = hadamard(k)
    s -= np.diag(np.diag(s))
    assert np.all(np.abs(b) == 1)
    assert np.all(b.sum(axis=0) == 0)
    assert np.all(b == b.T)
    assert np.all(s == s.T) and np.all(np.diag(s) == 0)
    assert np.all(np.abs(s[np.triu_indices(k, 1)]) == 1)
    off_norm = np.linalg.norm(s, 2) * math.sqrt(2 * ell)
    assert off_norm < ell
    rng = np.random.default_rng(k)
    coherent_max = (n * ell - n) // 2
    for _ in range(24):
        x = rng.choice([-1, 1], size=(k, ell))
        means = x.sum(axis=1)
        delta = (k * ell**2 - means @ means) / 2
        energy = (means @ means - n) / 2
        energy += np.sum(x * (s @ x @ b)) / 2
        assert energy <= coherent_max - delta + off_norm * delta / ell + 1e-7
        assert energy >= -n / 2 - n * off_norm / 2 - 1e-7
    coefficient_bound = (
        math.log(2) / 2
        + (math.sqrt(2) * (1 + 1 / math.sqrt(k))
           - 8 * (math.sqrt(2) - 1) + math.sqrt(2 / n)) / 4
    )
    return dict(k=k, ell=ell, n=n, off_norm_over_sqrt_n=off_norm / math.sqrt(n),
                payment_gap_upper_per_vertex=coefficient_bound, checks="PASS")


def flat_eigenvector(n, negative=False):
    assert n >= 4
    v = np.array([-1, 1, 1, 1] if negative else [1, 1, 1, -1])
    while len(v) < n:
        v = np.kron(v, np.array([1, 1, 1, -1]))
    assert len(v) == n
    assert np.array_equal(hadamard(n) @ v, (-1 if negative else 1) * int(math.sqrt(n)) * v)
    return v


def sharper_check(k):
    ell, n = 4 * k, 4 * k * k
    inner = flat_eigenvector(ell // 2, negative=True)
    h = hadamard(ell // 2)
    b = np.block([[h, -h], [-h, h]])
    z = np.concatenate([inner, -inner])
    assert z.sum() == 0
    assert np.array_equal(b @ z, -int(math.sqrt(2 * ell)) * z)
    d = n // 2
    for sign in [-1, 1]:
        s = flat_eigenvector(k // 2, negative=sign < 0)
        outer = sign * hadamard(k // 2)
        outer -= np.diag(np.diag(outer))
        energy = -d / 2 + (s @ outer @ s) * (z @ b @ z) / 2
        claimed = -d / 2 - math.sqrt(2) * d**1.5 / 2
        assert abs(energy - claimed) < 1e-5
    bound = math.log(2) - 2 * (math.sqrt(2) - 1) + math.sqrt(2 / k)
    return dict(k=k, ell=ell, n=n, sharper_gap_upper_per_vertex=bound, checks="PASS")


if __name__ == "__main__":
    for k in [2, 4, 8]:
        print(json.dumps(check(k)), flush=True)
    limiting = math.log(2) / 2 + (math.sqrt(2) - 8 * (math.sqrt(2) - 1)) / 4
    assert limiting < -0.1283
    print(json.dumps(dict(limiting_upper_coefficient=limiting, checks="PASS")))
    for k in [8, 32, 128]:
        print(json.dumps(sharper_check(k)), flush=True)
