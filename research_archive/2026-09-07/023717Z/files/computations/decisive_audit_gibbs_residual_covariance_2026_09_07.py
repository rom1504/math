"""Exact-state diagnostics of the actual-optimizer covariance theorem.

All switching classes at orders 3..6 are enumerated at each beta. The
Gaussian orientation field is zero here. Matrix inequalities use floating
eigenvalues; the proof, not these diagnostics, establishes the theorem.
"""
from itertools import combinations, product
import math
import numpy as np


def states(n):
    x = np.asarray(list(product((-1, 1), repeat=n)), dtype=float)
    return np.concatenate((x, x)), np.concatenate((np.ones(len(x)), -np.ones(len(x))))


def law(a, beta, h=0.0):
    n = len(a)
    x, s = states(n)
    j = beta*a/math.sqrt(n)
    fields = x @ j
    energy = np.einsum("ij,ij->i", x, fields)/2
    logw = s*energy
    if h:
        nodes, quadw = np.polynomial.hermite.hermgauss(64)
        logits = logw[None, :]+h*math.sqrt(2)*nodes[:, None]*s[None, :]
        peaks = logits.max(axis=1)
        unnorm = np.exp(logits-peaks[:, None])
        totals = unnorm.sum(axis=1)
        quadw /= math.sqrt(math.pi)
        logz = float(quadw @ (peaks+np.log(totals)))
        weights = quadw @ (unnorm/totals[:, None])
    else:
        peak = logw.max()
        weights = np.exp(logw-peak)
        logz = float(peak+math.log(weights.sum()))
        weights /= weights.sum()
    return logz, x, s, j, fields, weights


def gauge_signings(n):
    pairs = list(combinations(range(1, n), 2))
    for signs in product((-1, 1), repeat=len(pairs)):
        a = np.ones((n, n))-np.eye(n)
        for (i, j), sign in zip(pairs, signs):
            a[i, j] = a[j, i] = sign
        yield a


def score_check(a, beta, optimizer, h=0.0):
    n = len(a)
    phi, x, s, j, fields, p = law(a, beta, h)
    local_v = 1/np.cosh(fields)**2
    eta = x-s[:, None]*np.tanh(fields)
    v = p @ local_v
    cross = eta.T @ (p[:, None]*x)
    assert np.max(np.abs(cross-np.diag(v))) < 2e-11
    cov = x.T @ (p[:, None]*x)
    scores = eta.T @ (p[:, None]*eta)
    leading = -j*((local_v.T*(p*s)) @ local_v)
    error = scores-np.diag(v)-leading
    assert np.max(np.abs(error)-2*j*j) < 2e-11
    d = float(np.linalg.norm(j, 2))
    upper = 1+d+2*np.sum(j*j, axis=1)
    assert np.linalg.eigvalsh(np.diag(upper)-scores).min() > -2e-11
    lower = v*v/upper
    assert np.linalg.eigvalsh(cov-np.diag(lower)).min() > -2e-11
    if optimizer:
        assert v.min() >= math.exp(-beta*beta)-2e-11
        for i in range(n):
            keep = [z for z in range(n) if z != i]
            # The inherited cavity uses the ORIGINAL edge normalization.
            sub = a[np.ix_(keep, keep)]
            child_phi, _, _, _, _, _ = law(sub, beta*math.sqrt((n-1)/n), h)
            expected_log_c = phi-math.log(2)-child_phi
            assert expected_log_c <= (n-1)*math.log(math.cosh(beta/math.sqrt(n)))+2e-10
    return cov, lower


def main():
    classes = optima = generic = 0
    saved = []
    for n in (3, 4, 5, 6):
        for beta in (0.3, 0.9, 1.8):
            rows = list(gauge_signings(n))
            values = [law(a, beta)[0] for a in rows]
            minimum = min(values)
            for a, value in zip(rows, values):
                classes += 1
                # Score bound requires no optimization at all.
                is_opt = abs(value-minimum) < 1e-11
                result = score_check(a, beta, is_opt)
                generic += 1
                if is_opt:
                    optima += 1
                    saved.append((n, result))
    quenched = 0
    for h in (0.7, 3.0):
        for beta in (0.9, 1.8):
            rows = list(gauge_signings(4))
            values = [law(a, beta, h)[0] for a in rows]
            minimum = min(values)
            for a, value in zip(rows, values):
                score_check(a, beta, abs(value-minimum) < 1e-11, h)
                quenched += 1
    # Global optimization is essential for the uniform row-variance floor.
    bad = np.ones((6, 6))-np.eye(6)
    _, _, _, _, fields, p = law(bad, 1.8)
    bad_v = p @ (1/np.cosh(fields)**2)
    assert bad_v.min() < math.exp(-1.8**2)
    bridge_checks = 0
    rng = np.random.RandomState(20260907)
    for n, (cov, lower) in saved[:30]:
        for m, (other, other_lower) in saved[-10:]:
            b = rng.choice((-1., 1.), size=(n, m))
            variance = np.trace(b @ other @ b.T @ cov)
            assert variance >= lower.sum()*other_lower.sum()-2e-10
            bridge_checks += 1
    print("PASS: %d switching classes, %d global minima, %d score tests, %d bridges; "
          "%d quenched checks; nonoptimizer row floor fails as expected"
          % (classes, optima, generic, bridge_checks, quenched))


if __name__ == "__main__":
    main()
