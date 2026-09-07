"""Independent finite replay of optimized width's log-cosh interpolation.

Only numpy and the archived order-14 witness are inputs.  No asymptotic
conclusion is inferred from the finite replay.
"""
import itertools
import json
from pathlib import Path
import numpy as np


def signs(k):
    return np.array(list(itertools.product((-1., 1.), repeat=k)))


def problem(n):
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    x = np.column_stack((np.ones(2 ** (n - 1)), signs(n - 1)))
    v = np.array([x[:, i] * x[:, j] for i, j in edges]).T
    free = [e for e, (i, _) in enumerate(edges) if i > 0]
    a = np.ones((2 ** len(free), len(edges)))
    a[:, free] = signs(len(free))
    return edges, v, a


def branch_statistics(v, a, lam):
    e = (a * lam) @ v.T
    logs, correlations = [], []
    for s in (1., -1.):
        z = s * e
        peak = z.max(axis=1)
        weights = np.exp(z - peak[:, None])
        den = weights.sum(axis=1)
        logs.append(peak + np.log(den))
        correlations.append(weights @ v / den[:, None])
    return .5 * (logs[0] + logs[1]), correlations


def replay():
    cases = edge_checks = 0
    for n in range(3, 7):
        edges, v, a = problem(n)
        for m in range(1, n // 2 + 1):
            z = np.array([np.sqrt((n - m) / m)] * m
                         + [-np.sqrt(m / (n - m))] * (n - m))
            D = np.array([z[i] * z[j] for i, j in edges])
            assert abs(z.sum()) < 1e-12
            for beta in (.2, .8, 1.7):
                tau = np.log(np.cosh(2 * beta / np.sqrt(n)))
                for u in (.03, .27, .61, .94):
                    lam = .5 * np.arccosh(np.exp(tau * (1 + u * D)))
                    f, cs = branch_statistics(v, a, lam)
                    best = int(f.argmin())
                    A = a[best]
                    cp, cm = cs[0][best], cs[1][best]
                    p = cp * cm
                    d = .5 * A * (cp - cm)
                    slack = 1 - p - 2 * d / np.tanh(2 * lam)
                    assert slack.min() > -2e-10
                    assert np.max(d - np.tanh(lam)) < 2e-10
                    P = np.eye(n)
                    for e, (i, j) in enumerate(edges):
                        P[i, j] = P[j, i] = p[e]
                    assert np.linalg.eigvalsh(P).min() > -2e-10
                    lamprime = tau * D / (2 * np.tanh(2 * lam))
                    derivative = d @ lamprime
                    rhs = -tau * (z @ P @ z) / 8 - tau * (D @ slack) / 4
                    assert abs(derivative - rhs) < 2e-10
                    t = np.tanh(lam)
                    c0p = (cp - A * t) / (1 - A * t * cp)
                    c0m = (cm + A * t) / (1 + A * t * cm)
                    b = abs(c0p - c0m)
                    assert np.max(abs(A * (c0p - c0m) + b)) < 2e-9
                    q = 1 - t * b - t * t * c0p * c0m
                    formula = b * (1 - t * t) ** 2 / (2 * t * q)
                    assert np.max(abs(formula - slack)) < 2e-9
                    for e in range(len(edges)):
                        flipped = A.copy()
                        flipped[e] *= -1
                        ff, _ = branch_statistics(v, flipped[None, :], lam)
                        delta = ff[0] - f[best]
                        exact = np.expm1(2 * delta) / np.sinh(2 * lam[e]) ** 2
                        assert delta > -2e-10
                        assert abs(exact - slack[e]) < 2e-9
                        edge_checks += 1
                    eps = 1e-5
                    values = []
                    for uu in (u - eps, u + eps):
                        ll = .5 * np.arccosh(np.exp(tau * (1 + uu * D)))
                        values.append(branch_statistics(v, a, ll)[0].min())
                    numeric = (values[1] - values[0]) / (2 * eps)
                    assert abs(numeric - derivative) < 2e-7
                    if n == 4 and m == 2:
                        xx, yy = lam.max(), lam.min()
                        tt, vv = np.tanh(xx), np.tanh(yy)
                        # This replay sums only projective spins.  Add log 2
                        # to recover the full-spin formula in the audit.
                        closed = (4 * np.log(2) + 2 * np.log(np.cosh(xx))
                                  + 4 * np.log(np.cosh(yy))
                                  + np.log(1 + vv ** 4 - 2 * tt * tt * vv * vv))
                        assert abs(f[best] + np.log(2) - closed) < 2e-10
                    cases += 1
    path = Path(__file__).resolve().parents[1] / 'computations/results/conference_completion_m13.json'
    record = json.loads(path.read_text())
    matrix = np.array(record['conference_matrix'], dtype=int)
    x = np.column_stack((np.ones(2 ** 13, dtype=int), signs(13).astype(int)))
    energies = np.einsum('bi,ij,bj->b', x, matrix, x) // 2
    assert matrix.shape == (14, 14)
    assert np.all((matrix * matrix).sum(axis=1) == 13)
    assert abs(energies).max() == 21
    print(json.dumps({'path_cases': cases, 'exact_edge_flip_checks': edge_checks,
                      'order14_projective_spins': len(x), 'order14_cap': 21,
                      'status': 'PASS'}, sort_keys=True))


if __name__ == '__main__':
    replay()
