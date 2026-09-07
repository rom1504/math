"""Exploratory lower witnesses for fixed 7^e Gauss phase limits."""
import argparse
import numpy as np
from scipy.optimize import minimize_scalar


def phase_operator(e, theta):
    v = 7**e
    lam = np.zeros(v, dtype=complex)
    qr = {1, 2, 4}
    for j in range(1, v):
        jj, valuation = j, 0
        while jj % 7 == 0:
            jj //= 7
            valuation += 1
        eps = 1 if jj % 7 in qr else -1
        lam[j] = np.exp(1j * eps * 7**valuation * theta)
    # The opposite theta convention gives the same optimized family.
    kernel = np.fft.ifft(lam).real
    indices = np.arange(v)
    return kernel[(indices[:, None] + indices[None, :]) % v]


def objective(U, B, f):
    return np.abs(U @ f @ B).sum() / (2 * len(U))


def alternate(U, B, f, iterations=100):
    previous = -1
    for _ in range(iterations):
        g = np.where(U @ f @ B >= 0, 1., -1.)
        f = np.where(U @ g @ B >= 0, 1., -1.)
        value = objective(U, B, f)
        if value <= previous + 1e-12:
            break
        previous = value
    return value, f


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--e', type=int, default=2)
    parser.add_argument('--restarts', type=int, default=100)
    parser.add_argument('--angles', type=int, default=151)
    parser.add_argument('--seed', type=int, default=719)
    parser.add_argument('--ratio', type=float, default=2.)
    args = parser.parse_args()
    v = 7**args.e
    r = args.ratio
    B = np.array([[-1., r], [r, r*r]])
    rng = np.random.default_rng(args.seed)
    best = (-1, None, None)
    for theta in np.linspace(0, np.pi, args.angles):
        U = phase_operator(args.e, theta)
        assert np.max(np.abs(U@U - np.eye(v) + np.ones((v,v))/v)) < 1e-9
        local = (-1, None)
        seeds = [rng.choice([-1., 1.], (v,2)) for _ in range(args.restarts)]
        if best[1] is not None:
            seeds.append(best[1])
        for f in seeds:
            score, candidate = alternate(U, B, f)
            if score > local[0]:
                local = score, candidate
        if local[0] > best[0] + 1e-10:
            best = local[0], local[1], theta
            print('record', best[0], theta, flush=True)
    f, theta = best[1], best[2]
    for _ in range(30):
        fit = minimize_scalar(lambda th: -objective(phase_operator(args.e, th), B, f),
                              bounds=(theta-.08, theta+.08), method='bounded',
                              options={'xatol': 1e-14})
        theta = fit.x
        value, f = alternate(phase_operator(args.e, theta), B, f)
        if value <= best[0] + 1e-12:
            break
        best = value, f, theta
    print('BEST', best[0], best[2], flush=True)
    print(best[1].astype(int).tolist(), flush=True)
    print('means', best[1].mean(axis=0), flush=True)


if __name__ == '__main__':
    main()
