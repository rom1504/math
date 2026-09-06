#!/usr/bin/env python3
"""Finite checks of the near-equator switching theorem.

Self-contained; no input files or generated outputs. Finite tests are
not certificates for the asymptotic statements in the companion proof.
"""
from itertools import product
import math
import numpy as np


def cube(n):
    return np.array(list(product((-1, 1), repeat=n)), dtype=np.int64)


def energies(a, spins):
    return np.einsum('bi,ij,bj->b', spins, a, spins) / 2


def small_cube_checks(rng):
    moments = 0
    balances = 0
    for n in range(3, 10):
        spins = cube(n)
        for _ in range(15):
            raw = rng.choice((-1, 1), size=(n, n))
            a = np.triu(raw, 1)
            a = a + a.T
            h = energies(a, spins)
            sigma2 = float(np.sum(a*a)/2)
            g4 = 3*sigma2*sigma2 + 3*float(np.trace(a @ a @ a @ a))
            assert float(np.mean(h**4)) <= g4 + 1e-7
            assert g4 <= 15*sigma2*sigma2 + 1e-7
            moments += 1
            x = spins[np.argmin(np.abs(h))]
            switched = a * x[:, None] * x[None, :]
            initial_sum = int(np.sum(switched)//2)
            edits = 0
            while abs(np.sum(switched)//2) > 1:
                sign = 1 if np.sum(switched) > 0 else -1
                row, col = next((i, j) for i in range(n) for j in range(i+1, n)
                                if switched[i, j] == sign)
                switched[row, col] *= -1
                switched[col, row] *= -1
                edits += 1
            assert abs(np.sum(switched)//2) <= 1
            assert 2*edits <= abs(initial_sum)
            assert np.max(np.abs(energies(switched, spins))) <= np.max(np.abs(h)) + 2*edits
            balances += 1
            weighted = rng.normal(size=(n, n))
            weighted = np.triu(weighted, 1)
            weighted += weighted.T.copy()
            wh = energies(weighted, spins)
            ws2 = float(np.sum(weighted*weighted)/2)
            wg4 = 3*ws2*ws2 + 3*float(np.trace(weighted @ weighted @ weighted @ weighted))
            assert float(np.mean(wh**4)) <= wg4 + 1e-7
            assert wg4 <= 15*ws2*ws2 + 1e-7
            moments += 1
    return moments, balances


def path_checks(rng):
    cases = 0
    for n in (24, 48, 96):
        for _ in range(12):
            a = np.triu(rng.choice((-1, 1), size=(n, n)), 1)
            a += a.T.copy()
            bound = math.sqrt(2*(n-1)*math.log(7200*n))
            for attempt in range(10000):
                x = rng.choice((-1, 1), size=n)
                y = rng.choice((-1, 1), size=n)
                hx = int(x @ a @ x // 2)
                hy = int(y @ a @ y // 2)
                if hx*hy >= 0:
                    continue
                z = x.copy()
                h = hx
                best = abs(h)
                maximum_field = 0
                for i in range(n):
                    field = int(a[i] @ z)
                    maximum_field = max(maximum_field, abs(field))
                    h += int(y[i]-z[i])*field
                    z[i] = y[i]
                    best = min(best, abs(h))
                assert h == hy
                if maximum_field <= bound:
                    assert best <= bound
                    cases += 1
                    break
            else:
                raise AssertionError('finite randomized diagnostic did not find a certified path')
    return cases


def main():
    rng = np.random.default_rng(202609062244)
    moments, balances = small_cube_checks(rng)
    print('fourth_moment_cases', moments)
    print('exact_balance_and_cap_cases', balances)
    print('certified_coordinate_paths', path_checks(rng))
    print('PASS: finite algebra and witnessed paths, not asymptotic proof')


if __name__ == '__main__':
    main()
