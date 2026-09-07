"""Finite normalization checks; not a numerical proof of nondisplaceability."""
import itertools
import numpy as np


def spins(n):
    return np.array(list(itertools.product((-1, 1), repeat=n)), dtype=float)


def cap(a):
    x = spins(len(a))
    return np.max(np.abs(np.einsum("bi,ij,bj->b", x, a, x))) / 2


def main():
    rng = np.random.default_rng(20260907)
    h = np.array([[1., 1.], [1., -1.]])
    z = np.exp(1j * np.pi * np.array([1., -3.]) / 8)
    assert np.allclose(h @ z, np.sqrt(2) * z.conj())
    assert np.allclose(z @ h @ z, 2 * np.sqrt(2))
    x = spins(2)
    assert np.max(np.abs(np.einsum("bi,ij,bj->b", x, h, x))) == 2
    phase_checks = lift_checks = 0
    for n in range(2, 9):
        for _ in range(10):
            a = np.triu(rng.choice((-1, 1), size=(n, n)), 1)
            a = a + a.T
            q = cap(a)
            for _ in range(20):
                z = np.exp(2j * np.pi * rng.random(n))
                assert abs(z @ a @ z) <= 4 * q + 1e-9
                # Explicit four-phase mean representation: put mass
                # |Re mu| and |Im mu| at their signed axis vertices;
                # put the remaining mass equally at +1 and -1.
                mu = z / np.sqrt(2)
                probabilities = np.zeros((n, 4))
                for j, w in enumerate(mu):
                    probabilities[j, 0 if w.real >= 0 else 2] += abs(w.real)
                    probabilities[j, 1 if w.imag >= 0 else 3] += abs(w.imag)
                    rest = 1 - abs(w.real) - abs(w.imag)
                    assert rest >= -1e-12
                    probabilities[j, 0] += rest / 2
                    probabilities[j, 2] += rest / 2
                assert np.allclose(probabilities @ np.array([1, 1j, -1, -1j]), mu)
                assert np.allclose(mu @ a @ mu, (z @ a @ z) / 2)
                phase_checks += 1
            ell = np.block([[a, a + np.eye(n)], [a + np.eye(n), -a]])
            for _ in range(40):
                w = rng.choice(np.array([1, 1j, -1, -1j]), size=n)
                u, v = w.real, w.imag
                xy = np.r_[u + v, u - v]
                lhs = xy @ ell @ xy / 2
                value = w @ a @ w
                rhs = value.real + value.imag + np.sum(u*u-v*v)
                assert abs(lhs-rhs) < 1e-9
                lift_checks += 1
    print({"phase_and_rounding_checks": phase_checks,
           "exact_lift_identity_checks": lift_checks,
           "H2_flat_conjugate_witness": "PASS"})


if __name__ == "__main__":
    main()
