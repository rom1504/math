"""Exact noncommuting example and numerical matrix precision diagnostics."""

from fractions import Fraction as F
import numpy as np


def inverse2(a):
    det = a[0, 0] * a[1, 1] - a[0, 1] * a[1, 0]
    return np.array([[a[1, 1], -a[0, 1]], [-a[1, 0], a[0, 0]]], dtype=object) / det


def exact_case():
    l1 = np.array([[F(1), F(0)], [F(0), F(1, 3)]], dtype=object)
    l2 = np.array([[F(2, 3), F(1, 6)], [F(1, 6), F(2, 3)]], dtype=object)
    assert np.any(l1 @ l2 != l2 @ l1)
    p, d = (l1 + l2) / 2, (l1 - l2) / 2
    s = p - d @ inverse2(p) @ d
    assert np.all(s == 2 * inverse2(inverse2(l1) + inverse2(l2)))
    for x in range(-3, 4):
        ea = np.array([F(x, 2), F(2, 3)], dtype=object)
        eb = np.array([F(-1, 3), F(x, 5)], dtype=object)
        lhs = ((ea + eb) @ l1 @ (ea + eb) + (ea - eb) @ l2 @ (ea - eb)) / 2
        z = ea + inverse2(p) @ d @ eb
        assert lhs == z @ p @ z + eb @ s @ eb


def main():
    exact_case()
    rng = np.random.default_rng(2026090624)
    min_curvature = float("inf")
    max_harmonic_error = 0.0
    max_majorization_error = 0.0
    count = 0
    for k in (1, 2, 3, 4, 8):
        for _ in range(100):
            t = float(np.exp(rng.uniform(-2, 2)))
            def precision():
                q, _ = np.linalg.qr(rng.normal(size=(k, k)))
                return (q * (t * np.exp(rng.uniform(-5, 0, size=k)))) @ q.T
            l1, l2 = precision(), precision()
            p, d = (l1 + l2) / 2, (l1 - l2) / 2
            s = p - d @ np.linalg.solve(p, d)
            harmonic = 2 * np.linalg.inv(np.linalg.inv(l1) + np.linalg.inv(l2))
            error = np.linalg.norm(s-harmonic) / t
            max_harmonic_error = max(max_harmonic_error, error)
            assert error < 1e-12
            piv = np.r_[np.linalg.eigvalsh(p), np.linalg.eigvalsh(s)]
            eigen = np.r_[np.linalg.eigvalsh(l1), np.linalg.eigvalsh(l2)]
            assert piv.min() > 0 and piv.max() <= t * (1 + 1e-12)
            sums = np.cumsum(np.sort(np.log(piv))[::-1]) - np.cumsum(np.sort(np.log(eigen))[::-1])
            max_majorization_error = max(max_majorization_error, float(sums.max()))
            assert sums.max() < 1e-10 and abs(sums[-1]) < 1e-10
            c = lambda vals: float(np.log(vals * (2*t-vals) / t**2).sum() / 4)
            slack = c(piv) - c(eigen)
            min_curvature = min(min_curvature, slack)
            assert slack >= -1e-12
            count += 1
    print({"status": "PASS", "exact_noncommuting_square_tests": 7,
           "random_matrix_cases": count, "minimum_curvature_slack": min_curvature,
           "maximum_harmonic_error": max_harmonic_error,
           "maximum_log_majorization_error": max_majorization_error})


if __name__ == "__main__":
    main()
