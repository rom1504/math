#!/usr/bin/env python3
"""Exact checks for the Wave 51 spectral-excess coarea audit."""

from fractions import Fraction
from math import comb


def falling(x: int, j: int) -> int:
    ans = 1
    for i in range(j):
        ans *= x - i
    return ans


def eig(n: int, m: int, ell: int, j: int) -> Fraction:
    if ell < j:
        return Fraction(0)
    return Fraction(
        falling(ell, j) * falling(n - m, j),
        falling(m, j) * falling(n - ell, j),
    )


def cap10_table() -> None:
    n, m, family_size = 10, 6, 5
    # Ordered intersection histogram of the two identical cap-ten families.
    hist = {2: 10, 4: 10, 6: 5}
    actual_degree = Fraction(1, 42)
    expected_gaps = {1: Fraction(1, 42), 2: Fraction(1, 75)}
    expected_bounds = {1: Fraction(1, 70), 2: Fraction(14, 2475)}
    rows = []
    for ell in range(1, m):
        denom = comb(m, ell) * comb(n - ell, m - ell)
        p = Fraction(
            sum(count * comb(overlap, ell)
                for overlap, count in hist.items()),
            family_size * denom,
        )
        lam1, lam2 = eig(n, m, ell, 1), eig(n, m, ell, 2)
        gap = p - lam2
        spectral_denom = 1 - lam1 + n * (lam1 - lam2)
        bound = max(gap, 0) / spectral_denom
        assert bound <= actual_degree
        if ell in expected_gaps:
            assert gap == expected_gaps[ell]
            assert bound == expected_bounds[ell]
        else:
            assert gap < 0
        rows.append((ell, lam1, lam2, p, gap, spectral_denom, bound))

    print("cap-10 spectral-excess table")
    print("ell | lambda1 | lambda2 | P | P-lambda2 | denominator | degree bound")
    for row in rows:
        print(" | ".join(map(str, row)))
    print(f"actual degree={actual_degree}")


def k0_endpoint_check() -> None:
    # Exact cap-ten ell=2 values.  Every rational theta on this grid must lie
    # between the two endpoint degree bounds.
    r = Fraction(1, 42)
    p = Fraction(29, 1050)
    lam1, lam2 = Fraction(1, 6), Fraction(1, 70)
    endpoint_pure = (p - lam2) / (1 - lam1 + 10 * (lam1 - lam2))
    for q in range(1, 101):
        for j in range(q + 1):
            theta = Fraction(j, q)
            numerator = (1 - theta) * r + theta * (p - lam2)
            denominator = ((1 - theta)
                           + theta * (1 - lam1 + 10 * (lam1 - lam2)))
            value = numerator / denominator
            assert min(r, endpoint_pure) <= value <= max(r, endpoint_pure)
    assert endpoint_pure == Fraction(14, 2475)
    print(f"K0 endpoint bounds: tautology={r}, pure={endpoint_pure}; grid PASS")


def singleton_obstruction() -> None:
    print("singleton fixed-density audit")
    display = {20, 30, 40, 50, 100, 200}
    checked = 0
    for n in range(20, 301, 5):
        m = 3 * n // 5
        a = Fraction(1, comb(n, m))
        positive = []
        for ell in range(1, m):
            h = Fraction(1, comb(n - ell, m - ell))
            lam2 = eig(n, m, ell, 2)
            if h > lam2:
                positive.append((ell, h - lam2))
            if ell >= 2:
                assert h <= lam2
                ratio_formula = Fraction(
                    m * (m - 1),
                    ell * (ell - 1) * comb(n - ell - 2, n - m - 2),
                )
                assert ratio_formula == h / lam2
        assert positive == [(1, Fraction(1, comb(n - 1, m - 1)))]
        h1 = positive[0][1]
        assert h1 / eig(n, m, 1, 1) == Fraction(
            m, comb(n - 2, n - m - 1)
        )
        checked += 1
        if n in display:
            print(
                f"n={n}, m={m}, degree=1/C(n,m)={a}, "
                f"only positive scale ell=1 with excess={h1}"
            )
    print(f"exhausted n=20,25,...,300: sizes={checked}, every ell PASS")


def stronger_kappa_optimization() -> None:
    # Verify (R51C.6) on the cap-ten ell=2 data.
    r = Fraction(1, 42)
    p = Fraction(29, 1050)
    lam1, lam2 = Fraction(1, 6), Fraction(1, 70)
    d = lam1 - p
    theta_star = r / (r + d)
    kappa_max = (
        r * (lam1 - lam2) / (d + r * (1 - lam1))
    )
    retention = (1 - theta_star) * r + theta_star * p
    assert retention == theta_star * lam1
    direct_kappa = theta_star * (lam1 - lam2) / (1 - theta_star * lam1)
    assert direct_kappa == kappa_max
    print(f"strong-target optimum: theta*={theta_star}, kappa_max={kappa_max}")


if __name__ == "__main__":
    cap10_table()
    k0_endpoint_check()
    singleton_obstruction()
    stronger_kappa_optimization()
    print("PASS: exact spectral-excess and repair identities reproduced")
