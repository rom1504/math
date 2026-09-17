"""Exact centered-energy scope counterexample; not a full-sign construction."""
from __future__ import annotations

from fractions import Fraction
from itertools import product
from math import comb, floor, pi, sqrt
from pathlib import Path
import json


def balanced_response(n: int, plus: int) -> Fraction:
    r = n // 2
    denominator = comb(n, r)
    value = 0
    for j in range(max(0, r-(n-plus)), min(plus, r)+1):
        # X has r plus signs, j of which fall in h's plus set.
        overlap = 4*j-2*plus
        value += abs(overlap)*comb(plus, j)*comb(n-plus, r-j)
    return Fraction(value, denominator)


def main() -> None:
    reports = []
    exact_classes = 0
    exhaustive_columns = 0
    two_block_states = 0
    for n in range(2, 102, 2):
        r = n//2
        f = balanced_response(n, r)
        assert f == Fraction(4*r*comb(r-1, floor(r/2))**2, comb(2*r, r))
        value = n*f/(n+f)
        alpha = f/(n+f)
        assert alpha*n == (1-alpha)*f == value
        for plus in range(n+1):
            t = abs(2*plus-n)
            response = balanced_response(n, plus)
            assert response >= f*(1-Fraction(t, n))
            assert alpha*t+(1-alpha)*response >= value
            exact_classes += 1
        if n <= 10:
            words = list(product((-1, 1), repeat=n))
            balanced = [x for x in words if sum(x) == 0]
            for h in words:
                response = Fraction(sum(abs(sum(a*b for a,b in zip(h,x))) for x in balanced), len(balanced))
                assert response == balanced_response(n, (n+sum(h))//2)
                assert max(abs(sum(h)), float(response)) >= float(value)
                # The complete absolute energy ground code is exact.
                twice_scaled_energy = 2*sum(h)**2-n*n
                is_ground = abs(twice_scaled_energy) == n*n
                assert is_ground == (abs(sum(h)) == n or sum(h) == 0)
                exhaustive_columns += 1
        if n <= 12 or n in (20, 50, 100):
            reports.append({"n": n, "f_exact": str(f), "game_exact": str(value),
                            "normalized_game": float(value)/sqrt(n),
                            "sqrt_n_gap": sqrt(2/pi)*sqrt(n)-float(value)})
    for m in (4, 6, 8):
        for ell in (m, m+2, m+4):
            a = Fraction(m*(m-2), ell*(ell-2))
            cap = Fraction(m*(m-1), 2)+a*ell/2
            assert cap == Fraction(m, 2)+a*ell*(ell-1)/2
            for sx in range(-m, m+1, 2):
                for sy in range(-ell, ell+1, 2):
                    energy = (sx*sx-m-a*(sy*sy-ell))/2
                    assert abs(energy) <= cap
                    assert (abs(energy) == cap) == ((abs(sx) == m and sy == 0) or (sx == 0 and abs(sy) == ell))
                    two_block_states += 1
    result = {"status": "PASS", "exact_magnetization_classes": exact_classes,
              "exhaustive_physical_columns": exhaustive_columns,
              "exact_two_block_magnetization_states": two_block_states,
              "scope": "weighted quadratics, with and without a constant; NOT full-sign",
              "reports": reports}
    target = Path(__file__).resolve().parents[1]/"tmp/paper_portfolio_2026_09_17/bernoulli/centered_energy_scope.json"
    target.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
