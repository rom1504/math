"""Outward-rounded scalar premises for the stability/selector improvement.

This certifies only elementary entropy inequalities. The positive improvement
itself is a symbolic constant from the analytic orthant theorem; no floating
evaluation of its extraordinarily small magnitude is substituted for proof.
"""
from fractions import Fraction
import json
from pathlib import Path

import mpmath as mp


def rational(n, d=1):
    return mp.iv.mpf(n) / d


def entropy(x):
    return -x * mp.iv.ln(x) - (1-x) * mp.iv.ln(1-x)


def main():
    mp.iv.dps = 70
    p = rational(24, 25)
    eps = rational(1, 100000)
    q = rational(1, 10000)
    gap = entropy(p) + p * mp.iv.ln(2) - rational(5151, 6250)
    xi = (entropy(q) + q * mp.iv.ln(1+2/mp.iv.sqrt(eps))
          + entropy(16*eps) + 16*eps*mp.iv.ln(2))
    assert bool(gap > rational(9, 1000))
    assert bool(xi < rational(4, 1000))
    assert bool(xi < gap/2)
    assert Fraction(100000)**2 >= 4*100**2/(Fraction(1, 2)*Fraction(1, 100000))
    assert Fraction(10**10) >= 4*100**2/(Fraction(1, 2)*Fraction(1, 100000))
    uniform_eps = rational(1, 1000)
    uniform_q = rational(1, 400)
    uniform_xi = (entropy(uniform_q)
                  + uniform_q*mp.iv.ln(1+2/mp.iv.sqrt(uniform_eps))
                  + entropy(16*uniform_eps)+16*uniform_eps*mp.iv.ln(2))
    assert bool(uniform_xi < rational(122, 1000))
    assert bool(entropy(p) > rational(167, 1000))
    assert Fraction(2000)**2 >= 4*20**2/Fraction(1, 1000)
    assert Fraction(2000000) >= 4*20**2/Fraction(1, 1000)
    payload = {
        "status": "OUTWARD_ROUNDED_ENTROPY_GAP_PASS", "precision": 70,
        "p": "24/25", "t": "97/20", "epsilon": "1/100000", "V": 100,
        "delta": "1/2", "U": 100000, "C": 10**10,
        "gap_interval": str(gap), "xi_interval": str(xi),
        "proved": ["gap>9/1000", "xi<4/1000", "xi<gap/2"],
        "theta": "1/16000000000",
        "kappa": "sech(97000000)^2/1600000",
        "uniform_good_selector": {
            "epsilon": "1/1000", "V": 20, "delta": 1,
            "U": 2000, "C": 2000000,
            "xi_interval": str(uniform_xi), "selector_entropy": str(entropy(p)),
            "proved": ["xi<.122", "h(p)>.167", "xi<h(p)"],
            "theta": "1/3200000", "kappa": "sech(388000)^2/8000",
            "scope": "No phase-gap assumption: each selector excludes ALL concentrated spin words."
        },
        "scope": "Entropy premises only; stability gain is positive symbolically, not a reported decimal."
    }
    output = Path(__file__).resolve().parent / "results/principle_director_2026_09_07_selector_entropy_certificate.json"
    output.write_text(json.dumps(payload, indent=2)+"\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
