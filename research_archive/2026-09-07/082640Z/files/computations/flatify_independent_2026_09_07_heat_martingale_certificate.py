"""Exact rational check of the heat-martingale SK lower threshold.

Analytic dependencies: concavity of sqrt(cos(pi*s/2)), pi<22/7,
and the alternating cosine series. No floating-point pruning is used.
"""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import json


def check():
    numerators = [990, 961, 911, 840, 745, 618, 441]
    nodes = []
    for j, numerator in enumerate(numerators, 1):
        argument_upper = F(22*j, 7*16)
        cosine_lower = sum(
            (-1)**ell * argument_upper**(2*ell) / factorial(2*ell)
            for ell in range(6)
        )
        ordinate_lower = F(numerator, 1000)
        margin = cosine_lower - ordinate_lower**2
        assert margin > 0
        nodes.append(dict(
            node=str(F(j, 8)), argument_upper=str(argument_upper),
            cosine_lower=str(cosine_lower),
            ordinate_lower=str(ordinate_lower),
            square_margin=str(margin),
        ))
    integral_lower = F(1, 16) + sum(F(z, 8000) for z in numerators)
    assert integral_lower == F(3003, 4000)
    assert integral_lower > F(3, 4)
    # sqrt(2)<99/70; this is deliberately a coarse exact bridge floor.
    assert F(99, 70)**2 > 2
    bridge_lower = integral_lower * F(70, 99)
    assert bridge_lower > F(53, 100)
    return dict(
        status="PASS: EXACT RATIONAL THRESHOLD CERTIFICATE",
        nodes=nodes, integral_lower=str(integral_lower),
        bridge_lower=str(bridge_lower),
        assertions=["SK constant > 3/4", "iid balanced bridge floor > 53/100"],
        scope="Depends on the analytic heat-martingale and comparison proofs; not original M_n lower bound",
    )


if __name__ == "__main__":
    result = check()
    target = Path(__file__).resolve().parent / "results" / "flatify_independent_2026_09_07_heat_martingale_certificate.json"
    target.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: val for key, val in result.items() if key != "nodes"}, indent=2))
