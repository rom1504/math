"""Independent exact interval check for the root's clustered cap extraction."""

import json
from fractions import Fraction as F
from math import isqrt


def log_two_interval(terms=40):
    z = F(1, 3)
    partial = 2 * sum((z ** (2 * k + 1) / (2 * k + 1) for k in range(terms)), F(0))
    tail = 2 * z ** (2 * terms + 1) / ((2 * terms + 1) * (1 - z * z))
    return partial, partial + tail


def main():
    p = F(31, 32)
    a = F(91470529542342299, 20460000000000000000)
    delta = F(1, 2 ** 24)
    rho_squared = (1 - 2 * delta) ** 2
    scale = 10 ** 30
    root_floor = isqrt(p.numerator * scale * scale // p.denominator)
    sqrt_low, sqrt_high = F(root_floor, scale), F(root_floor + 1, scale)
    assert sqrt_low ** 2 <= p <= sqrt_high ** 2
    log_low, log_high = log_two_interval()
    # -log(1-delta)=sum_{j>=1}delta^j/j, with a geometric tail bound.
    negative_log_partial = sum((delta ** j / j for j in range(1, 5)), F(0))
    negative_log_tail = delta ** 5 / (5 * (1 - delta))
    entropy_low = 24 * delta * log_low + (1 - delta) * negative_log_partial
    entropy_high = 24 * delta * log_high + (1 - delta) * (negative_log_partial + negative_log_tail)
    old_low = F(1, 2) - a / (8 * sqrt_low)
    old_high = F(1, 2) - a / (8 * sqrt_high)
    new_low = (old_low - sqrt_high * entropy_high / 8) / rho_squared
    new_high = (old_high - sqrt_low * entropy_low / 8) / rho_squared
    assert new_high < old_low
    outward_decimal = F(499432211, 10 ** 9)
    assert new_high < outward_decimal
    print(json.dumps({"delta": str(delta), "old_interval": [float(old_low), float(old_high)],
                      "entropy_interval": [float(entropy_low), float(entropy_high)],
                      "cluster_interval": [float(new_low), float(new_high)],
                      "strict_improvement_lower": float(old_low - new_high),
                      "outward_decimal": str(outward_decimal),
                      "exact_outward_check": new_high < outward_decimal}))


if __name__ == "__main__":
    main()
