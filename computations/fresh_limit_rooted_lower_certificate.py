"""Exact rational interval certificate for the proposed rooted lower bound.

No floating-point arithmetic is used.  Every operation rounds outwards to a
fixed rational grid.  Pi comes from Machin's identity and alternating arctan
series; exp and the Gaussian integral use finite series with explicit tails.
Run: .venv/bin/python -B computations/fresh_limit_rooted_lower_certificate.py
"""

from fractions import Fraction as F
from math import factorial, isqrt
import json


DIGITS = 60
GRID = 10 ** DIGITS


def floor_grid(x):
    return F((x.numerator * GRID) // x.denominator, GRID)


def ceil_grid(x):
    return -floor_grid(-x)


class I:
    def __init__(self, lo, hi=None):
        lo = F(lo)
        hi = lo if hi is None else F(hi)
        assert lo <= hi
        self.lo = floor_grid(lo)
        self.hi = ceil_grid(hi)

    @staticmethod
    def get(x):
        return x if isinstance(x, I) else I(x)

    def __add__(self, other):
        other = I.get(other)
        return I(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-I.get(other))

    def __rsub__(self, other):
        return I.get(other) - self

    def __mul__(self, other):
        other = I.get(other)
        vals = [self.lo * other.lo, self.lo * other.hi,
                self.hi * other.lo, self.hi * other.hi]
        return I(min(vals), max(vals))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = I.get(other)
        assert other.lo * other.hi > 0
        return self * I(1 / other.hi, 1 / other.lo)

    def __rtruediv__(self, other):
        return I.get(other) / self

    def sqrt(self):
        assert self.lo >= 0
        lower_int = isqrt((self.lo.numerator * GRID * GRID)
                          // self.lo.denominator)
        upper_int = isqrt((self.hi.numerator * GRID * GRID)
                          // self.hi.denominator) + 1
        return I(F(lower_int, GRID), F(upper_int, GRID))

    def json(self):
        return {"lower": decimal(self.lo), "upper": decimal(self.hi)}


def decimal(x):
    scaled = x * GRID
    assert scaled.denominator == 1
    val = scaled.numerator
    sign = "-" if val < 0 else ""
    val = abs(val)
    return sign + str(val // GRID) + "." + str(val % GRID).zfill(DIGITS)


def arctan_reciprocal(q, degree=50):
    # Alternating terms decrease for q>1.  The first omitted term bounds
    # the absolute remainder, without assumptions about machine arithmetic.
    s = sum((F((-1) ** k, (2 * k + 1) * q ** (2 * k + 1))
             for k in range(degree + 1)), F(0))
    rem = F(1, (2 * degree + 3) * q ** (2 * degree + 3))
    return I(s - rem, s + rem)


PI = 16 * arctan_reciprocal(5) - 4 * arctan_reciprocal(239)
INV_SQRT_2PI = 1 / (2 * PI).sqrt()


def exp_negative(u, degree=50):
    u = I.get(u)
    assert 0 <= u.lo <= u.hi <= 1
    term = I(1)
    total = term
    for k in range(1, degree + 1):
        term = term * (-u) / k
        total = total + term
    # Lagrange absolute remainder <= exp(1) u^(d+1)/(d+1)!;
    # exp(1)<3 follows directly from its factorial series.
    rem = 3 * u.hi ** (degree + 1) / factorial(degree + 1)
    return total + I(-rem, rem)


def phi(x):
    x = I.get(x)
    return INV_SQRT_2PI * exp_negative(x * x / 2)


def Phi(x, degree=50):
    x = I.get(x)
    assert 0 <= x.lo <= x.hi <= 1
    u = x * x / 2
    term = I(1)
    total = term
    for k in range(1, degree + 1):
        term = term * (-u) / k
        total = total + term / (2 * k + 1)
    # Integrate the alternating exp(-z²/2) series on [0,x].
    rem = (x.hi * u.hi ** (degree + 1)
           / (factorial(degree + 1) * (2 * degree + 3)))
    integral = x * total + I(-rem, rem)
    return I(F(1, 2)) + INV_SQRT_2PI * integral


def main():
    t = I(F(7, 8))
    p = I(F(8, 125))
    a = 2 * phi(t)
    b = 2 * Phi(t) - 1
    sigma = b.sqrt()
    z = a / sigma
    ell = a * (2 * Phi(z) - 1) + 2 * sigma * phi(z)
    result = ((1 - p) * (1 - p) * a * b
              + p * (1 - p) * ell) / (1 + p * p)
    target = F(849, 2500)  # 0.3396 exactly.
    assert result.lo > target
    print(json.dumps({
        "method": "exact_fraction_outward_intervals",
        "grid_decimal_digits": DIGITS,
        "series_degree": 50,
        "t": "7/8", "p": "8/125",
        "pi": PI.json(), "a": a.json(), "b": b.json(),
        "ell": ell.json(), "lower_bound": result.json(),
        "target": "849/2500",
        "margin_above_target": (result - target).json(),
        "verified": True,
    }, indent=2))


if __name__ == "__main__":
    main()
