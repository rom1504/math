"""Exact rational audit of the affine-Gram ceiling for low-cap signings."""

from fractions import Fraction as Q
from decimal import Decimal, getcontext


def atan_bounds(denominator, terms):
    partial = sum((Q((-1) ** k, (2*k+1)*denominator ** (2*k+1))
                   for k in range(terms)), Q(0))
    error = Q(1, (2*terms+1)*denominator ** (2*terms+1))
    return (partial, partial+error) if terms % 2 == 0 else (partial-error, partial)


def decimal(q):
    return str(Decimal(q.numerator) / Decimal(q.denominator))


def main():
    # Machin's identity, with the elementary alternating-series enclosure.
    lo5, hi5 = atan_bounds(5, 60)
    lo239, hi239 = atan_bounds(239, 20)
    pi_lo, pi_hi = 16*lo5-4*hi239, 16*hi5-4*lo239
    assert Q('3.14159265358979323846264338327950288419716939937510') < pi_lo
    assert pi_hi < Q('3.14159265358979323846264338327950288419716939937511')

    c = Q('0.4333221116640807')  # Downward rounding of the certified lower bound.
    u = Q('0.494515125')
    ell = 2*c-u
    ceiling_lo = 1/(2*pi_hi*ell)
    ceiling_hi = 1/(2*pi_lo*ell)
    threshold_lo = 2*c-1/(2*pi_lo*c)
    threshold_hi = 2*c-1/(2*pi_hi*c)
    assert ceiling_hi < Q('0.427687444510305')
    assert c-ceiling_hi > Q('0.005634667153775')
    assert u < threshold_lo
    assert threshold_hi < Q('0.499432220485404')

    getcontext().prec = 55
    print('PASS exact Machin/alternating-series rational enclosure')
    print('opposite cap lower:', decimal(ell))
    print('signed-energy ceiling interval:', decimal(ceiling_lo), decimal(ceiling_hi))
    print('width gap lower:', decimal(c-ceiling_hi))
    print('sufficient upper threshold interval:', decimal(threshold_lo), decimal(threshold_hi))
    print('PASS new upper passes; old upper misses this sufficient threshold')


if __name__ == '__main__':
    main()
