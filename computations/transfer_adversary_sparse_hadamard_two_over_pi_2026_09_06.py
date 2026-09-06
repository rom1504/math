"""Integer semicircle moments for the fixed-degree 2/pi witness family."""

from fractions import Fraction
import json
import math


def add(a, b, scale=1):
    out = [0]*max(len(a), len(b))
    for i, value in enumerate(a):
        out[i] += value
    for i, value in enumerate(b):
        out[i] += scale*value
    return out


def multiply(a, b):
    out = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out


def moment(k):
    return 0 if k % 2 else math.comb(k, k//2)//(k//2+1)


def expectation(coefficients, shift=0):
    return sum(value*moment(k+shift) for k, value in enumerate(coefficients))


def main():
    previous, current = [0], [1]
    total = [0]
    records = []
    for q in range(33):
        assert expectation(multiply(current, current)) == 1
        total = add(total, current)
        covariance = multiply(total, total)
        covariance[0] += 1
        assert expectation(covariance) == q+2
        assert expectation(covariance, shift=1) == 2*q
        second = expectation(multiply(covariance, covariance))
        assert second >= (q+2)**2
        if q in (0, 1, 2, 4, 8, 16, 32):
            records.append({"degree": q, "P_coefficients_low_to_high": total,
                            "R_mean": q+2, "LR_mean": 2*q,
                            "R_squared_mean": second,
                            "coefficient_before_division_by_pi": str(Fraction(2*q, q+2))})
        previous, current = current, add([0]+current, previous, scale=-1)
    print(json.dumps({"degrees_checked": 33, "witnesses": records,
                      "limiting_coefficient_before_division_by_pi": 2,
                      "status": "all exact orthogonality and moment identities pass"}, indent=2))


if __name__ == "__main__":
    main()
