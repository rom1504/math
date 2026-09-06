#!/usr/bin/env python3
"""Exact integer replay for the fixed-mask iteration obstruction."""

import argparse
import json
from fractions import Fraction
from pathlib import Path


def decode(text):
    return [{"+": 1, "-": -1, "0": 0}[x] for x in text if not x.isspace()]


def hadamard_entry(i, j, factors):
    value = 1
    for _ in range(factors):
        if i % 4 == j % 4:
            value = -value
        i //= 4
        j //= 4
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    f = decode("0000+000 --+00++- -0--0++0 -+0+-0-+ +--+000+ -00+0+-+ 00+0-0-0 -+0+-0--")
    c = decode("+-+-0-+- 000-+000 0-00+00+ 00+00+00 0000+--0 0--0+000 ++0+0-0- 00+00+00")
    n = 64
    assert len(f) == len(c) == n
    assert all(abs(f[i]) + abs(c[i]) == 1 for i in range(n))
    matrix = [[hadamard_entry(i, j, 3) for j in range(n)] for i in range(n)]
    assert all(matrix[i][j] == matrix[j][i] for i in range(n) for j in range(n))
    assert all(sum(matrix[i][j] * matrix[j][k] for j in range(n)) == (n if i == k else 0)
               for i in range(n) for k in range(n))
    assert all(sum(row) == 8 for row in matrix)
    rf = [sum(row[j] * f[j] for j in range(n)) for row in matrix]
    rc = [sum(row[j] * c[j] for j in range(n)) for row in matrix]
    min_center = min(c[i] * rf[i] for i in range(n) if c[i])
    min_response = min(f[i] * rc[i] for i in range(n) if f[i])
    assert min_center > 0 and min_response > 0
    cross = sum(f[i] * rc[i] for i in range(n))
    ff = sum(f[i] * rf[i] for i in range(n))
    cc = sum(c[i] * rc[i] for i in range(n))
    assert (cross, ff, cc) == (223, 35, -5)
    cross_limit = Fraction(cross, 512)
    endpoint_limit = Fraction(abs(ff + cc) + 2 * cross, 1024)
    assert cross_limit > Fraction(43, 100)
    assert endpoint_limit < Fraction(1, 2)
    x = decode("-+++-+++--+---+-----+++++-+-+-+-++-+----+++++----+-+-+-+----++++")
    y = decode("-++-+--++-+--+-+---------+----+-+++++++++++++++++--+-++---++++--")
    assert len(x) == len(y) == n
    rx = [sum(row[j] * x[j] for j in range(n)) for row in matrix]
    ry = [sum(row[j] * y[j] for j in range(n)) for row in matrix]
    xrx = sum(x[i] * rx[i] for i in range(n))
    yry = sum(y[i] * ry[i] for i in range(n))
    xmargin = min(x[i] * rx[i] for i in range(n))
    ymargin = min(-y[i] * ry[i] for i in range(n))
    assert (xrx, yry, xmargin, ymargin) == (448, -464, 4, 4)
    variable_cross = Fraction(xrx - yry, 2048)
    variable_endpoint = Fraction(max(abs(xrx), abs(yry)), 1024)
    assert variable_cross == Fraction(57, 128) > Fraction(43, 100)
    assert variable_endpoint == Fraction(29, 64) < Fraction(1, 2)
    result = {
        "matrix": "(J4-2I4)^tensor3",
        "f": f,
        "c": c,
        "mask_mass": sum(abs(x) for x in c),
        "R_f": rf,
        "R_c": rc,
        "minimum_signed_center_field": min_center,
        "minimum_signed_response_field": min_response,
        "f_R_c": cross,
        "f_R_f": ff,
        "c_R_c": cc,
        "tensor_lift_bilinear_limit": str(cross_limit),
        "tensor_lift_best_endpoint_limit": str(endpoint_limit),
        "tensor_lift_actual_optimum_limit": "1/2",
        "variable_mask_strict_plateau": {
            "x": x,
            "y": y,
            "R_x": rx,
            "R_y": ry,
            "x_R_x": xrx,
            "y_R_y": yry,
            "minimum_positive_local_margin": xmargin,
            "minimum_negative_local_margin": ymargin,
            "tensor_lift_bilinear_limit": str(variable_cross),
            "tensor_lift_best_endpoint_limit": str(variable_endpoint),
        },
        "verified": True,
    }
    rendered = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
