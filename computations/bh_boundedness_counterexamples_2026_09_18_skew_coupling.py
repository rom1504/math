#!/usr/bin/env python3
"""Exact moments and small skew-conference coupling replay.

All coefficient, moment, and exhaustive four-phase assertions are exact.
The analytic entropy/support argument is proved in the companion text.
"""
from collections import Counter
from fractions import Fraction
import importlib.util
import itertools
import json
from pathlib import Path


def main():
    source = Path(__file__).with_name("bh_boundedness_counterexamples_2026_09_18_tangent.py")
    spec = importlib.util.spec_from_file_location("tangent", str(source))
    t = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(t)
    roots = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
    masks = [0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 16, 17, 18, 20, 24]
    powers = [0, 2, 0, 5, 1, 4, 2, 4, 3, 1, 4, 1, 2, 2, 3, 2]
    f4 = dict(zip(masks, [roots[j] for j in powers]))
    a = {0: (1, -1), 4: (0, 1), 8: (1, 0)}
    histogram = Counter()
    for x in range(32):
        fx4 = t.evaluate(f4, x)
        fx = (fx4[0]//4, fx4[1]//4)
        histogram[t.mul(t.evaluate(a, x), t.conj(fx))] += 1
    assert histogram == {(0,0):8,(2,-2):5,(0,2):5,(2,0):5,
                         (0,-2):3,(-2,0):3,(-2,2):3}
    mean = (sum(z[0]*p for z,p in histogram.items()),
            sum(z[1]*p for z,p in histogram.items()))
    assert mean == (8,0)
    assert sum(t.norm2(z)*p for z,p in histogram.items()) == 96
    second = (0,0)
    for z,p in histogram.items():
        zz = t.mul(z,z)
        second = t.add(second,(p*zz[0],p*zz[1]))
    assert second == (0,0)
    u4 = t.convolution(a,{i:t.conj(z) for i,z in f4.items()})
    assert len(u4) == 24
    assert max(bin(i).count("1") for i in u4) == 3
    K = [[0,1,1,1],[-1,0,1,-1],[-1,-1,0,1],[-1,1,-1,0]]
    for i in range(4):
        for j in range(4):
            assert K[i][j] == -K[j][i]
            assert sum(K[i][l]*K[j][l] for l in range(4)) == (3 if i==j else 0)
    total_second = 0
    maximum = 0
    for indices in itertools.product(list(histogram),repeat=4):
        q = (0,0)
        probability_numerator = 1
        for z in indices:
            probability_numerator *= histogram[z]
        for i in range(4):
            for j in range(i+1,4):
                z = t.mul(indices[i],t.conj(indices[j]))
                zc = t.conj(z)
                q = t.add(q,(K[i][j]*(z[0]-zc[0]),K[i][j]*(z[1]-zc[1])))
        nn = t.norm2(q)
        maximum = max(maximum,nn)
        total_second += probability_numerator*nn
    assert Fraction(total_second,32**4) == 108
    # Divide K by sqrt(3) for the exact orthogonal matrix.
    assert Fraction(total_second,3*32**4) == 36
    assert Fraction(maximum,3) >= 4**2
    print(json.dumps({"status":"PASS","u_support":24,"u_degree":3,
                      "u_mean":"1/4","u_absolute_second_moment":"3",
                      "u_complex_second_moment":"0","phase_tuples_checked":7**4,
                      "orthogonal_N4_squared_L2":"36",
                      "orthogonal_N4_squared_cap":str(Fraction(maximum,3))},indent=2))


if __name__ == "__main__":
    main()
