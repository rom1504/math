#!/usr/bin/env python3
"""Independent exact checks of band Rayleigh and trace-four formulas."""

from __future__ import annotations

import itertools
import json
import math
import random

import numpy as np


def choose(n: int, k: int) -> int:
    return math.comb(n, k) if 0 <= k <= n else 0


def operator(a: np.ndarray, levels: tuple[int, ...]) -> tuple[np.ndarray, list[tuple[int, ...]]]:
    n = len(a)
    sets = [s for k in levels for s in itertools.combinations(range(n), k)]
    pos = {s: i for i, s in enumerate(sets)}
    t = np.zeros((len(sets), len(sets)), dtype=np.int64)
    for row, ss in enumerate(sets):
        s = set(ss)
        for i in range(n):
            for j in range(i + 1, n):
                u = tuple(sorted(s.symmetric_difference({i, j})))
                if u in pos:
                    t[row, pos[u]] = int(a[i, j])
    return t, sets


def energy(a: np.ndarray, x: np.ndarray) -> int:
    return int(x @ a @ x) // 2


def cycle_sum(a: np.ndarray) -> int:
    result = 0
    for q in itertools.combinations(range(len(a)), 4):
        i, j, k, l = q
        result += int(a[i,j]*a[j,k]*a[k,l]*a[l,i])
        result += int(a[i,j]*a[j,l]*a[l,k]*a[k,i])
        result += int(a[i,k]*a[k,j]*a[j,l]*a[l,i])
    return result


def trace4_formula(a: np.ndarray, k: int) -> int:
    n = len(a)
    c = 32*choose(n-4,k-2)+8*choose(n-4,k-1)+8*choose(n-4,k-3)
    johnson = 0
    for j in range(min(k,n-k)+1):
        multiplicity = choose(n,j)-choose(n,j-1)
        theta = (k-j)*(n-k-j)-j
        johnson += multiplicity*theta**4
    return johnson-3*c*choose(n,4)+c*cycle_sum(a)


def band_rho(n: int, levels: tuple[int, ...]) -> tuple[int, int]:
    d = sum(choose(n,k) for k in levels)
    boundary = choose(n-2,levels[0]-2)+choose(n-2,levels[-1])
    return d-boundary,d


def main() -> None:
    rng = random.Random(20260813)
    trace_tests=[]
    band_tests=[]
    for n in range(4,10):
        for sample in range(4):
            a=np.zeros((n,n),dtype=np.int64)
            for i in range(n):
                for j in range(i+1,n):
                    a[i,j]=a[j,i]=rng.choice((-1,1))
            for k in range(n+1):
                t,_=operator(a,(k,))
                actual=int(np.trace(t@t@t@t))
                predicted=trace4_formula(a,k)
                assert actual==predicted,(n,sample,k,actual,predicted)
                trace_tests.append([n,sample,k,actual])
            for parity in (0,1):
                available=[k for k in range(n+1) if k%2==parity]
                for lo in range(len(available)):
                    for hi in range(lo,len(available)):
                        levels=tuple(available[lo:hi+1])
                        t,sets=operator(a,levels)
                        numerator,d=band_rho(n,levels)
                        for spin_sample in range(3):
                            x=np.asarray([rng.choice((-1,1)) for _ in range(n)],dtype=np.int64)
                            v=np.asarray([math.prod(int(x[i]) for i in s) for s in sets],dtype=np.int64)
                            lhs=int(v@t@v)
                            rhs=numerator*energy(a,x)
                            assert lhs==rhs,(n,sample,levels,spin_sample,lhs,rhs)
                        band_tests.append([n,sample,list(levels),numerator,d])
    output={
        "schema":"global-layer-formulas-independent-check-v1",
        "seed":20260813,
        "trace_four_cases":len(trace_tests),
        "band_rho_cases":len(band_tests),
        "spin_checks_per_band_case":3,
        "orders":"4 through 9",
        "status":"all exact integer identities passed",
    }
    print(json.dumps(output,indent=2,sort_keys=True))


if __name__ == "__main__":
    main()
