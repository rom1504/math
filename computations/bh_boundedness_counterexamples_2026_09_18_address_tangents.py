#!/usr/bin/env python3
"""Exact sparse-coefficient audit of the addressed tangent obstruction."""
from fractions import Fraction
import importlib.util
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
    f_squared16 = t.convolution(f4, f4)
    assert len(f_squared16) == 16
    assert max(bin(i).count("1") for i in f_squared16) == 4
    a = {0: (1, -1), 4: (0, 1), 8: (1, 0)}
    assert sum(t.norm2(z) for z in a.values()) == 3
    assert t.convolution(a, a).get(0, (0, 0)) == (0, 0)
    # Verify all six nonzero A phases occur, and the same maximum tangent
    # displacement is attainable for every phase of the selected F.
    Avalues = set()
    for x in range(32):
        v = t.evaluate(a, x)
        Avalues.add(v)
        Avalues.add((-v[0], -v[1]))
    assert Avalues == {(0, 0)} | {(2*z[0], 2*z[1]) for z in roots}
    for phase in roots:
        displacements = []
        for value in Avalues:
            u = t.mul(value, t.conj(phase))
            cu = t.conj(u)
            displacements.append((u[0]-cu[0], u[1]-cu[1]))
        assert max(t.norm2(z) for z in displacements) == 12
        assert (-2, 4) in displacements and (2, -4) in displacements

    reports = []
    for k in [1, 2, 3]:
        N = 1 << k
        F4N = {}
        FF16N = {}
        A = []
        distinguished = []
        for leaf in range(N):
            shift = k + 6 * leaf
            sign_mask = 1 << (shift + 5)
            distinguished.append(sign_mask)
            A.append({(i << shift) | sign_mask: z for i, z in a.items()})
            for address_mask in range(N):
                sign = (-1) ** bin(address_mask & leaf).count("1")
                for i, z in f4.items():
                    key = address_mask | (i << shift) | sign_mask
                    F4N[key] = t.add(F4N.get(key, (0, 0)), (sign*z[0], sign*z[1]))
                for i, z in f_squared16.items():
                    key = address_mask | (i << shift)
                    FF16N[key] = t.add(FF16N.get(key, (0, 0)), (sign*z[0], sign*z[1]))
        F4N = {i:z for i,z in F4N.items() if z != (0,0)}
        FF16N = {i:z for i,z in FF16N.items() if z != (0,0)}
        assert len(F4N) == 16 * N * N
        assert all(t.norm2(z) == 1 for z in F4N.values())
        assert max(bin(i).count("1") for i in F4N) == k+3
        assert len(FF16N) == 1+15*N*N
        all_distinguished = sum(distinguished)
        accumulated = set()
        support_counts = []
        for j in range(N):
            product = t.convolution(FF16N, {i:t.conj(z) for i,z in A[j].items()})
            Hj16N = {i:(-z[0],-z[1]) for i,z in product.items()}
            for i,z in A[j].items():
                Hj16N[i] = t.add(Hj16N.get(i,(0,0)), (16*N*z[0],16*N*z[1]))
            Hj16N = {i:z for i,z in Hj16N.items() if z != (0,0)}
            assert max(bin(i).count("1") for i in Hj16N) == k+6
            assert all(i & all_distinguished == distinguished[j] for i in Hj16N)
            assert accumulated.isdisjoint(Hj16N)
            accumulated.update(Hj16N)
            denominator = (16*N)**2
            mass = Fraction(sum(t.norm2(z) for z in Hj16N.values()), denominator)
            outside = Fraction(sum(t.norm2(z) for i,z in Hj16N.items() if i not in F4N), denominator)
            assert mass == 6
            assert outside == (Fraction(189,1)-Fraction(69,N))/64
            support_counts.append(len(Hj16N))
        reports.append({"address_bits":k,"leaves":N,"F_support":len(F4N),
                        "Hj_supports":support_counts,"Hj_squared_norm":"6",
                        "outside_squared_norm":str(outside),"actual_degree":k+6})
    print(json.dumps({"status":"PASS","local_phase_alignment_exact":True,
                      "cases":reports},indent=2))


if __name__ == "__main__":
    main()
