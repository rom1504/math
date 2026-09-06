"""Exact rational cosquare seed certificate. No optimizer is required.

All 4096 Boolean inputs and both rational positive-definiteness claims
are checked. The matrix P was suggested by a numerical SDP; this replay
uses that suggestion only as fixed rational data.
"""
from fractions import Fraction
import itertools
import json
import numpy as np


C = np.array([
 [1,-1,1,1,1,-1,1,-1,-1,-1,1,1],
 [-1,1,1,-1,-1,-1,1,-1,1,1,1,1],
 [1,1,1,-1,-1,-1,-1,1,1,-1,1,1],
 [1,-1,-1,1,1,-1,-1,1,1,-1,1,1],
 [1,-1,-1,1,1,1,-1,1,-1,-1,1,1],
 [-1,-1,-1,-1,1,1,1,-1,1,1,1,1],
 [1,1,-1,-1,-1,1,1,-1,-1,1,1,1],
 [-1,-1,1,1,1,-1,-1,1,-1,1,1,1],
 [-1,1,1,1,-1,1,-1,-1,1,-1,1,1],
 [-1,1,-1,-1,-1,1,1,1,-1,1,1,1],
 [1,1,1,1,1,1,1,1,1,1,-1,-1],
 [1,1,1,1,1,1,1,1,1,1,-1,1]], dtype=np.int64)

SCALE = 1000000
P_SCALED = np.array([
 [3870935,39998,91463,-129392,423645,-207876,-659,841234,81808,-307209,-394412,79750],
 [39998,2781750,917884,-886054,-971177,935126,255952,39999,912778,255952,-283910,136699],
 [91463,917884,2405333,837965,-47421,-634482,-63524,91463,817116,-63524,-217618,217618],
 [-129392,-886054,837965,4268849,881005,539171,-335893,-129392,33849,-335892,-360198,121807],
 [423645,-971177,-47421,881005,3919002,-2095,307852,423645,-44007,307851,-418891,85803],
 [-207876,935126,-634482,539171,-2095,4047163,23981,-207876,150395,23981,-422201,44539],
 [-659,255952,-63524,-335893,307852,23981,4302683,-307209,-46300,1066402,-431048,93083],
 [841234,39999,91463,-129392,423645,-207876,-307209,3870935,81808,-659,-394412,79750],
 [81808,912778,817116,33849,-44007,150395,-46300,81808,2402537,-46300,-217187,217187],
 [-307209,255952,-63524,-335892,307851,23981,1066402,-659,-46300,4302684,-431048,93083],
 [-394412,-283910,-217618,-360198,-418891,-422201,-431048,-394412,-217187,-431048,3059240,1361263],
 [79750,136699,217618,121807,85803,44539,93083,79750,217187,93083,1361263,2819078]], dtype=np.int64)


def exact_positive_ldl(matrix):
    n = len(matrix)
    assert np.array_equal(matrix, matrix.T)
    a = [[Fraction(int(z)) for z in row] for row in matrix]
    pivots = []
    for k in range(n):
        pivot = a[k][k]
        assert pivot > 0
        pivots.append(pivot)
        for i in range(k+1,n):
            for j in range(i,n):
                a[j][i] = a[i][j] = a[i][j] - a[i][k]*a[j][k]/pivot
    return pivots


def verify():
    cp = C.copy()
    cp[:10,:10] *= -1
    assert np.array_equal(C,C.T) and np.array_equal(cp,cp.T)
    assert np.all(np.abs(C)==1) and np.all(C[:10,:10].sum(axis=0)==0)
    assert np.array_equal(C @ C, cp @ cp)
    spins = np.array(list(itertools.product((-1,1),repeat=12)),dtype=np.int64)
    caps = []
    for a in (C,cp):
        energies = np.sum((spins @ a)*spins,axis=1)
        q_twice = int(np.max(np.abs(energies)))
        beta = int(np.max(np.sum(np.abs(spins @ a),axis=1)))
        witness = spins[int(np.argmax(np.abs(energies)))].tolist()
        caps.append(dict(q_twice=q_twice,beta=beta,quadratic_witness=witness))
    assert caps[0]['q_twice']==54 and caps[0]['beta']==54
    assert caps[1]['q_twice']==58 and caps[1]['beta']==58
    exact_positive_ldl(P_SCALED-SCALE*C)
    exact_positive_ldl(P_SCALED+SCALE*C)
    majorant_max = int(np.max(np.sum((spins @ P_SCALED)*spins,axis=1)))
    upper = Fraction(majorant_max,2*SCALE)
    assert upper == Fraction(60806691,2000000)
    assert upper < Fraction(1946,64)
    print(json.dumps(dict(order=12,spin_inputs=4096,cosquare=True,
                         full_caps=caps,majorant_P_plus_minus_C_positive=True,
                         T_upper=str(upper),T_upper_decimal=str(float(upper)),
                         needed_outer16_q=1946,
                         scalable_separation_established=False),sort_keys=True))


if __name__ == '__main__':
    verify()
