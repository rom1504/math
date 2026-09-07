"""Exact finite gauge check and elementary multiplicative-tail diagnostics."""
import json
import math
from fractions import Fraction
from pathlib import Path

H = [[1,1,1,1],[1,-1,1,-1],[1,1,-1,-1],[1,-1,-1,1]]
x = [1,1,1,-1]
h = [sum(H[i][j]*x[i] for i in range(4)) for j in range(4)]
gauge = [-1,1,1,1]
hg = [a*b for a,b in zip(h,gauge)]
y = [sum(Fraction(H[i][j]*hg[j],4) for j in range(4)) for i in range(4)]
assert h == [2,2,2,-2]
assert y == [0,0,0,-2]
# Every H4 sign vector with flat spectrum has odd input parity. Two such
# vectors therefore have even Hamming distance and inner product 0 or +/-4.
flat = []
for mask in range(16):
    v = [1 if mask >> j & 1 else -1 for j in range(4)]
    hv = [sum(H[i][j]*v[i] for i in range(4)) for j in range(4)]
    if all(abs(a)==2 for a in hv):
        flat.append(v)
overlaps = sorted({sum(a*b for a,b in zip(u,v)) for u in flat for v in flat})
assert overlaps == [-4,0,4]
result = dict(h=h, gauged_spectrum=hg, inverse_gauged=[str(a) for a in y],
              flat_vectors=len(flat), flat_pair_overlaps=overlaps,
              incompatible_spectral_pair_overlap=sum(a*b for a,b in zip(h,hg)))
print(json.dumps(result,indent=2))
