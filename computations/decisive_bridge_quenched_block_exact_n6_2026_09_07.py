"""Exact globally optimized positive block derivative, for EVERY branch-field h."""
from fractions import Fraction
import numpy as np
from decisive_bridge_actual_cavity_psd_checks_2026_09_07 import objects

edges,chars,signs=objects(6)
chars=chars.astype(np.int64);signs=signs.astype(np.int64)
weights=np.array([2 if (i<3)==(j<3) else 1 for i,j in edges])
E=(signs*weights)@chars.T
# lambda_internal=log2, lambda_cross=log2/2. All E are odd, between -21 and21.
# Multiply each branch partition by 2^(21/2) to obtain the following integers.
wp=np.left_shift(np.int64(1),(E+21)//2)
wm=np.left_shift(np.int64(1),(-E+21)//2)
zp=wp.sum(axis=1);zm=wm.sum(axis=1)
products=zp*zm
ids=np.flatnonzero(products==products.min())
assert len(ids)==12
assert np.all(zp[ids]==zm[ids])
direction=np.array([1 if w==2 else -2 for w in weights])
derivatives=[]
for ix in ids:
    numer=((wp[ix]-wm[ix])@chars)*signs[ix]
    den=int(zp[ix]+zm[ix])
    derivative=Fraction(int(numer@direction),den)
    assert derivative==Fraction(58,119)
    derivatives.append(derivative)
print('Global minimum integer branch product:',int(products.min()))
print('Exactly',len(ids),'global product minimizers; every one is branch-balanced.')
print('All active directional correlations:',str(derivatives[0]))
print('Normalized derivative coefficient multiplying log2:',Fraction(5,96)*derivatives[0])
print('PASS: globally minimized quenched pressure has derivative 145 log(2)/5712>0 for every h>=0.')
