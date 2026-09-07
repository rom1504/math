"""Independent exact arithmetic for the mixed pinned G/P sector."""
from fractions import Fraction as Q
from pathlib import Path
import json
import sympy as S
import mpmath as mp
from flatify_independent_2026_09_07_ternary_interval import endpoints,rational_iv

t=S.Rational(10,3); lam=S.Rational(35,18); s2=S.Rational(9,25)
O=S.Matrix([[1,1],[1,-1]])/S.sqrt(2)
K=S.zeros(4); K[:2,2:]=O; K[2:,:2]=O.T
penalty=lam*S.eye(4)-t*K/2
factors={}
for left,right in [('G','G'),('G','P'),('P','P')]:
    lv=[s2,s2] if left=='G' else [2*s2,0]
    rv=[s2,s2] if right=='G' else [2*s2,0]
    det=S.factor((S.eye(4)+2*S.diag(*(lv+rv))*penalty).det())
    factors[left+right]=det
assert factors=={'GG':S.Rational(108,25)**2,
                 'GP':S.Rational(1872,125),'PP':S.Rational(289,25)}
a=1+2*lam*s2
assert a==S.Rational(12,5)
mean_coefficient=S.factor(t*t*s2/(2*a)-lam)
assert mean_coefficient==-S.Rational(10,9)
assert 2*lam+2*S.Rational(16,25)*mean_coefficient==S.Rational(37,15)

mp.iv.dps=45
D=rational_iv(Q(108,25)); E=rational_iv(Q(1872,125)); F=rational_iv(Q(289,25))
A=mp.iv.ln(D)-mp.iv.ln(E)/2
B=-mp.iv.ln(D)/2+mp.iv.ln(E)/2-mp.iv.ln(F)/4
d=rational_iv(Q(1,10)); h=-d*mp.iv.ln(d)-(1-d)*mp.iv.ln(1-d)
margin=mp.iv.ln(D)/2-2*h
assert endpoints(B)[0]>0
assert endpoints(h-A-B)[0]>0
assert endpoints(margin)[0]>Q(8,100)
reduced=margin-rational_iv(Q(1,15))
assert endpoints(reduced)[0]>Q(1,75)
parent_gap=rational_iv(Q(94,100))*(mp.iv.sqrt(2)-rational_iv(Q(16,25)))-rational_iv(Q(18,25))
assert endpoints(parent_gap)[0]>Q(7,1000)
for m in [2,3,4,8]:
    aa,bb,cc,dd=0,1,0,1
    phi=[dd if i==bb else cc for i in range(m)]
    psi=[bb if j==cc else aa for j in range(m)]
    assert all(psi[phi[i]]!=i for i in range(m))
result={'status':'PASS: EXACT DETERMINANTS AND DIRECTED UNIFORM TYPE MARGINS',
        'precision_determinants':{key:str(val) for key,val in factors.items()},
        'mean_coefficient':str(mean_coefficient),
        'A':list(map(str,endpoints(A))),'B':list(map(str,endpoints(B))),
        'heterogeneous_type_margin_coefficient':list(map(str,endpoints(h-A-B))),
        'margin_at_bridge_18_over_25':list(map(str,endpoints(reduced))),
        'parent_gap_at_seed_47_over_100':list(map(str,endpoints(parent_gap)))}
Path('computations/results/flatify_adversary_2026_09_07_mixed_pinned_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
