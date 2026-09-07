"""Exact precision identity and directed entropy margin for an actual sign sector."""
from fractions import Fraction as F
from pathlib import Path
import json
import sympy as sp
import mpmath as mp
from flatify_independent_2026_09_07_ternary_interval import endpoints, rational_iv

G=sp.Matrix([[1,sp.Rational(4,5)],[sp.Rational(4,5),1]])
K=sp.Matrix([[0,sp.Rational(2,5)],[sp.Rational(2,5),sp.Rational(37,50)]])
Sigma=G.row_join(K).col_join(K.row_join(G))
P=Sigma.inv()
T=-P[:2,2:]
L=P[:2,:2]-G.inv()
assert T==sp.diag(-sp.Rational(8,5),sp.Rational(10,3))
assert L==sp.Matrix([[sp.Rational(64,45),-sp.Rational(16,9)],[-sp.Rational(16,9),sp.Rational(35,9)]])
assert P[2:,2:]-G.inv()==L
assert Sigma.det()/G.det()**2==sp.Rational(25,108)
assert all(Sigma[:j,:j].det()>0 for j in range(1,5))
center_block=sp.Matrix([[sp.Rational(21,5),sp.Rational(8,5)],[sp.Rational(8,5),sp.Rational(21,5)]])
noise_block=sp.Rational(10,3)*sp.Matrix([[2,-1],[-1,2]])
coupling_block=-4*sp.eye(2)
assert center_block-coupling_block*noise_block.inv()*coupling_block.T==sp.eye(2)
assert (sp.Rational(25,9)*sp.eye(2)).det()/noise_block.det()==sp.Rational(25,108)
mp.iv.dps=50
d=rational_iv(F(1,10))
rate=mp.iv.ln(rational_iv(F(108,25)))/2
entropy=2*(-d*mp.iv.ln(d)-(1-d)*mp.iv.ln(1-d))
margin=rate-entropy
assert endpoints(margin)[0]>F(8,100)
upper=rational_iv(F(493608094,10**9))
parent_margin=2*mp.iv.sqrt(2)*upper-rational_iv(F(37,50))-rational_iv(F(32,25))*upper
assert endpoints(parent_margin)[0]>F(2,100)
result={'status':'PASS: EXACT MATRIX AND DIRECTED INTERVAL CERTIFICATE',
        'alpha':'4/5','bridge_threshold':'37/50',
        'cross_source':str(T),'marginal_gram_penalty':str(L),
        'determinant_ratio':'25/108',
        'arbitrary_center_schur_complement':'EXACT I_2',
        'conditional_rate_interval':list(map(str,endpoints(rate))),
        'rate_minus_two_cloud_entropy':list(map(str,endpoints(margin))),
        'restricted_parent_margin_at_current_upper':list(map(str,endpoints(parent_margin))),
        'display_only_rate':float(endpoints(rate)[0])}
Path('computations/results/flatify_independent_2026_09_07_joint_gaussian_bridge_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
