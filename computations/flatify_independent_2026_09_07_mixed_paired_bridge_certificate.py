"""Exact mixed G/P Gaussian kernels and directed uniform entropy margin."""
from fractions import Fraction as F
from pathlib import Path
import json
import sympy as sp
import mpmath as mp
from flatify_independent_2026_09_07_ternary_interval import endpoints, rational_iv

s2=sp.Rational(9,25); t=sp.Rational(10,3); lam=sp.Rational(35,18)
a=1+2*lam*s2; ap=1+4*lam*s2
D=a*a-t*t*s2*s2
E=a*ap-2*t*t*s2*s2
DP=ap*ap-2*t*t*s2*s2
assert a==sp.Rational(12,5) and ap==sp.Rational(19,5)
assert D==sp.Rational(108,25)
assert a*E==sp.Rational(1872,125)
assert DP==sp.Rational(289,25)
assert t*t*s2/(2*a)-lam==-sp.Rational(10,9)
assert 2*lam+2*sp.Rational(16,25)*(t*t*s2/(2*a)-lam)==sp.Rational(37,15)
assert D/a==sp.Rational(9,5) and E/a==sp.Rational(13,5)

mp.iv.dps=50
logD=mp.iv.ln(rational_iv(F(108,25)))
logGP=mp.iv.ln(rational_iv(F(1872,125)))
logPP=mp.iv.ln(rational_iv(F(289,25)))
d=rational_iv(F(1,10))
h=-d*mp.iv.ln(d)-(1-d)*mp.iv.ln(1-d)
A=logD-logGP/2
B=-logD/2+logGP/2-logPP/4
base=logD/2-2*h
assert endpoints(base)[0]>F(8,100)
base_072=base-rational_iv(F(1,15))
assert endpoints(base_072)[0]>F(1,75)
assert endpoints(B)[0]>0
assert endpoints(h-A-B)[0]>F(1,5)
result={'status':'PASS: EXACT KERNELS AND DIRECTED UNIFORM MIXED-SECTOR MARGIN',
        'alpha':'4/5','t':'10/3','lambda':'35/18','bridge_threshold':'18/25',
        'GG_factor':'(108/25)^(-1)','GP_factor':'(1872/125)^(-1/2)',
        'PP_factor':'(289/25)^(-1/2)',
        'special_edge_mean_coefficient':'-10/9',
        'special_edge_determinant_ratio_upper':'sqrt(13/5)',
        'linear_moment_excess_A':list(map(str,endpoints(A))),
        'quadratic_moment_excess_B':list(map(str,endpoints(B))),
        'base_margin':list(map(str,endpoints(base))),
        'base_margin_at_072':list(map(str,endpoints(base_072))),
        'typed_entropy_minus_A_minus_B':list(map(str,endpoints(h-A-B))),
        'uniform_margin_lower':'1/75+(epsilon_left+epsilon_right)/10',
        'scope':'Same ensemble over all G/P labels outside four reserved G fibres; fixed translations; pinned paired-profile sector only'}
Path('computations/results/flatify_independent_2026_09_07_mixed_paired_bridge_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
