"""Directed arithmetic for the 2:1 regular-support bridge corollary."""
from fractions import Fraction as F
import json
from pathlib import Path
import mpmath as mp

from flatify_independent_2026_09_07_ternary_interval import endpoints, rational_iv

mp.iv.dps=50
pL,pR,t=map(rational_iv,[F(24,25),F(12,25),F(97,20)])
entropy=-pR*mp.iv.ln(pR)-(1-pR)*mp.iv.ln(1-pR)+pR*mp.iv.ln(2)
E_right_upper=-(entropy-mp.iv.ln(1+2*mp.iv.exp(-t/pR)))/2
coefficient=(1+((pL+pR)*mp.iv.ln(2)-rational_iv(F(5151,6250))+E_right_upper)/(2*t))/mp.iv.sqrt(pL)
upper=endpoints(coefficient)[1]
assert upper < F(99,100)
result={'status':'PASS: DIRECTED INTERVAL ARITHMETIC','p_left':'24/25','p_right':'12/25','t':'97/20','aspect_ratio':'2','E_left_upper':'-5151/6250','E_right_Fano_upper_interval':list(map(str,endpoints(E_right_upper))),'Bennett_coefficient_interval':list(map(str,endpoints(coefficient))),'certified_less_than':'99/100','display_only_coefficient_upper':float(upper)}
Path('computations/results/flatify_independent_2026_09_07_rectangular_bridge_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
