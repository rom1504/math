"""Independent exact/interval certificate for root's arbitrary-P upgrade."""
from fractions import Fraction as Q
from pathlib import Path
import json
import mpmath as mp
from flatify_independent_2026_09_07_ternary_interval import endpoints,rational_iv

t=Q(10,3); lg=Q(35,18); lp=Q(5,4); v=Q(9,25); a=1+2*lg*v
assert a==Q(12,5)
assert t*t<8*lp*lp  # t/(2 sqrt2)<lambda_P
assert lp-t*t*v/(2*a)==Q(5,12)
assert (lg-lp)*v==Q(1,4)
mp.iv.dps=45
d=rational_iv(Q(1,10)); h=-d*mp.iv.ln(d)-(1-d)*mp.iv.ln(1-d)
gp=mp.iv.ln(rational_iv(Q(9,5)))-rational_iv(Q(1,2))-h
pp=mp.iv.ln(rational_iv(Q(108,25)))-1-2*h
assert endpoints(gp)[1]<-Q(23,100)
assert endpoints(pp)[1]<-Q(18,100)
margin=mp.iv.ln(rational_iv(Q(108,25)))/2-2*h-rational_iv(Q(1,15))
assert endpoints(margin)[0]>Q(1,75)
result={'status':'PASS: POINTWISE ARBITRARY-P AND DIRECTED ENTROPY MARGINS',
        'lambda_G':str(lg),'lambda_P':str(lp),'t':str(t),
        'GP_constant':'5/12','GP_residual_coefficient':'-5/12',
        'P_row_norm_mean_saving':'1/4',
        'GP_entropy_moment_change':list(map(str,endpoints(gp))),
        'PP_entropy_moment_change':list(map(str,endpoints(pp))),
        'uniform_margin_at_bridge_18_over_25':list(map(str,endpoints(margin)))}
Path('computations/results/flatify_adversary_2026_09_07_arbitrary_p_sector_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
