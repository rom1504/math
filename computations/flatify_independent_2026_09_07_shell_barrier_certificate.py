"""Directed numerical margin for the actual-child shell certificate barrier."""
from fractions import Fraction as F
from pathlib import Path
import json
import mpmath as mp
from flatify_independent_2026_09_07_ternary_interval import endpoints, rational_iv

mp.iv.dps=50
upper=rational_iv(F(493608094,10**9))
width=rational_iv(F(4333221116640807,10**16))
delta=rational_iv(F(1,10))
alpha=1-2*delta
slack=2*mp.iv.sqrt(2)*upper-2*alpha**2*width
rate=-mp.iv.ln(1-slack**2)/2
entropy=2*(-delta*mp.iv.ln(delta)-(1-delta)*mp.iv.ln(1-delta))
margin=entropy-rate
assert endpoints(margin)[0]>F(3,100)
result={'status':'PASS: DIRECTED INTERVAL ARITHMETIC',
        'upper':'493608094/1000000000',
        'width_lower':'4333221116640807/10000000000000000',
        'noise_probability':'1/10',
        'bridge_slack_interval':list(map(str,endpoints(slack))),
        'rate_upper_interval':list(map(str,endpoints(rate))),
        'two_cloud_entropy_interval':list(map(str,endpoints(entropy))),
        'positive_entropy_minus_rate_margin':list(map(str,endpoints(margin))),
        'certified_margin_greater_than':'3/100',
        'display_only':{k:float(endpoints(v)[1]) for k,v in [('slack',slack),('rate',rate),('entropy',entropy)]}}
Path('computations/results/flatify_independent_2026_09_07_shell_barrier_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
