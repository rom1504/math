"""Outward-rounded constants for the safe Sylvester mixed cone."""
import json
from pathlib import Path
import mpmath as mp

mp.iv.dps = 65
iv = mp.iv.mpf
p,t,b = iv(31)/32,iv(97)/20,iv(499)/1000
fee = -p*mp.iv.ln(p)-(1-p)*mp.iv.ln(1-p)-mp.iv.ln(32)/32
A = p*mp.iv.ln(2)-iv(4)/5+fee
d = 1-2*b*mp.iv.sqrt(p)
sigma = -(A+t*d)
xi,alpha,r0 = iv(1)/200,iv(1)/400,iv(1)/1000000
risk = p*(-r0*mp.iv.ln(r0)-(1-r0)*mp.iv.ln(1-r0))+4*t*d*r0
slack = (sigma-xi)*alpha-(1-alpha)*risk
assert sigma.a > xi.b
assert slack.a > 0
source = Path('computations/results/principle_director_sylvester_stratified_certificate_2026_09_07.json')
certificate = json.loads(source.read_text())
assert certificate['p'] == '31/32' and certificate['t'] == '97/20'
assert certificate['E_upper'] == '-4/5'
def bounds(x):
    return [str(x.a),str(x.b)]
result = dict(status='PASS', precision=65, source=str(source),
              p='31/32', t='97/20', cap='499/1000', xi='1/200',
              balanced_fraction='1/400', other_minority_max='1/1000000',
              sigma_interval=bounds(sigma), risk_interval=bounds(risk),
              cone_slack_interval=bounds(slack),
              scope='Constraint-set variance cap; not all-profile cap or convergence.')
Path('computations/results/principle_construct_2026_09_07_mixed_cone_constants.json').write_text(
    json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
