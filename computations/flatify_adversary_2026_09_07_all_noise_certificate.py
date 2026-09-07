"""Directed all-noise-sector certificate at seed constant 47/100."""
from fractions import Fraction as Q
import json
import math
from pathlib import Path
import mpmath as mp
from scipy.optimize import brentq
from flatify_independent_2026_09_07_ternary_interval import endpoints, rational_iv

mp.iv.dps = 45
SQRT2_LO = Q(1414213562373095,10**15)
assert SQRT2_LO**2 < 2

def enclose(lo, hi):
    return mp.iv.mpf([rational_iv(lo).a, rational_iv(hi).b])

def entropy(x):
    return -x*mp.iv.ln(x)-(1-x)*mp.iv.ln(1-x)

def leaf(lo, hi):
    mid = float((lo+hi)/2)
    delta = (1-mid)/2
    target = 2*(-delta*math.log(delta)-(1-delta)*math.log1p(-delta))+.005
    d = 1-mid**2
    def rate(z):
        return -math.log1p(-z*z/d)-.5*math.log1p(-z*z/(mid*mid))
    z_float = brentq(lambda z: rate(z)-target, 0,
                     min(mid,math.sqrt(d))*(1-1e-10))
    zq = Q(round(z_float*10**12),10**12)
    a, z = enclose(lo,hi), rational_iv(zq)
    d = 1-a*a
    if endpoints(d-z*z)[0] <= 0 or endpoints(a*a-z*z)[0] <= 0:
        return None
    rate_iv = -mp.iv.ln(1-z*z/d)-mp.iv.ln(1-z*z/(a*a))/2
    delta_iv = (1-a)/2
    rate_gap = endpoints(rate_iv-2*entropy(delta_iv))[0]
    bridge = z*(1+a*a-z*z)/a
    budget = rational_iv(Q(94,100))*(rational_iv(SQRT2_LO)-a*a)
    cap_gap = endpoints(budget-bridge)[0]
    if min(rate_gap,cap_gap) <= 0:
        return None
    return {'lo':str(lo),'hi':str(hi),'z':str(zq),
            'rate_gap_lower':str(rate_gap),'cap_gap_lower':str(cap_gap)}

leaves=[]
pending=[(Q(59,100),Q(99,100))]
nodes=0
while pending:
    lo,hi=pending.pop()
    nodes+=1
    out=leaf(lo,hi)
    if out is None:
        assert hi-lo > Q(1,10**8), 'unexpected unresolved cell'
        mid=(lo+hi)/2
        pending.extend([(mid,hi),(lo,mid)])
    else:
        leaves.append(out)

leaves.sort(key=lambda item: Q(item['lo']))
assert Q(leaves[0]['lo'])==Q(59,100)
assert Q(leaves[-1]['hi'])==Q(99,100)
assert all(Q(a['hi'])==Q(b['lo']) for a,b in zip(leaves,leaves[1:]))
low_gap=Q(94,100)*(SQRT2_LO-Q(59,100)**2)-1
high_bound=2*mp.iv.sqrt(rational_iv(1-Q(99,100)**2))+rational_iv(1-Q(99,100)**2)
high_gap=endpoints(rational_iv(Q(94,100)*(SQRT2_LO-1))-high_bound)[0]
assert min(low_gap,high_gap)>0
result={'status':'PASS: DIRECTED UNIFORM SAME-NOISE-LEVEL SECTOR CERTIFICATE',
        'seed_constant':'47/100','nodes':nodes,'leaf_count':len(leaves),
        'low_endpoint_gap':str(low_gap),'high_endpoint_gap':str(high_gap),
        'minimum_middle_rate_gap':str(min(Q(x['rate_gap_lower']) for x in leaves)),
        'minimum_middle_cap_gap':str(min(Q(x['cap_gap_lower']) for x in leaves)),
        'leaves':leaves}
path=Path('computations/results/flatify_adversary_2026_09_07_all_noise_certificate.json')
path.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='leaves'},indent=2))
