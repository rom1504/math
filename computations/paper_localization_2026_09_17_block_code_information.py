"""Finite binomial replay for the exact shared-phase block-code theorem."""
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp
from scipy.stats import binom

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"tmp/paper_portfolio_2026_09_17/localization/block_code_information.json"

local_checks=0
for b in range(16,401):
    kk=np.arange(b+1)
    ss=(2*kk-b)/math.sqrt(b)
    keep=np.abs(ss)<=math.sqrt(b)/2+1e-12
    lp=binom.logpmf(kk[keep],b,.5)
    lower=-1-math.log(4)-.5*math.log(b)-5*ss[keep]**2
    assert np.all(lp>=lower-1e-10)
    local_checks+=int(keep.sum())


def band(b,lo,hi):
    # Evaluate only the positive band and its reflected partners.
    rt=math.sqrt(b)
    kk=np.arange(max(0,int(math.floor((b+lo*rt)/2))-2),
                 min(b,int(math.ceil((b+hi*rt)/2))+2)+1)
    kk=np.unique(np.concatenate([kk,b-kk]))
    ss=(2*kk-b)/rt
    keep=(np.abs(ss)>=lo-1e-10)&(np.abs(ss)<=hi+1e-10)
    kk=kk[keep]; ss=ss[keep]
    assert len(kk)>0
    ll=binom.logpmf(kk,b,.5)
    logp=float(logsumexp(ll))
    weights=np.exp(ll-logp)
    return logp,float(weights@(ss*ss))


rows=[]
for eps in (.9,.5,.25,.125,.0625):
    for odd_offset in (0,1):
        b=int(math.ceil(4096/eps**2))+odd_offset
        a=eps/4; t=4/eps
        lp0,v0=band(b,0,a)
        lp1,v1=band(b,t,t+1)
        pi=(1-v0)/(v1-v0)
        assert 0<pi<1
        assert v0<=a*a+1e-8 and v1>=t*t-1e-8
        assert pi<=2/(t*t)
        assert lp0>=math.log(eps)-21/16-math.log(32)-1e-8
        assert lp1>=-1-math.log(16)-5*(t+1)**2-1e-8
        cov=(1-pi)*v0+pi*v1
        assert abs(cov-1)<1e-12
        response=(1-pi)*math.sqrt(v0)+pi*math.sqrt(v1)
        assert response<=3*eps/4+1e-12
        delta=eps/8
        final_response=(1-delta)*response+delta
        assert final_response<=7*eps/8+1e-12
        cost_per_block=-(1-pi)*lp0-pi*lp1
        assert cost_per_block<=math.log(1/eps)+23
        for r in (1,2,17,1000000):
            assert (1-delta)*r*cost_per_block<=r*math.log(1/eps)+24*r
        rows.append(dict(epsilon=eps,block_size=b,hot_probability=pi,
                         cold_variance=v0,hot_variance=v1,variance=cov,
                         full_support_response_bound=final_response,
                         entropy_upper_per_block=(1-delta)*cost_per_block))

OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(dict(status="PASS",local_binomial_checks=local_checks,cases=rows),indent=2)+"\n")
print(json.dumps(dict(status="PASS",local_binomial_checks=local_checks,shared_phase_cases=len(rows),rank_checks=4*len(rows))))
