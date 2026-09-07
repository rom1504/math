"""Exact integer tests for arbitrary-block Huang energy composition."""
import itertools
import json
import numpy as np
from decisive_bridge_huang_flat_eigenvectors_2026_09_07 import apply_h

H4 = np.array([1,1,-1,1,1,-1,-1,1,1,1,-1,-1,-1,-1,-1,1], dtype=np.int64)
assert np.array_equal(apply_h(H4,4),2*H4)

def compose(fs, outer=H4):
    ds=[len(f).bit_length()-1 for f in fs]
    out=[]
    for xs in itertools.product(*[range(len(f)) for f in fs]):
        parity=0
        value=1
        for f,x in zip(fs,xs):
            parity=2*parity+bin(x).count('1')%2
            value*=int(f[x])
        out.append(value*int(outer[parity]))
    return np.array(out,dtype=np.int64),sum(ds)

rng=np.random.default_rng(1307)
trials=0
for ds in [(1,1,1,1),(1,2,3,2),(3,2,1,4),(4,4,4,4)]:
    for _ in range(5):
        fs=[rng.choice([-1,1],size=1<<d) for d in ds]
        f,d=compose(fs)
        # Rayleigh = one half the sum of inner Rayleigh quotients.
        lhs=2*int(f@apply_h(f,d))
        rhs=sum(int(g@apply_h(g,di))*(1<<(d-di)) for g,di in zip(fs,ds))
        assert lhs==rhs,(ds,lhs,rhs)
        for g,di in zip(fs,ds):
            local=g*apply_h(g,di)
            totals=[sum(int(local[x]) for x in range(len(g)) if bin(x).count('1')%2==p) for p in [0,1]]
            assert totals[0]==totals[1]
        trials+=1

f16,d=compose([H4]*4)
assert d==16
assert np.array_equal(apply_h(f16,16),4*f16)
f5=np.repeat(H4,2)
assert int(f5@apply_h(f5,5))==2*int(H4@apply_h(H4,4))
product_trials=0
for a,b in [(1,2),(2,2),(2,3),(3,3)]:
    for _ in range(5):
        inner=rng.choice([-1,1],size=1<<a)
        outer=rng.choice([-1,1],size=1<<b)
        out,d=compose([inner]*b,outer)
        lhs=int(out@apply_h(out,d))*(1<<a)*(1<<b)
        rhs=int(inner@apply_h(inner,a))*int(outer@apply_h(outer,b))*(1<<d)
        assert lhs==rhs
        product_trials+=1
for s in [1,2,4,8,16]:
    d=s*s+1
    for u in range(d%2,d+1,2):
        assert 2*s*u<=u*u+s*s-1
print(json.dumps(dict(status='PASS',arbitrary_inner_trials=trials,arbitrary_product_trials=product_trials,flat_dimensions=[4,16],ignored_coordinate_check=True,plateau_lattice_inequality=True,h9_search='UNKNOWN after 60 seconds; arithmetic proof rules out flat eigenvector')))
