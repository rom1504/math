"""Solver-free verification and exact rank-two inflation of a rare cloud witness."""
from pathlib import Path
from fractions import Fraction as Fq
import json
import numpy as np
import mpmath as mp
from scipy.linalg import hadamard
from flatify_independent_2026_09_07_ternary_interval import endpoints, rational_iv

pilot=json.loads(Path('computations/results/flatify_independent_2026_09_07_adversarial_cloud_pilot.json').read_text())
r=next(r for r in pilot if r['n']==32 and 'C' in r)
C=np.array(r['C'],dtype=np.int64)
Fs=[np.array(a,dtype=np.int64) for a in r['F']]
Gs=[np.array(a,dtype=np.int64) for a in r['G']]
x,y,x0,y0=[np.array(r[key],dtype=np.int64) for key in ['x','y','x0','y0']]
H2=hadamard(2).astype(np.int64)
for A in Fs+Gs:
    assert np.all(abs(A)==1) and np.array_equal(A@A.T,8*np.eye(8,dtype=np.int64))
reconstructed=np.block([[Fs[i][:,2*j:2*j+2]@H2@Gs[j][:,2*i:2*i+2].T//2 for j in range(4)] for i in range(4)])
assert np.array_equal(C,reconstructed)
assert np.all(abs(C)==1) and np.array_equal(C@C.T,32*np.eye(32,dtype=np.int64))
assert int(x0@C@y0)==0 and int(x@C@y)==152
assert np.count_nonzero(x!=x0)==np.count_nonzero(y!=y0)==3
a,b=9,2
assert x[a]==x0[a] and y[b]==y0[b]
assert x[a]*(C@y)[a]==0 and y[b]*(C.T@x)[b]==2
assert x[a]*C[a,b]*y[b]==1
inflations=[]
for q in [2,4]:
    H=hadamard(q).astype(np.int64)
    R=np.array([[H[t,u]*H[v,s] for u in range(q) for v in range(q)] for s in range(q) for t in range(q)],dtype=np.int64)
    w=np.array([H[t,s] for s in range(q) for t in range(q)],dtype=np.int64)
    assert np.array_equal(R@R.T,q*q*np.eye(q*q,dtype=np.int64))
    assert np.array_equal(R@w,q*w)
    inflated=np.kron(C,R)
    Fbig=[np.kron(A,H) for A in Fs]; Gbig=[np.kron(A,H) for A in Gs]
    blocks=[]
    for i in range(4):
        for s in range(q):
            row=[]
            for j in range(4):
                for u in range(q):
                    row.append(Fbig[i][:,[2*j*q+u,(2*j+1)*q+u]]@H2@Gbig[j][:,[2*i*q+s,(2*i+1)*q+s]].T//2)
            blocks.append(row)
    actual=np.block(blocks)
    perm=np.array([(i*8+aa)*q*q+s*q+t for i in range(4) for s in range(q) for aa in range(8) for t in range(q)])
    assert np.array_equal(actual,inflated[np.ix_(perm,perm)])
    xx,yy,xx0,yy0=[np.kron(z,w) for z in [x,y,x0,y0]]
    n=32*q*q; extra=q*q//5
    rng=np.random.default_rng(20260907+q)
    for attempt in range(10000):
        subset=rng.choice(q*q,extra,replace=False)
        quadratic=int(w[subset]@R[np.ix_(subset,subset)]@w[subset])
        if quadratic>=0: break
    else: raise AssertionError('No nonnegative subset found')
    xx[a*q*q+subset]*=-1; yy[b*q*q+subset]*=-1
    center=int(xx0@inflated@yy0); target=int(xx@inflated@yy)
    assert center==0 and target>=152*q**3-4*q*extra
    assert np.count_nonzero(xx!=xx0)==np.count_nonzero(yy!=yy0)==n//10
    inflations.append({'q':q,'n':n,'center':center,'target':target,'distance':n//10,'subset':subset.tolist(),'subset_quadratic':quadratic,'verified_actual_rank_two_blocks':True})
mp.iv.dps=50
lower=rational_iv(Fq(189,160))/mp.iv.sqrt(2)
assert endpoints(lower)[0]>Fq(83,100)
upper=rational_iv(Fq(493608094,10**9))
child_excess=lower-(2*mp.iv.sqrt(2)-rational_iv(Fq(32,25)))*upper
assert endpoints(child_excess)[0]>Fq(7,100)
result={'status':'PASS: EXACT INTEGER WITNESS AND ARCHITECTURE INFLATIONS',
        'base':{key:r[key] for key in ['C','F','G','x','y','x0','y0']},
        'base_center':0,'base_target':152,'base_distance':3,
        'additional_flip_base_coordinates':[a,b],
        'inflations':inflations,
        'all_scale_normalized_lower':'189/(160 sqrt(2))',
        'directed_lower_interval':list(map(str,endpoints(lower))),
        'actual_optimal_child_excess_interval':list(map(str,endpoints(child_excess))),
        'actual_optimal_child_excess_greater_than':'7/100',
        'certified_greater_than':'83/100'}
Path('computations/results/flatify_independent_2026_09_07_exceptional_cloud_exact.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({key:value for key,value in result.items() if key!='base'},indent=2))
