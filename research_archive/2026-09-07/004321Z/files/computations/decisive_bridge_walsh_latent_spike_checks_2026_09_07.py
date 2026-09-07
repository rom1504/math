"""Finite common-mask conditional-type/Fourier-spike check; diagnostic only."""
import itertools
import json
import math
import numpy as np

k=4
r=1
dimension=18
m=2**dimension
label_count=2**(k+r)
p=.75
ell=int(p*m)
q=[.5,1.0]
b=[.25,.75]
rng=np.random.default_rng(73271)
Y=np.zeros((k,m),dtype=np.int64)
spins=list(itertools.product([-1,1],repeat=k))
for label in range(label_count):
    amp=label&1
    eps=np.array([1 if label&(1<<(i+1)) else -1 for i in range(k)])
    positions=np.arange(label,m,label_count)
    vectors=[]
    zero_count=round(len(positions)*(1-q[amp]))
    if zero_count:
        vectors.append(np.zeros((zero_count,k),dtype=np.int64))
    for signs in spins:
        prob=q[amp]
        for s in signs:
            prob *= (1+b[amp]*s)/2
        count=round(len(positions)*prob)
        assert abs(count-len(positions)*prob)<1e-9
        vectors.append(np.tile(np.array(signs)*eps,(count,1)))
    cell=np.concatenate(vectors)
    assert len(cell)==len(positions)
    rng.shuffle(cell)
    Y[:,positions]=cell.T
assert np.all(np.count_nonzero(Y,axis=1)==ell)
assert np.all(np.any(Y!=0,axis=0)==np.all(Y!=0,axis=0))

fourier=Y.copy()
step=1
while step<m:
    view=fourier.reshape(k,-1,2*step)
    left=view[:,:,:step].copy()
    right=view[:,:,step:].copy()
    view[:,:,:step]=left+right
    view[:,:,step:]=left-right
    step*=2

support=[[1<<(i+1),(1<<(i+1))|1] for i in range(k)]
# epsilon above is negative at zero bit, so the two coefficients are
# minus the usual positive-character coefficients.
coeff=[-(q[0]*b[0]+q[1]*b[1])/2,
       -(q[0]*b[0]-q[1]*b[1])/2]
for i in range(k):
    for j in range(label_count):
        expected=m*coeff[support[i].index(j)] if j in support[i] else 0
        assert fourier[i,j]==expected
spike_blocks=[]
for s in range(2):
    columns=[support[i][s] for i in range(k)]
    block=fourier[:,columns]/math.sqrt(ell)
    assert np.max(np.abs(block-np.diag(np.diag(block))))==0
    spike_blocks.append(block)

reserved={j for row in support for j in row}
self_columns=[j for j in range(label_count) if j not in reserved][:k]
assert len(self_columns)==k
assert np.max(np.abs(fourier[:,self_columns]))==0
remaining=np.array([j for j in range(m) if j not in reserved])
rng.shuffle(remaining)
blocks=(fourier[:,remaining].reshape(k,-1,k).transpose(1,0,2)/math.sqrt(ell)).reshape(-1,k*k)
physical_distortion=p-sum((q[s]*b[s])**2 for s in range(2))/2
v=physical_distortion/p
second=np.mean(np.sum(blocks*blocks,axis=1))
expected_total=(k*m-k*m*(p-physical_distortion)/p)/(m//k-2)
assert abs(second-expected_total)<1e-10
centered=blocks-blocks.mean(axis=0)
cov=centered.T@centered/len(blocks)
result={"m":m,"common_mask_size":ell,"rows":k,
        "nonzero_mean_frequencies_per_row":2,"mean_coefficients":coeff,
        "zero_self_block_columns":self_columns,
        "max_offdiagonal_spike_entry":max(float(np.max(np.abs(B-np.diag(np.diag(B))))) for B in spike_blocks),
        "predicted_bulk_coordinate_variance":v,
        "observed_mean_bulk_coordinate_variance":float(np.trace(cov)/(k*k)),
        "bulk_covariance_frobenius_error":float(np.linalg.norm(cov-v*np.eye(k*k))),
        "bulk_mean_norm":float(np.linalg.norm(blocks.mean(axis=0))),
        "exact_bulk_second_moment":float(second),
        "limiting_bulk_second_moment":k*k*v}
print(json.dumps(result,indent=2))
