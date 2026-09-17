"""Exact finite replay of the non-Gaussian balanced-slice repair."""
from fractions import Fraction as F
import itertools
import json
import math
import numpy as np

n=12
words=np.asarray(list(itertools.product((-1,1),repeat=n)),dtype=np.int64)
sums=words.sum(axis=1)
weights=[]
for magnetization in sums:
    if magnetization==0:
        weight=F(2,3*math.comb(n,n//2))
    elif abs(magnetization)==6:
        weight=F(1,6*math.comb(n,9))
    else:
        weight=F(0)
    weights.append(weight)
den=math.lcm(*(weight.denominator for weight in weights))
nums=np.asarray([int(weight*den) for weight in weights],dtype=object)
assert sum(nums)==den
assert all(v==0 for v in words.astype(object).T@nums)
cov_num=words.astype(object).T@(nums[:,None]*words)
assert np.array_equal(cov_num,den*np.eye(n,dtype=object))
abs_response=F(int(np.abs(sums).astype(object)@nums),den)
assert abs_response==2
cube_response=F(int(abs(sums).sum()),2**n)
mixed_response=F(7,8)*abs_response+cube_response/8
assert mixed_response==F(4277,2048)
mixed_weights=[F(7,8)*weight+F(1,8*2**n) for weight in weights]
assert all(weight>0 for weight in mixed_weights)
# Exact three-point moment domination.
normal_moment=1
for k in range(1,101):
    normal_moment*=2*k-1
    assert 3**(k-1)<=normal_moment
# Verify regression and conditional-variance constants on every n6 slice,
# with a nontrivial rational-free centered test vector.
test=np.asarray([-5,-3,-1,1,3,5],dtype=np.int64)
swap_count=0
for word in itertools.product((-1,1),repeat=6):
    value=sum(int(a)*b for a,b in zip(test,word))
    increments=[(int(test[i])-int(test[j]))*(word[i]-word[j])
                for i in range(6) for j in range(6) if i!=j]
    assert F(sum(increments),30)==F(2,5)*value
    assert F(sum(d*d for d in increments),30)<=F(8,5)*int(test@test)
    swap_count+=len(increments)
# Exact block-query response formula for equal-size blocks.
block_values=[]
for r in range(1,17):
    a_r=F(sum(math.comb(r,j)*abs(r-2*j) for j in range(r+1)),2**r)
    assert a_r*a_r<=r
    block_values.append({"blocks":r,"rademacher_mean_abs_exact":str(a_r),
                         "normalized_response_diagnostic":float(a_r)/math.sqrt(3*r)})
print(json.dumps({"status":"PASS exact full-cube covariance, support, response, moments and swap identities",
                  "n":n,"slice_support":sum(weight>0 for weight in weights),
                  "full_support_after_repair":len(mixed_weights),
                  "slice_center_response_exact":str(abs_response),
                  "full_support_center_response_exact":str(mixed_response),
                  "center_response_normalized_diagnostic":float(mixed_response)/math.sqrt(n),
                  "swap_increments_checked":swap_count,
                  "three_point_moments_checked":100,
                  "block_responses":block_values},indent=2))
