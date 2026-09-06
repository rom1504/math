"""Balanced lower tests for vectorial-MM Walsh-conjugated power permutations."""
import argparse, math
import numpy as np
from scipy.linalg import hadamard

parser=argparse.ArgumentParser()
parser.add_argument('--m',type=int,default=6)
parser.add_argument('--restarts',type=int,default=100)
args=parser.parse_args()
m=args.m
q=2**m
polynomials={4:0b10011,5:0b100101,6:0b1000011,7:0b10000011,
             8:0b100011101,9:0b1000010001,10:0b10000001001}
modulus=polynomials[m]
power,logarithm=np.zeros(q-1,dtype=int),np.zeros(q,dtype=int)
a=1
for j in range(q-1):
    if j:
        assert a!=1
    power[j]=a
    logarithm[a]=j
    a<<=1
    if a>=q:
        a^=modulus
assert a==1
# Standard bit Walsh is a valid additive dual basis. Field power is then
# conjugated by that identification, equivalent to trace Walsh by a linear
# permutation; its real orthogonality and Boolean-profile norm are unchanged.
H=hadamard(q).astype(float)
B=np.array([[-1.,2.],[2.,4.]])
rng=np.random.default_rng(739)


def balanced_sign(fields):
    out=-np.ones_like(fields)
    indices=np.argsort(fields,axis=0)[q//2:]
    for i in range(2):
        out[indices[:,i],i]=1
    return out


used=set()
exponents=[]
for exponent in range(1,q-1):
    if math.gcd(exponent,q-1)>1 or exponent in used:
        continue
    orbit={(exponent*2**j)%(q-1) for j in range(m)}
    used.update(orbit)
    exponents.append(exponent)
record=(-1,None,None,None)
for exponent in exponents:
    permutation=np.zeros(q,dtype=int)
    permutation[1:]=power[(exponent*logarithm[1:])%(q-1)]
    U=H[:,permutation]@H/q
    assert np.max(np.abs(U.T@U-np.eye(q)))<1e-12
    for restart in range(args.restarts):
        F=balanced_sign(rng.normal(size=(q,2)))
        if restart%2:
            F[:,1]=F[:,0]
            # Flip equal plus/minus subsets to preserve balance and corr .75.
            for sign in [-1,1]:
                choices=np.flatnonzero(F[:,0]==sign)
                flips=rng.choice(choices,q//16,replace=False)
                F[flips,0]*=-1
        previous=-1
        for iteration in range(100):
            G=np.where(U@F@B>=0,1.,-1.)
            F=balanced_sign(U.T@G@B)
            value=np.abs(U@F@B).sum()/(2*q)
            if value<=previous+1e-12:
                break
            previous=value
        if value>record[0]+1e-10:
            record=(value,exponent,F.copy(),G.copy())
            print('record',value,'exponent',exponent,'restart',restart,flush=True)
print('BEST',record[0],record[1],flush=True)
print('F',record[2].astype(int).tolist(),flush=True)
