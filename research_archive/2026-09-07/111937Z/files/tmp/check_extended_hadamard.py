import itertools
import numpy as np


def sylvester(k):
    H = np.array([[1]], dtype=np.int64)
    while len(H) < k:
        H = np.block([[H, H], [H, -H]])
    return H


def build(k):
    H = sylvester(k)
    # nonconstant Hadamard rows indexed by other types, using a cyclic bijection
    vs = {}
    for a in range(k):
        others = [b for b in range(k) if b != a]
        for j, b in enumerate(others, 1):
            vs[a, b] = H[j].copy()
    N = k*k
    A = np.zeros((N+1, N+1), dtype=np.int64)
    for a in range(k):
        sl = slice(a*k, (a+1)*k)
        A[sl, sl] = 1
        A[np.arange(a*k,(a+1)*k),np.arange(a*k,(a+1)*k)] = 0
    for a in range(k):
        for b in range(a+1,k):
            block=np.outer(vs[a,b],vs[b,a])
            A[a*k:(a+1)*k,b*k:(b+1)*k]=block
            A[b*k:(b+1)*k,a*k:(a+1)*k]=block.T
    A[N,:N]=1; A[:N,N]=1
    return A,vs


def norm(A):
    n=len(A)
    bestp=-10**9; bestm=10**9; counts=[0,0]
    vals=[]
    for mask in range(1<<(n-1)):
        x=np.ones(n,dtype=np.int64)
        for i in range(n-1):
            if mask>>i&1:x[i]=-1
        e=int(x@A@x)//2
        if e>bestp:bestp=e;counts[0]=1
        elif e==bestp:counts[0]+=1
        if e<bestm:bestm=e;counts[1]=1
        elif e==bestm:counts[1]+=1
    return bestp,bestm,counts


A,vs=build(4)
print('parent',norm(A))
for i in range(len(A)):
    B=np.delete(np.delete(A,i,axis=0),i,axis=1)
    print(i,norm(B))
