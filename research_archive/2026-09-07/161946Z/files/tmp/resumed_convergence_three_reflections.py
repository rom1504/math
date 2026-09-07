"""Exploratory complete profile search at sampled binary atom weights."""
import numpy as np

B=np.array([[-1.,2.],[2.,4.]])
ids=np.arange(1<<15,dtype=np.uint64)
F=np.ones((len(ids),16))
F[:,1:]-=2*((ids[:,None]>>np.arange(15,dtype=np.uint64))&1)
F=F.reshape(-1,8,2)
rng=np.random.default_rng(20260906)
best=3.5
points=list(rng.uniform(.02,.5,size=(80,3)))+[[.25,.25,.25],[1/3]*3]
for ix,z in enumerate(points):
    mu=np.array([1.])
    U=np.array([[1.]])
    for p0 in z:
        p=np.array([p0,1-p0])
        mu=np.kron(mu,p)
        U=np.kron(U,2*np.tile(p,(2,1))-np.eye(2))
    fields=np.einsum('ab,tbi,ij->taj',U,F,B,optimize=True)
    vals=(np.abs(fields).sum(axis=2)*mu).sum(axis=1)/2
    j=vals.argmax()
    if vals[j]>best+1e-9:
        best=vals[j]
        print('improve',ix,z,best,F[j].tolist(),flush=True)
    if ix%20==0:
        print('progress',ix,best,flush=True)
print('final',best,flush=True)
