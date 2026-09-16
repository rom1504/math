"""Lower-witness search in the rigorously programmable F7^r phase family."""
import argparse
import numpy as np

parser=argparse.ArgumentParser()
parser.add_argument('--r',type=int,default=2)
parser.add_argument('--restarts',type=int,default=200)
parser.add_argument('--iterations',type=int,default=200)
parser.add_argument('--ratio',type=float,default=2.)
args=parser.parse_args()
r=args.r
shape=(7,)*r
q=7**r
indices=np.array(np.unravel_index(np.arange(q),shape)).T
flat_weights=7**np.arange(r-1,-1,-1)
unused=set(range(1,q))
plus,minus=[],[]
while unused:
    a=min(unused)
    pp=[int(((c*indices[a])%7)@flat_weights) for c in [1,2,4]]
    mm=[int(((c*indices[a])%7)@flat_weights) for c in [3,5,6]]
    plus.append(pp)
    minus.append(mm)
    unused.difference_update(pp+mm)
plus,minus=np.array(plus),np.array(minus)
rr=args.ratio
B=np.array([[-1.,rr],[rr,rr*rr]])
rng=np.random.default_rng(733)
axes=tuple(range(r))


def ft(f):
    return np.fft.fftn(f.reshape(shape+(2,)),axes=axes).reshape((q,2))


def ift(fhat):
    return np.fft.ifftn(fhat.reshape(shape+(2,)),axes=axes).real.reshape((q,2))


record=(-1,None,None,None)
for restart in range(args.restarts):
    F=rng.choice([-1.,1.],(q,2))
    G=rng.choice([-1.,1.],(q,2))
    lam=np.ones(q,dtype=complex)
    previous=-1
    for iteration in range(args.iterations):
        FF,GG=ft(F),ft(G)
        cross=np.sum(GG.conj()*(FF@B),axis=1)
        s=cross[plus].sum(axis=1)
        z=np.divide(s.conj(),np.abs(s),out=np.ones_like(s),where=np.abs(s)>1e-14)
        lam[plus]=z[:,None]
        lam[minus]=z.conj()[:,None]
        G=np.where(ift(lam[:,None]*(FF@B))>=0,1.,-1.)
        F=np.where(ift(lam.conj()[:,None]*(ft(G)@B))>=0,1.,-1.)
        value=np.sum(G*ift(lam[:,None]*(ft(F)@B)))/(2*q)
        if value<=previous+1e-12:
            break
        previous=value
    if value>record[0]+1e-10:
        record=(value,F.copy(),G.copy(),lam.copy())
        print('record',value,'restart',restart,'iterations',iteration,flush=True)
print('BEST',record[0],flush=True)
print('F',record[1].astype(int).tolist(),flush=True)
print('G',record[2].astype(int).tolist(),flush=True)
print('phases',np.angle(record[3][plus[:,0]]).tolist(),flush=True)
