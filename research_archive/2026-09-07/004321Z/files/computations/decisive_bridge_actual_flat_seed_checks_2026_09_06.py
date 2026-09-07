"""Integer verification of regular-Hadamard flat sectors and L2 thickening."""
import itertools
import json
import numpy as np

k=4
d=4
m=k*d
base=np.ones((4,4),dtype=np.int64)-2*np.eye(4,dtype=np.int64)
H=np.kron(base,base)
assert np.array_equal(H.T@H,m*np.eye(m,dtype=np.int64))
assert np.all(H.sum(axis=0)==4)
seeds=[np.ones((k,k),dtype=np.int64),np.ones((k,k),dtype=np.int64)]
seeds[1][0,1]=seeds[1][1,0]=-1
caps=[max(abs(np.array(c)@A@np.array(c)) for c in itertools.product([-1,1],repeat=k))
      for A in seeds]
assert caps==[16,12]
rng=np.random.default_rng(7234)
t=1.0
checks=[]
for A in seeds:
    for trial in range(8):
        spins=np.ones((d,k,m),dtype=np.int64)
        for a in range(d):
            for i in range(k):
                spins[a,i,rng.integers(m)]=-1
        # Columns are ordered (j,beta), matching k-by-d spectral blocks.
        spectra=spins@H
        blocks=spectra.reshape(d,k,k,d).transpose(0,3,1,2)/4
        energy=sum(float(np.sum(blocks[a,b]**2)) for a in range(d) for b in range(d) if a!=b)
        delta=sum(float(np.sum((blocks[a,b]-1)**2)) for a in range(d) for b in range(d) if a!=b)
        assert energy <=k*k*d*d+1e-12
        assert delta <=4*(1/m)*k*k*d*d+1e-12
        val=0.0
        for a in range(d):
            for b in range(a+1,d):
                U,V=blocks[a,b],blocks[b,a]
                z=2*t*float(np.sum(U*A*V.T))
                val += -t*(float(np.sum(U*U))+float(np.sum(V*V)))+float(np.logaddexp(z,-z))-np.log(2)
        z=2*t*float(A.sum())
        flat=d*(d-1)/2*(-2*t*k*k+float(np.logaddexp(z,-z))-np.log(2))
        bound=8*t*k*k*np.sqrt(1/m)*d*d
        assert abs(val-flat)<=bound+1e-10
        checks.append(abs(val-flat)/bound)
print(json.dumps({"seed_full_caps":[int(c) for c in caps],"tests":len(checks),
                  "max_fraction_of_L2_bound":max(checks)},indent=2))
