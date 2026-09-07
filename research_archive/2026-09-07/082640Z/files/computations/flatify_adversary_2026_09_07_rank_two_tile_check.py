"""Integer reconstruction of arbitrary-frame rank-two tiles and Walsh lift."""
import json
import itertools
import math
from pathlib import Path
import numpy as np
from flatify_adversary_2026_09_07_hadamard_bridge_phase_search import hadamard,spins


def tiles(frames):
    m=len(frames); k=len(frames[0]); s=np.zeros((m*k,m*k),dtype=np.int64)
    h=hadamard(2).astype(np.int64)
    for i in range(m):
        for j in range(m):
            block=frames[i][:,2*j:2*j+2]@h@frames[j][:,2*i:2*i+2].T
            assert np.all(np.abs(block)==2)
            s[i*k:(i+1)*k,j*k:(j+1)*k]=block//2
    return s


def canonical(m):
    h=hadamard(m).astype(np.int64)
    r=np.einsum('aj,ci->iajc',h,h).reshape(m*m,m*m)
    return np.kron(r,hadamard(2)),r


def run():
    rng=np.random.default_rng(9070815); records=[]
    for m in (2,4,8):
        k=2*m; n=m*k
        for trial in range(4):
            frames=[]
            for i in range(m):
                f=hadamard(k).astype(np.int64)[rng.permutation(k)][:,rng.permutation(k)]
                f=f*rng.choice((-1,1),size=k)[:,None]*rng.choice((-1,1),size=k)[None,:]
                frames.append(f)
            s=tiles(frames)
            assert np.array_equal(s,s.T) and np.all(np.abs(s)==1)
            assert np.array_equal(s@s,n*np.eye(n,dtype=np.int64))
            assert int(np.trace(s))==0
            cap=None
            if n==8:
                x=spins(n).astype(np.int64)
                cap=int(np.max(np.abs(np.sum((x@s)*x,axis=1)//2)))
                assert cap==10
            records.append(dict(m=m,trial=trial,cap=cap,frames=[f.tolist() for f in frames],
                                full_sign_matrix=s.tolist()))
    small,_=canonical(2); x=spins(8).astype(np.int64)
    e=np.sum((x@small)*x,axis=1)//2; witness=x[int(np.argmax(e))]
    assert int(e.max())==10
    lifts=[]
    for length in (1,2,4):
        big,_=canonical(2*length); _,r=canonical(length); h=hadamard(length).T.ravel().astype(np.int64)
        assert np.array_equal(r@h,length*h)
        tensor=np.kron(r,small); big_witness=np.kron(h,witness)
        # Tensor index (p,q,i,a,d) corresponds to canonical ((p,i),(q,a),d).
        permutation=[]
        for p in range(length):
            for q in range(length):
                for i in range(2):
                    for a in range(2):
                        for d in range(2):
                            permutation.append(((p*2+i)*(2*length)+(q*2+a))*2+d)
        assert np.array_equal(big[np.ix_(permutation,permutation)],tensor)
        energy=int(big_witness@tensor@big_witness//2)
        assert energy==10*length**3
        lifts.append(dict(length=length,energy=energy,witness=big_witness.tolist(),permutation=permutation))
    # Entire scalar domains relevant to the N32 quantized spectral estimate.
    assert all(abs(z)*12<=z*z+32 for z in range(-32,33,4))
    assert all(abs(z)*8<=z*z+12 for z in range(-30,31,4))
    kernel_checks=[]
    for a,b,c,d in itertools.product(range(3),repeat=4):
        vals=np.array([a*c*e*f+a*d*e*g+b*c*h*f-b*d*h*g
                       for e,h,f,g in itertools.product((-1,1),repeat=4)],dtype=np.int64)
        r2=a*a+b*b; s2=c*c+d*d; p=r2*s2
        assert int(np.sum(vals*vals))==16*p
        if p==0: continue
        # Z=sqrt(2)*vals. Fourth moment identity checked rationally.
        un=(a*a-b*b)**2; vn=(c*c-d*d)**2
        dn=r2*r2*s2*s2-un*s2*s2-vn*r2*r2+2*un*vn
        assert 4*int(np.sum(vals**4))==16*(8*p*p-4*dn)
        defect=dn/(p*p)
        for t in (.05,.25):
            direct=float(np.exp(math.sqrt(2)*t*vals).mean())
            args=math.sqrt(2)*t*np.array([a*c,a*d,b*c,b*d],dtype=float)
            exact=float(np.prod(np.cosh(args))-np.prod(np.sinh(args)))
            u=t*math.sqrt(p)
            angular=(math.cosh(u*math.sqrt(2*defect))+(1-defect)*math.cosh(2*u))/(2-defect)
            assert abs(direct-exact)<1e-10*max(1,direct)
            assert direct<=angular+1e-10*max(1,direct)
            kernel_checks.append(dict(magnitudes=[a,b,c,d],t=t,direct=direct,
                                      exact_formula=exact,angular_upper=angular,D=defect))
    out=dict(status='PASS: INTEGER ALGEBRA; COMPLETE N8 CUBES; EXACT TENSOR LIFTS',
             arbitrary_frame_cases=records,canonical_seed=witness.tolist(),tensor_lifts=lifts,
             universal_N32_cap_upper=88,kernel_checks=kernel_checks)
    target=Path(__file__).resolve().parent/'results'/'flatify_adversary_2026_09_07_rank_two_tile_check.json'
    target.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(status=out['status'],cases=len(records),lifts=[(z['length'],z['energy']) for z in lifts])),flush=True)


if __name__=='__main__': run()
