"""Actual cross/fill witnesses and their row profiles. Lower bounds ONLY."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.special import ndtri
from flatify_construct_2026_09_07_rank_two_weave import sylvester,climb


def main():
    p=argparse.ArgumentParser(); p.add_argument('--m',type=int,default=16)
    p.add_argument('--instances',type=int,default=4)
    p.add_argument('--restarts',type=int,default=128)
    args=p.parse_args(); m=args.m; k=2*m; n=m*k
    rng=np.random.default_rng(2026090700+m); h=sylvester(k)
    z=ndtri((1+np.arange(k+1)/k)/2)
    density=np.exp(-z*z/2)/np.sqrt(2*np.pi)
    weights=2*(density[:-1]-density[1:])
    records=[]
    for it in range(args.instances):
        frames=[h[:,rng.permutation(k)]*rng.choice([-1,1],k) for _ in range(m)]
        a=np.zeros((n,n),dtype=np.int64)
        for i in range(m):
            for j in range(i):
                tile=(frames[i][:,2*j:2*j+2]@np.array([[1,1],[1,-1]])
                      @frames[j][:,2*i:2*i+2].T)//2
                assert np.all(np.abs(tile)==1)
                a[i*k:(i+1)*k,j*k:(j+1)*k]=tile
                a[j*k:(j+1)*k,i*k:(i+1)*k]=tile.T
        cross=a.copy()
        filler=h.copy(); np.fill_diagonal(filler,0)
        for i in range(m): a[i*k:(i+1)*k,i*k:(i+1)*k]=filler
        assert np.all(np.abs(a[np.triu_indices(n,1)])==1)
        assert np.array_equal(a,a.T)
        filled_lower,filled_x=climb(a,rng,args.restarts)
        lower,x=climb(cross,rng,args.restarts); x=np.array(x,dtype=np.int64)
        profiles=[]
        for i in range(m):
            coeff=frames[i].T@x[i*k:(i+1)*k]
            mag=np.sort(np.abs(coeff)/np.sqrt(k))
            dgauss=max(0.,2-2*np.dot(mag,weights))
            dbent=float(np.mean((mag-1)**2))
            flat=[]; r=1
            while r<=k:
                flat.append(2-2*float(mag[-r:].sum())/np.sqrt(k*r));r*=4
            profiles.append(dict(coefficients=coeff.tolist(),gaussian_W2_squared=float(dgauss),
                                 fully_flat_W2_squared=dbent,exact_flat_comparison_squared=min(flat)))
        row=dict(index=it,cross_cap_lower=lower,normalized_cross_lower=lower/n**1.5,
                 filled_cap_lower=filled_lower,filled_witness=filled_x,
                 cross_witness_energy=int(x@cross@x)//2,
                 gaussian_delta=float(np.sqrt(np.mean([z['gaussian_W2_squared'] for z in profiles]))),
                 mixed_delta=float(np.sqrt(np.mean([min(z['gaussian_W2_squared'],z['fully_flat_W2_squared']) for z in profiles]))),
                 frames=[f.tolist() for f in frames],matrix=a.tolist(),witness=x.tolist(),profiles=profiles)
        records.append(row)
        print(json.dumps({key:row[key] for key in ['index','cross_cap_lower','normalized_cross_lower','filled_cap_lower','cross_witness_energy','gaussian_delta','mixed_delta']}),flush=True)
    dest=Path(f'computations/results/flatify_director_rank_two_witness_profiles_2026_09_07_m{m}.json')
    dest.write_text(json.dumps(dict(status='heuristic search; exact stored witness energies, NO cap upper certificate',
                                    m=m,k=k,n=n,seed=2026090700+m,records=records),indent=2)+'\n')


if __name__=='__main__':main()
