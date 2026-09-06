"""Independent bilinear search, with exact H4 quadratic embedding."""
import argparse
import json
import time
import numpy as np
from fresh_cosquare_lift_spin_search import pair


def main(power, seconds, seed, batch, target):
    core=np.array(pair()['flipped'], dtype=np.int16)
    h4=np.ones((4,4),dtype=np.int16)-2*np.eye(4,dtype=np.int16)
    outer=np.ones((1,1),dtype=np.int16)
    for _ in range(power):
        outer=np.kron(outer,h4)
    a=np.kron(outer,core)
    n=len(a)
    rng=np.random.default_rng(seed)
    best=-1
    bx=by=None
    start=time.monotonic()
    tries=0
    while time.monotonic()-start<seconds and best<target:
        x=rng.choice(np.array([-1,1],dtype=np.int16),size=(batch,n))
        if bx is not None:
            x[:batch//2]=bx
            flips=rng.random((batch//2,n))<rng.choice([.01,.03,.08,.15,.3])
            x[:batch//2]*=1-2*flips.astype(np.int16)
        for iteration in range(60):
            field=x@a
            y=np.where(field>=0,1,-1).astype(np.int16)
            vals=np.sum(np.abs(field).astype(np.int64),axis=1)
            idx=int(np.argmax(vals))
            if int(vals[idx])>best:
                best=int(vals[idx]);bx=x[idx].copy();by=y[idx].copy()
                print(json.dumps(dict(elapsed=round(time.monotonic()-start,3),
                                      n=n,beta=best,normalized_seed_q=best/(2*len(outer)**1.5))),flush=True)
            updated=np.where(y@a>=0,1,-1).astype(np.int16)
            if np.array_equal(updated,x) or best>=target:
                break
            x=updated
        tries+=batch
    assert int(bx.astype(np.int64)@a@by.astype(np.int64))==best
    # z=(x,x,y,y) and z'=(x,-x,y,-y) provide signed cross terms;
    # the two useful choices below have opposite self-energy sums.
    candidates=[]
    lifted=np.kron(h4,a)
    for perm in ((bx,bx,by,by),(bx,-bx,by,-by),(bx,by,bx,by),
                 (bx,by,-bx,-by),(bx,-by,bx,-by),(-bx,by,-bx,by)):
        z=np.concatenate(perm)
        q2=int(z.astype(np.int64)@lifted@z.astype(np.int64))
        candidates.append((abs(q2),q2,z))
    qabs,q2,z=max(candidates,key=lambda item:item[0])
    # Direct ±y versions of the repeated blocks always certify q>=4 beta.
    assert qabs>=8*best
    print(json.dumps(dict(method='alternating_bilinear_exact_quadratic_embedding',
                          power=power,n=n,elapsed=time.monotonic()-start,tries=tries,
                          beta=best,x=bx.tolist(),y=by.tolist(),
                          lifted_order=4*n,lifted_q_twice=q2,lifted_spin=z.tolist(),
                          target=target,target_reached=best>=target)),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--power',type=int,default=2)
    parser.add_argument('--seconds',type=float,default=120)
    parser.add_argument('--seed',type=int,default=2026090599)
    parser.add_argument('--batch',type=int,default=128)
    parser.add_argument('--target',type=int,default=3892)
    args=parser.parse_args()
    main(args.power,args.seconds,args.seed,args.batch,args.target)
