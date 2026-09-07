"""Bounded exact-cap local search inside the Hadamard BRIDGE orbit.

Both actual minimizing order-eight children remain frozen throughout.
Moves are bridge row/column signs, permutations, and valid closed-quad
Hadamard trades. This is distinct from trading a whole doubled matrix.
"""
import itertools
import json
import math
from pathlib import Path
import numpy as np
from flatify_adversary_2026_09_07_hadamard_bridge_phase_search import spins


def run():
    base=Path(__file__).resolve().parent/'results'
    data=json.loads((base/'flatify_adversary_2026_09_07_hadamard_bridge_phase_search.json').read_text())
    seed=np.array(data['records'][1]['best']['parent_matrix'],dtype=np.int16)
    n=8; x=spins(n); a=seed[:n,:n]; d=seed[n:,n:]
    ea=np.sum((x@a)*x,axis=1)//2; ed=np.sum((x@d)*x,axis=1)//2
    assert max(np.max(np.abs(ea)),np.max(np.abs(ed)))==10
    internal=np.abs(ea[:,None]+ed[None,:])
    rng=np.random.default_rng(9070748)
    quads=list(itertools.combinations(range(n),4))

    def score(b):
        costs=internal+np.abs(x@b@x.T)
        cap=int(costs.max())
        soft=cap+.5*math.log(float(np.exp(2*(costs-cap)).sum()))
        return cap,soft,costs

    best_b=seed[:n,n:].copy(); best_cap,best_soft,_=score(best_b)
    history=[]; accepted=0; attempts=0
    for restart in range(8):
        b=best_b.copy(); cap,soft,_=score(b)
        for iteration in range(6000):
            attempts+=1
            proposal=b.copy(); move=int(rng.integers(5))
            if move<2:
                i=int(rng.integers(n))
                if move==0: proposal[i,:]*=-1
                else: proposal[:,i]*=-1
            elif move<4:
                i,j=rng.choice(n,2,replace=False)
                if move==2: proposal[[i,j],:]=proposal[[j,i],:]
                else: proposal[:,[i,j]]=proposal[:,[j,i]]
            else:
                rows=np.array(quads[int(rng.integers(len(quads)))])
                products=np.prod(b[rows,:],axis=0)
                if np.any(products!=products[0]): continue
                normalized=b[rows,:]*b[rows[0],:]
                col=int(rng.integers(n))
                cols=np.flatnonzero(np.all(normalized==normalized[:,col,None],axis=0))
                proposal[np.ix_(rows,cols)]*=-1
                assert np.array_equal(proposal@proposal.T,n*np.eye(n,dtype=np.int16))
            nc,ns,_=score(proposal)
            temperature=.35*(1-iteration/6000)+.02
            if ns<=soft or rng.random()<math.exp(min(0,(soft-ns)/temperature)):
                b=proposal; cap=nc; soft=ns; accepted+=1
            if (cap,soft)<(best_cap,best_soft):
                best_cap=cap; best_soft=soft; best_b=b.copy()
                row=dict(restart=restart,iteration=iteration,move=move,
                         cap=cap,soft=soft,attempt=attempts)
                history.append(row); print(json.dumps(row),flush=True)
            if best_cap<=28: break
        if best_cap<=28: break
    parent=np.block([[a,best_b],[best_b.T,d]])
    xp=spins(16)
    energies=np.sum((xp@parent)*xp,axis=1)//2
    assert int(np.max(np.abs(energies)))==best_cap
    assert np.array_equal(best_b@best_b.T,n*np.eye(n,dtype=np.int16))
    values,counts=np.unique(energies,return_counts=True)
    result=dict(status='BOUNDED SEARCH; EXACT FULL-CUBE CAP OF SAVED OUTPUT',
        attempts=attempts,accepted=accepted,cap=best_cap,soft=best_soft,
        children_caps=[10,10],bridge=best_b.tolist(),parent=parent.tolist(),
        energy_histogram={int(v):int(c) for v,c in zip(values,counts)},history=history)
    print(json.dumps({k:result[k] for k in ['status','attempts','accepted','cap']}),flush=True)
    return result


if __name__=='__main__':
    result=run()
    target=Path(__file__).resolve().parent/'results'/'flatify_adversary_2026_09_07_frozen_hadamard_trades.json'
    target.write_text(json.dumps(result,indent=2)+'\n')
