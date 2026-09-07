"""Explore cap-nonincreasing two-edge trades with frozen actual children.

The cap-30 Hadamard witness has four neutral two-edge trades, all leaving
the orthogonal bridge class. Explore their component, stopping on cap28.
No reachability or infeasibility claim is made after a state-count cutoff.
"""
import itertools
import json
from pathlib import Path
import numpy as np
from flatify_adversary_2026_09_07_hadamard_bridge_phase_search import spins


def run():
    base=Path(__file__).resolve().parent/'results'
    source=json.loads((base/'flatify_adversary_2026_09_07_hadamard_bridge_phase_search.json').read_text())
    c=np.array(source['records'][1]['best']['parent_matrix'],dtype=np.int16)
    x=spins(8); a=c[:8,:8]; d=c[8:,8:]; initial=c[:8,8:].ravel()
    ha=np.sum((x@a)*x,axis=1)//2; hd=np.sum((x@d)*x,axis=1)//2
    fixed=np.abs(ha[:,None]+hd[None,:]).ravel()
    f=np.array([np.outer(x[:,i],x[:,j]).ravel() for i in range(8) for j in range(8)],dtype=np.int16)
    pairs=list(itertools.combinations(range(64),2)); visited={0:None}; queue=[0]
    witness=None; processed=0; neutral_edges=0; max_states=5000
    while processed<len(queue) and processed<max_states and witness is None:
        code=queue[processed]; processed+=1
        b=initial*np.array([1-2*((code>>i)&1) for i in range(64)],dtype=np.int16)
        bridge=b@f; score=fixed+np.abs(bridge)
        active=np.flatnonzero(score>=26); fa=f[:,active]; ea=bridge[active]; ia=fixed[active]
        changes=-2*b[:,None]*fa
        for i,j in pairs:
            newcode=code^(1<<i)^(1<<j)
            if newcode in visited: continue
            partial=int(np.max(ia+np.abs(ea+changes[i]+changes[j])))
            if partial>30: continue
            full=bridge-2*b[i]*f[i]-2*b[j]*f[j]
            cap=int(np.max(fixed+np.abs(full)))
            assert cap<=30
            visited[newcode]=(code,i,j); queue.append(newcode); neutral_edges+=1
            if cap<=28:
                witness=(newcode,cap); break
        if processed%100==0:
            print(json.dumps(dict(processed=processed,discovered=len(queue))),flush=True)
    path=[]; bestcode=0
    if witness is not None:
        bestcode,cap=witness; p=bestcode
        while visited[p] is not None:
            p0,i,j=visited[p]; path.append([i,j]); p=p0
        path.reverse()
    else: cap=30
    b=(initial*np.array([1-2*((bestcode>>i)&1) for i in range(64)],dtype=np.int16)).reshape(8,8)
    parent=np.block([[a,b],[b.T,d]]); y=spins(16)
    assert int(np.max(np.abs(np.sum((y@parent)*y,axis=1)//2)))==cap
    out=dict(status='COMPONENT EXHAUSTED' if processed==len(queue) else 'BOUNDED SEARCH',
             processed=processed,discovered=len(queue),cap=cap,trade_path=path,
             parent=parent.tolist(),bridge=b.tolist(),bridge_gram=(b@b.T).tolist(),
             visited_predecessors={str(k):v for k,v in visited.items()})
    (base/'flatify_adversary_2026_09_07_neutral_bridge_trade_graph.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('parent','bridge','bridge_gram','visited_predecessors')}),flush=True)


if __name__=='__main__': run()
