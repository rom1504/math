"""Probe a precise candidate: subhalf optimal caps forbid isotropic grounds.

LP floating feasibility is reconnaissance only. Positive answers must be
converted to exact rational certificates before mathematical use.
"""
import json
from pathlib import Path
import networkx as nx
import numpy as np
from scipy.optimize import linprog

root = Path(__file__).resolve().parents[1]
n = 8
x = np.column_stack((np.ones(1 << (n-1), dtype=np.int64),
                     1-2*((np.arange(1 << (n-1))[:,None] >>
                           np.arange(n-1)) & 1)))
edges = np.triu_indices(n,1)
products = x[:,edges[0]]*x[:,edges[1]]
records=[]
for label,g in enumerate(nx.graph_atlas_g()):
    if len(g)!=7:
        continue
    a=np.ones((n,n), dtype=np.int64)-np.eye(n,dtype=np.int64)
    adjacency=nx.to_numpy_array(g,dtype=np.int64)
    a[1:,1:]=np.ones((7,7),dtype=np.int64)-np.eye(7,dtype=np.int64)-2*adjacency
    energy=products@a[edges]
    if np.max(np.abs(energy))!=10:
        continue
    ground=np.flatnonzero(np.abs(energy)==10)
    constraints=np.vstack((np.ones(len(ground)),products[ground].T))
    result=linprog(np.zeros(len(ground)),A_eq=constraints,
                   b_eq=np.r_[1,np.zeros(len(edges[0]))],
                   bounds=(0,None),method='highs')
    orthogonal=(x[ground]@x[ground].T==0)
    graph=nx.from_numpy_array(orthogonal)
    clique=max(nx.find_cliques(graph),key=len,default=[])
    rec=dict(atlas=label,ground_count=len(ground),
             floating_isotropy_feasible=bool(result.success),
             orthogonal_ground_clique=len(clique))
    if result.success:
        support=np.flatnonzero(result.x>1e-9)
        rec.update(matrix=a.tolist(),ground_spins=x[ground[support]].tolist(),
                   floating_weights=result.x[support].tolist())
    records.append(rec)
stored = []
for order in range(3, 15):
    path = root/f'computations/results/exact_m{order}.json'
    if not path.exists():
        continue
    raw=json.loads(path.read_text())
    if 'matrix' not in raw:
        continue
    a=np.array(raw['matrix'],dtype=np.int64)
    words=1-2*((np.arange(1 << (order-1))[:,None] >> np.arange(order))&1)
    ei=np.triu_indices(order,1)
    en=(words[:,ei[0]]*words[:,ei[1]])@a[ei]
    cap=int(max(abs(en)))
    grounds=words[abs(en)==cap]
    gram=grounds.T@grounds
    exact_uniform=np.array_equal(gram,len(grounds)*np.eye(order,dtype=np.int64))
    item=dict(order=order,cap=cap,ground_count=len(grounds),
              exact_uniform_isotropy=bool(exact_uniform))
    if exact_uniform:
        item.update(matrix=a.tolist(),ground_spins=grounds.tolist(),
                    ground_gram=gram.tolist())
    stored.append(item)
target=root/'computations/results/principle_director_ground_isotropy_probe_2026_09_07.json'
target.write_text(json.dumps(dict(status='MIXED: LP reconnaissance; integer Gram certificates exact',
                                records=records,stored_signings=stored),indent=2)+'\n')
print(json.dumps(dict(optimal_caps_checked=len(records),
    feasible=sum(r['floating_isotropy_feasible'] for r in records),
    max_orthogonal_clique=max(r['orthogonal_ground_clique'] for r in records),
    exact_isotropic_stored_orders=[r['order'] for r in stored if r['exact_uniform_isotropy']])))
