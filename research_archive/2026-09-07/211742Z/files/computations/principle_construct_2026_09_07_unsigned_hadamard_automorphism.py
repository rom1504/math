"""Exact unsigned input automorphisms of the normalized Paley12 frame."""
import json
from pathlib import Path
import networkx as nx
import numpy as np
from principle_director_stratified_selector_check_2026_09_07 import paley12

h = paley12()
g = nx.Graph()
for i in range(12):
    g.add_node(i, kind='point')
for i in range(1, 12):
    for sign in (-1, 1):
        node = 12 + 2*(i-1) + (sign+1)//2
        g.add_node(node, kind='block')
        g.add_edges_from((node,j) for j in range(12) if h[i,j] == sign)
g1, g2 = g.copy(), g.copy()
g1.nodes[0]['kind'] = 'marked'
g2.nodes[1]['kind'] = 'marked'
gm = nx.algorithms.isomorphism.GraphMatcher(
    g1, g2, node_match=lambda a,b: a['kind'] == b['kind'])
assert gm.is_isomorphic()
perm = [gm.mapping[i] for i in range(12)]
num = h[:,perm] @ h.T
assert np.all(num % 12 == 0)
rowmap = num // 12
assert np.array_equal(rowmap @ rowmap.T, np.eye(12, dtype=int))
assert np.count_nonzero(rowmap) == 12 and rowmap[0,0] == 1
# Paley finite translations fix infinity and act transitively on 11 finite ports.
translation = [0] + [1 + (i+1)%11 for i in range(11)]
transnum = h[:,translation] @ h.T
assert np.all(transnum % 12 == 0)
assert np.count_nonzero(transnum) == 12
orbit = {0}
while True:
    grown = orbit | {perm[i] for i in orbit} | {translation[i] for i in orbit}
    if grown == orbit:
        break
    orbit = grown
assert len(orbit) == 12
h24 = np.kron(np.array([[1,1],[1,-1]], dtype=int), h)
generators24 = [list(range(12,24))+list(range(12)),
                perm+[12+i for i in perm],
                translation+[12+i for i in translation]]
for p in generators24:
    a = h24[:,p] @ h24.T
    assert np.all(a % 24 == 0)
    assert np.count_nonzero(a) == 24 and a[0,0] == 24
orbit24 = {0}
while True:
    grown = orbit24 | {p[i] for p in generators24 for i in orbit24}
    if grown == orbit24:
        break
    orbit24 = grown
assert len(orbit24) == 24
result = dict(status='PASS', matrix12=h.tolist(),
              infinity_to_zero=perm, output_signed_permutation=rowmap.tolist(),
              translation=translation, orbit12=sorted(orbit),
              generators24=generators24, orbit24=sorted(orbit24),
              scope='Unsigned input permutations; output signed row permutations fix positive DC row.')
path = Path('computations/results/principle_construct_2026_09_07_unsigned_hadamard_automorphism.json')
path.write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result))
