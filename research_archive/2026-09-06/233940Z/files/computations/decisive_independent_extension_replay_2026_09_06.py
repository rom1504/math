"""Independent integer replay of witnesses from the C++ all-minimizer scan."""
import json
from pathlib import Path
import numpy as np

CASES = [(3, 0, 3, 1), (4, 1, 4, 0), (5, 13, 4, 1),
         (6, 220, 5, 4), (7, 828, 9, 1), (7, 826, 9, 3),
         (8, 53014, 10, 2), (9, 898008, 12, 1), (9, 898023, 12, 3),
         (9, 6737136, 12, 3)]
out = []
for n, bits, claimed_m, claimed_e in CASES:
    a = np.ones((n, n), dtype=np.int64) - np.eye(n, dtype=np.int64)
    edges = [(i, j) for i in range(1, n) for j in range(i+1, n)]
    for k, (i, j) in enumerate(edges):
        a[i,j] = a[j,i] = 1 - 2 * ((bits >> k) & 1)
    x = np.ones((1 << (n-1), n), dtype=np.int64)
    for k in range(1, n):
        x[:, k] = 1 - 2 * ((np.arange(len(x)) >> (k-1)) & 1)
    q = np.einsum('bi,ij,bj->b', x, a, x) // 2
    m = int(abs(q).max())
    caps = np.max(abs(q)[None,:] + abs(x @ x.T), axis=1)
    ext = int(caps.min())
    assert (m, ext-m) == (claimed_m, claimed_e)
    row = x[int(caps.argmin())]
    parent = np.zeros((n+1, n+1), dtype=np.int64)
    parent[:n,:n] = a
    parent[n,:n] = parent[:n,n] = row
    full = np.ones((1 << n, n+1), dtype=np.int64)
    for k in range(1,n+1):
        full[:,k] = 1 - 2*((np.arange(len(full)) >> (k-1)) & 1)
    parent_m = int(abs(np.einsum('bi,ij,bj->b',full,parent,full)//2).max())
    assert parent_m == ext
    out.append(dict(n=n, free_edge_bits=bits, M=m, E=ext-m,
                    best_row=row.tolist(), matrix=a.tolist(),
                    parent_direct_M=parent_m,
                    energy_histogram={str(v):int((q==v).sum())
                                      for v in np.unique(q)},
                    row_cap_histogram={str(v):int((caps==v).sum())
                                       for v in np.unique(caps)}))
payload = dict(status='exact integer finite witness replay; no asymptotic claim',
               cases=out)
path=Path('computations/decisive_independent_extension_replay_2026_09_06.json')
path.write_text(json.dumps(payload,indent=2)+'\n')
print(json.dumps([dict(n=c['n'],M=c['M'],E=c['E']) for c in out]))
