"""Binary typed-row entropy diagnostic; finite-profile proxy only."""
import itertools
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp


def prepare(vi, qi, vj, qj, t):
    inds = np.array(list(itertools.product(range(len(vi)), range(len(vi)), range(len(vj)), range(len(vj)))))
    a, b = vi[inds[:, 0]], vi[inds[:, 1]]
    c, d = vj[inds[:, 2]], vj[inds[:, 3]]
    energies = []
    for s in itertools.product([-1, 1], repeat=4):
        aa, bb, cc, dd = a*s[0], b*s[1], c*s[2], d*s[3]
        energies.append((aa*cc+aa*dd+bb*cc-bb*dd)/np.sqrt(2))
    base = np.log(qi[inds[:, :2]]).sum(axis=1)+np.log(qj[inds[:, 2:]]).sum(axis=1)
    return inds, base+logsumexp(t*np.array(energies), axis=0)-np.log(16)


def mixed(vs, qs, weights, t):
    sizes = [len(v)-1 for v in vs]
    edges = [(i, j, weights[i]*weights[j]*(1 if i == j else 2),
              *prepare(vs[i], qs[i], vs[j], qs[j], t)) for i in range(2) for j in range(i, 2)]

    def objective(z):
        zs = [np.r_[z[:sizes[0]], 0.], np.r_[z[sizes[0]:], 0.]]
        gs = [-4*weights[i]*qs[i].copy() for i in range(2)]
        value = -4*sum(weights[i]*(qs[i]@zs[i]) for i in range(2))
        for i, j, weight, inds, base in edges:
            lw = base+zs[i][inds[:, :2]].sum(axis=1)+zs[j][inds[:, 2:]].sum(axis=1)
            lz = logsumexp(lw); prob = np.exp(lw-lz)
            value += weight*lz
            for col in range(4):
                kind = i if col < 2 else j
                gs[kind] += weight*np.bincount(inds[:, col], weights=prob, minlength=len(qs[kind]))
        return value, np.r_[gs[0][:-1], gs[1][:-1]]

    fit = minimize(objective, np.zeros(sum(sizes)), jac=True, method='L-BFGS-B',
                   options=dict(gtol=1e-9, ftol=1e-12, maxiter=1000))
    f, g = objective(fit.x)
    return float(f), float(np.max(np.abs(g)))


raw = json.loads(Path('computations/results/flatify_construct_2026_09_07_walsh_profile_census.json').read_text())['records'][-1]
k = raw['k']; profiles = raw['rows']; records = []; times = [4., 8., 16., 32.]
for i, j in itertools.combinations(range(len(profiles)), 2):
    vs = []; qs = []
    for index in [i, j]:
        hist = np.array(profiles[index]['histogram']); support = np.flatnonzero(hist)
        vs.append(support/np.sqrt(k)); qs.append(hist[support]/k)
    for p in [.25, .5, .75]:
        weights = [p, 1-p]
        entropy = sum(weights[l]*profiles[index]['entropy_per_coordinate'] for l, index in enumerate([i, j]))
        values = []
        for t in times:
            f, residual = mixed(vs, qs, weights, t)
            values.append(dict(t=t, free_energy=f, residual=residual,
                               cap_expression=(entropy+f/4)/t))
        records.append(dict(types=[i, j], p=p, entropy=entropy, values=values,
                            best_expression=min(v['cap_expression'] for v in values)))
    print(json.dumps(dict(pair=[i, j], worst_best=max(r['best_expression'] for r in records[-3:]))), flush=True)
common = [max(r['values'][j]['cap_expression'] for r in records) for j in range(len(times))]
Path('computations/results/flatify_construct_2026_09_07_rank_two_mixed_transport.json').write_text(
    json.dumps(dict(status='finite k16 two-type proxy; not asymptotic cap proof', times=times,
                    common_expressions=common, records=records), indent=2)+'\n')
print(json.dumps(dict(common_expressions=common, min_common=min(common))), flush=True)
