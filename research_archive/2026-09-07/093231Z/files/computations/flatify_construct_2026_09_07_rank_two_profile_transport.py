"""Homogeneous finite-profile entropy diagnostic; not an asymptotic cap proof."""
import itertools
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp


def free_energy(v, q, t):
    d = len(v)
    inds = np.array(list(itertools.product(range(d), repeat=4)))
    a, b, c, e = v[inds].T
    energies = []
    for signs in itertools.product([-1, 1], repeat=4):
        aa, bb, cc, ee = a*signs[0], b*signs[1], c*signs[2], e*signs[3]
        energies.append((aa*cc+aa*ee+bb*cc-bb*ee)/np.sqrt(2))
    logk = logsumexp(t*np.array(energies), axis=0)-np.log(16)
    logbase = np.log(q[inds]).sum(axis=1)+logk
    if d == 1:
        return float(logsumexp(logbase)), 0.

    def obj(z0):
        z = np.r_[z0, 0.]
        logw = logbase+z[inds].sum(axis=1)
        logz = logsumexp(logw)
        w = np.exp(logw-logz)
        marginal_sum = sum(np.bincount(inds[:, j], weights=w, minlength=d) for j in range(4))
        return logz-4*q@z, (marginal_sum-4*q)[:-1]

    fit = minimize(obj, np.zeros(d-1), jac=True, method='L-BFGS-B',
                   options=dict(gtol=1e-10, ftol=1e-13, maxiter=2000))
    value, grad = obj(fit.x)
    return float(value), float(np.max(np.abs(grad))) if len(grad) else 0.


def main():
    raw = json.loads(Path('computations/results/flatify_construct_2026_09_07_walsh_profile_census.json').read_text())
    times = [.25, .5, 1., 2., 4., 8., 16., 32., 64.]
    records = []
    for rec in raw['records']:
        k = rec['k']; rows = []
        for row in rec['rows']:
            hist = np.array(row['histogram'])
            inds = np.flatnonzero(hist)
            v = inds/np.sqrt(k)
            q = hist[inds]/k
            values = []
            for t in times:
                f, residual = free_energy(v, q, t)
                values.append(dict(t=t, free_energy=f, residual=residual,
                                   homogeneous_cap_expression=(row['entropy_per_coordinate']+f/4)/t))
            rows.append(dict(histogram=row['histogram'], count=row['count'], values=values))
        common = [max(row['values'][j]['homogeneous_cap_expression'] for row in rows) for j in range(len(times))]
        print(json.dumps(dict(k=k, best_common_t=times[int(np.argmin(common))],
                              min_common_expression=min(common))), flush=True)
        records.append(dict(k=k, rows=rows, common_expressions=common))
    Path('computations/results/flatify_construct_2026_09_07_rank_two_profile_transport.json').write_text(
        json.dumps(dict(status='finite homogeneous profile proxy only; heterogeneous rows and growing-order counts unproved', records=records), indent=2)+'\n')


if __name__ == '__main__':
    main()
