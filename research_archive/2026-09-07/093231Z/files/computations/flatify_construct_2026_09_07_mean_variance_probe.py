"""Two-label Gaussian-surrogate diagnostic, not a certified envelope bound."""
import collections
import json
from pathlib import Path
import numpy as np
from numpy.polynomial.hermite import hermgauss
from flatify_construct_2026_09_07_rank_two_profile_transport import free_energy


records = []
for order in [8, 16]:
    nodes, weights = hermgauss(order)
    nodes *= np.sqrt(2); weights /= np.sqrt(np.pi)
    for a in [0., .25, .5, .75, .9, 1.]:
        values = np.abs(a+np.sqrt(1-a*a)*nodes)
        support = collections.defaultdict(float)
        for v, q in zip(values, weights):
            support[round(float(v), 12)] += float(q)
        v = np.array(sorted(support)); q = np.array([support[x] for x in v])
        f, residual = free_energy(v, q, 8.)
        psi = -4+f/4
        p = (1+a)/2
        entropy = -p*np.log(p)-(1-p)*np.log(1-p) if a < 1 else 0.
        information = np.log(2)-entropy
        row = dict(order=order, posterior_mean=a, joint_potential=psi,
                   information=information, surrogate_reward=psi-information,
                   residual=residual)
        records.append(row); print(json.dumps(row), flush=True)
Path('computations/results/flatify_construct_2026_09_07_mean_variance_probe.json').write_text(
    json.dumps(dict(status='finite quadrature two-label diagnostic only', records=records), indent=2)+'\n')
