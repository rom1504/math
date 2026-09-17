"""Test the UNIFORM shifted-absolute binomial/normal comparison.

On each lattice interval the error is concave.  Its extrema are at
lattice endpoints or at the Gaussian quantile of the interval CDF.
Thus these finite candidates exhaust all real thresholds analytically.
Normal functions and binomial probabilities are evaluated numerically;
the displayed diagnostic does not replace the independent analytic proof.
"""

import json
import math

import numpy as np
from scipy.special import ndtr, ndtri
from scipy.stats import binom, norm


def check_order(q):
    positions = (2*np.arange(q+1)-q)/math.sqrt(q)
    probabilities = binom.pmf(np.arange(q+1),q,.5)
    cumulative = np.cumsum(probabilities)
    tail_probability = np.cumsum(probabilities[::-1])[::-1]
    tail_moment = np.cumsum((positions*probabilities)[::-1])[::-1]
    candidates = list(positions[positions>=0]) + [0.0]
    for i in range(q):
        probability = float(cumulative[i])
        if 0 < probability < 1:
            threshold = float(ndtri(probability))
            if positions[i] < threshold < positions[i+1] and threshold>=0:
                candidates.append(threshold)
    candidates = np.array(sorted(set(candidates)))
    first_above = np.searchsorted(positions,candidates,side='right')
    sf = np.zeros_like(candidates)
    first_moment = np.zeros_like(candidates)
    exists = first_above<=q
    sf[exists] = tail_probability[first_above[exists]]
    first_moment[exists] = tail_moment[first_above[exists]]
    binomial_abs = candidates+2*(first_moment-candidates*sf)
    gaussian_abs = candidates+2*(norm.pdf(candidates)-candidates*norm.sf(candidates))
    errors = binomial_abs-gaussian_abs
    at = int(np.argmax(np.abs(errors)))
    maximum = float(abs(errors[at]))
    assert maximum <= 3/q+1e-11
    # Conservative elementary small-ball estimate used by the proof.
    max_atom = float(np.max(binom.pmf(np.arange(q),q-1,.5)))
    assert max_atom <= math.sqrt(2/q)+1e-12
    return dict(q=q, candidates=len(candidates),
                maximum_error=maximum, q_times_error=q*maximum,
                maximizing_threshold=float(candidates[at]),
                sign_at_maximum=int(np.sign(errors[at])))


if __name__ == '__main__':
    orders = list(range(1,513)) + [1000,1024,2048,4096,10000]
    checks = [check_order(q) for q in orders]
    selected = [c for c in checks if c['q'] in (1,2,3,4,8,16,64,256,1024,10000)]
    print(json.dumps(dict(status='PASS',orders_checked=len(checks),
                          real_threshold_candidates=sum(c['candidates'] for c in checks),
                          largest_q_times_error=max(c['q_times_error'] for c in checks),
                          examples=selected,
                          scope='Finite candidates exhaust real thresholds; floating-point special functions'),indent=2))
