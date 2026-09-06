"""Reproducible statistical falsification test of weighted cubic transport.

NOT an asymptotic certificate: errors are Monte Carlo standard errors over
independent complete input vectors.  The seed has a genuinely nonuniform
limiting unmarked variance, including one value below one.
"""
import json
import math
import numpy as np


SEED = np.array([[-1,-1,-1,-1,-1],[-1,-1,-1,-1,1],
                 [-1,-1,1,1,-1],[-1,-1,1,1,1],[-1,1,-1,1,-1]],
                dtype=np.int64)


def summarize(values):
    values = np.concatenate(values)
    return {'mean': float(values.mean()),
            'standard_error_over_inputs': float(values.std(ddof=1)/math.sqrt(len(values)))}


def tensor_operator(depth):
    """Exact tensor algebra, evaluated in floating point; no dense n^2 arrays."""
    s, k = 4**depth, len(SEED)
    n, diagonal = k*s, (-1)**depth
    e = math.sqrt(s/(n-1))*SEED
    f = -diagonal*np.diag(np.diag(SEED))/math.sqrt(n-1)
    p, t = e@e+f@f, e@f+f@e
    ru = t**3/s
    ri = (p+diagonal*t/math.sqrt(s))**3-diagonal*ru/math.sqrt(s)
    vi = e@ri@e+e@ru@f+f@ru@e+f@ri@f
    vu = e@ri@f+f@ri@e+e@ru@e+f@ru@f
    variances = np.diag(vi+diagonal*vu/math.sqrt(s))
    operator_norm = max(np.abs(np.linalg.eigvalsh(e+f)).max(),
                        np.abs(np.linalg.eigvalsh(-e+f)).max())

    def apply_b(values):
        raw = values.reshape(k,s,-1)
        transformed = raw.copy()
        stride = 1
        for _ in range(depth):
            grouped = transformed.reshape(k,-1,4,stride,raw.shape[2])
            transformed = ((grouped.sum(axis=2,keepdims=True)-2*grouped)/2).reshape(raw.shape)
            stride *= 4
        return (np.einsum('ij,jab->iab',e,transformed)
                +np.diag(f)[:,None,None]*raw).reshape(n,-1)

    if depth == 1:
        outer = np.ones((4,4))-2*np.eye(4)
        dense = np.kron(SEED,outer)
        np.fill_diagonal(dense,0)
        dense /= math.sqrt(n-1)
        assert np.allclose(apply_b(np.eye(n)),dense,atol=1e-14)
        q = dense@dense
        expected = np.diag(dense@(q**3)@dense)
        assert np.allclose(np.repeat(variances,s),expected,atol=1e-12)
    return apply_b,np.repeat(variances,s),float(operator_norm)


def main():
    rng = np.random.default_rng(2026090523)
    q_seed = SEED@SEED
    variance_numerators = np.diag(SEED@(q_seed**3)@SEED)
    assert variance_numerators.tolist() == [733,733,733,733,533]
    reports = []
    samples, batch_size = 4096, 128
    mass = math.erf(1/math.sqrt(2))
    for depth in range(1,7):
        n = 5*4**depth
        apply_b,v,operator_norm = tensor_operator(depth)
        assert v.min()>-1e-10
        d = 1/np.sqrt(np.maximum(v,0.25))
        target = mass*float(np.mean(d*v))
        old_values, new_values = [], []
        for _ in range(samples//batch_size):
            spins = (2*rng.integers(0,2,size=(n,batch_size))-1).astype(float)
            g = apply_b(spins)
            h3 = (g**3-3*g)/math.sqrt(6)
            transported = apply_b(h3)
            raw_tree = apply_b(spins*(g*g-1)/math.sqrt(2))
            # Exact output-root collision subtraction for the injective
            # degree-three marked tree, with Rademacher inputs.
            old_tree = (raw_tree-math.sqrt(2)*spins*np.sum(spins*g,axis=0)/(n-1)
                        +2*math.sqrt(2)*g/(n-1))
            mask = (np.abs(old_tree)<1).astype(float)
            weighted = d[:,None]*mask*transported
            old_values.append(np.sum(weighted*apply_b(old_tree),axis=0)/n)
            new_values.append(np.sum(weighted*transported,axis=0)/n)
        old = summarize(old_values)
        new = summarize(new_values)
        reports.append({'n':n,'depth':depth,'samples':samples,
                        'operator_norm':operator_norm,
                        'minimum_local_variance':float(v.min()),
                        'mean_local_standard_deviation':float(np.sqrt(v).mean()),
                        'finite_gram_main_prediction':target,
                        'old_tree_channel_prediction':0,
                        'old_tree_channel':old,'cubic_channel':new,
                        'cubic_finite_size_error':new['mean']-target})
        print(json.dumps(reports[-1]),flush=True)
    report={'status':'heuristic_statistical_regression_not_asymptotic_proof',
            'rng_seed':2026090523,'full_sign_seed':SEED.tolist(),
            'limiting_variance_numerators':variance_numerators.tolist(),
            'limiting_variance_denominator':625,'orders':reports}
    with open('computations/results/fresh_weighted_projection_tensor_regression.json',
              'w',encoding='utf-8') as handle:
        json.dump(report,handle,indent=2)
        handle.write('\n')


if __name__=='__main__':
    main()
