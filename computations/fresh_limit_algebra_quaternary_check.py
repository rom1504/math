"""Exact integer checks of the restricted quadratic-phase no-go theorem."""
import itertools
import json
import numpy as np


def check(r):
    X = np.array(list(itertools.product((0,1), repeat=r)), dtype=np.int64)
    N = len(X)
    W = 1-2*((X @ X.T) % 2)
    pairs = list(itertools.combinations(range(r),2))
    P = np.array([X[:,i]*X[:,j] for i,j in pairs],dtype=np.int64)
    real = np.array([1,0,-1,0],dtype=np.int64)
    imag = np.array([0,1,0,-1],dtype=np.int64)
    maximum = -1
    count = 0
    for alpha in itertools.product(range(4),repeat=r):
        linear = X @ np.array(alpha,dtype=np.int64)
        for beta in itertools.product((0,1),repeat=len(pairs)):
            phase = (linear + 2*np.array(beta,dtype=np.int64) @ P) % 4
            value = int(np.abs(W @ real[phase]).sum()+np.abs(W @ imag[phase]).sum())
            maximum = max(maximum,value)
            count += 1
    bound = N*(2**(r//2))
    assert maximum <= bound
    return {'r':r,'phases_checked_modulo_global_fourth_root':count,
            'maximum_integer_statistic':maximum,'proved_bound':bound,
            'maximum_H2_seed_value_in_class':maximum/bound}


if __name__ == '__main__':
    print(json.dumps({'checks':[check(2),check(4)]}))
