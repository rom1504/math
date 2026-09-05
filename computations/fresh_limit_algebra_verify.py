"""Independent checks for the fixed outers and stored finite witnesses."""
import itertools
import json
from pathlib import Path
import numpy as np


def verify_outer():
    residues = {i*i % 11 for i in range(1, 11)}
    chi = lambda x: 0 if x % 11 == 0 else (1 if x % 11 in residues else -1)
    C = np.array([[chi(i-j) for j in range(11)] for i in range(11)], dtype=np.int64)
    F = np.block([[np.ones((1, 12), dtype=np.int64)],
                  [np.column_stack((-np.ones(11, dtype=np.int64), np.eye(11, dtype=np.int64)+C))]])
    assert np.array_equal(F @ F.T, 12*np.eye(12, dtype=np.int64))
    K = np.array([[F[i,q]*F[p,j] for p in range(12) for q in range(12)]
                  for i in range(12) for j in range(12)])
    d = F.reshape(-1)
    H = K*d[:, None]*d[None, :]
    assert np.array_equal(H, H.T)
    assert np.array_equal(H @ H, 144*np.eye(144, dtype=np.int64))
    assert np.all(H.sum(axis=1) == 12)
    assert np.all(np.diag(H) == 1)
    return {'F_order': 12, 'H_order': 144, 'H_row_sum': 12,
            'H_symmetric': True, 'H_Hadamard': True}


def seed(n):
    if n == 2:
        return np.array([[1,1],[1,-1]], dtype=np.int64)
    if n == 5:
        return np.array([[0,-1,1,-1,1],[-1,0,-1,1,1],[1,-1,0,1,-1],
                         [-1,1,1,0,-1],[1,1,-1,-1,0]], dtype=np.int64)
    A = np.ones((n,n), dtype=np.int64)-np.eye(n,dtype=np.int64)
    for bit, (i,j) in enumerate(itertools.combinations(range(1,n),2)):
        A[i,j] = A[j,i] = 1-2*((220>>bit)&1)
    return A


def check_witnesses():
    results = Path(__file__).parent/'results'
    runs = []
    for filename in ('fresh_limit_algebra_witnesses.json', 'fresh_limit_algebra_h2_witnesses.json'):
        runs.extend(json.loads((results/filename).read_text())['runs'])
    count = 0
    for run in runs:
        H = np.ones((1,1), dtype=np.int64)
        H4 = np.ones((4,4), dtype=np.int64)-2*np.eye(4,dtype=np.int64)
        while len(H) < run['k']:
            H = np.kron(H,H4)
        A = np.kron(H,seed(run['n']))
        for witness in run['witnesses']:
            x = np.array(witness['spins'],dtype=np.int64)
            assert len(x) == len(A) and np.all(np.abs(x)==1)
            assert int(x @ A @ x)//2 == witness['independently_checked_energy']
            count += 1
    return count


if __name__ == '__main__':
    print(json.dumps({'outer': verify_outer(), 'verified_witnesses': check_witnesses()}))
