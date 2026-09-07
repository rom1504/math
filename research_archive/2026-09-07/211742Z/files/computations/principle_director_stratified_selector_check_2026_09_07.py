"""Exact Hadamard clipping checks and outward-rounded selector constants.

The clipping theorem has an analytic proof; finite enumeration tests its
normalization. The ternary E bound is a separate directed certificate.
"""
import json
from pathlib import Path
import numpy as np
import mpmath as mp


def sylvester(n):
    h = np.ones((1, 1), dtype=np.int64)
    while len(h) < n:
        h = np.block([[h, h], [h, -h]])
    return h


def paley12():
    q = 11
    residues = {i*i % q for i in range(1, q)}
    h = np.ones((12, 12), dtype=np.int64)
    for i in range(q):
        for j in range(q):
            h[i+1, j+1] = -1 if i == j else (1 if (j-i) % q in residues else -1)
    assert np.array_equal(h@h.T, 12*np.eye(12, dtype=np.int64))
    return h


def check(h, exhaustive=True):
    n = len(h)
    assert np.array_equal(h@h.T, n*np.eye(n, dtype=np.int64))
    assert np.all(h[0] == 1)
    if exhaustive:
        numbers = np.arange(1, (1 << n)-1, dtype=np.int64)
        words = ((numbers[:, None] >> np.arange(n)) & 1).astype(np.int64)
    else:
        rng = np.random.default_rng(202609072045)
        words = rng.integers(0, 2, size=(20000, n), dtype=np.int64)
        words = words[(words.sum(axis=1) > 0) & (words.sum(axis=1) < n)]
        words = np.vstack([words, np.eye(n, dtype=np.int64), 1-np.eye(n, dtype=np.int64)])
    sizes = words.sum(axis=1)
    transforms = words@h.T
    sharp = np.minimum(transforms**2, n-1).sum(axis=1) - (n-1+sizes)
    older = np.minimum((n+1)*transforms**2, n*n).sum(axis=1) - (n*n+n*sizes)
    assert np.all(sharp >= 0) and np.all(older >= 0)
    return dict(order=n, cases=len(words), exhaustive=exhaustive,
                sharp_min_integer_slack=int(sharp.min()),
                older_min_integer_slack=int(older.min()), status='PASS')


def bounds(x):
    return [str(x.a), str(x.b)]


def main():
    results = [check(sylvester(n)) for n in [4, 8, 16]]
    results.append(check(paley12()))
    results.append(check(np.kron(sylvester(2), paley12()), False))
    mp.iv.dps = 65
    p = mp.iv.mpf(23)/24
    t = mp.iv.mpf(97)/20
    epsilon = -p*mp.iv.ln(p)-(1-p)*mp.iv.ln(1-p)-mp.iv.ln(24)/24
    balanced = (t+p*mp.iv.ln(2)-mp.iv.mpf(41)/50+epsilon)/(2*t*mp.iv.sqrt(p))
    rare = mp.iv.sqrt(p)/2
    assert balanced.b < mp.iv.mpf(1)/2
    assert rare.b < mp.iv.mpf(1)/2
    certificate_path = Path('computations/results/principle_director_stratified_balanced_certificate_2026_09_07.json')
    certificate = json.loads(certificate_path.read_text())
    assert certificate['status'] == 'directed interval certificate'
    assert certificate['p'] == '23/24' and certificate['t'] == '97/20'
    assert certificate['E_upper'] == '-41/50'
    sylvester_path = Path('computations/results/principle_director_sylvester_stratified_certificate_2026_09_07.json')
    sylvester_certificate = json.loads(sylvester_path.read_text())
    assert sylvester_certificate['status'] == 'directed interval certificate'
    assert sylvester_certificate['p'] == '31/32' and sylvester_certificate['t'] == '97/20'
    assert sylvester_certificate['E_upper'] == '-4/5'
    p32 = mp.iv.mpf(31)/32
    epsilon32 = -p32*mp.iv.ln(p32)-(1-p32)*mp.iv.ln(1-p32)-mp.iv.ln(32)/32
    balanced32 = (t+p32*mp.iv.ln(2)-mp.iv.mpf(4)/5+epsilon32)/(2*t*mp.iv.sqrt(p32))
    assert balanced32.b < mp.iv.mpf(1)/2
    result = dict(status='PASS', clipping=results, precision=65,
                  selector_entropy_interval=bounds(epsilon),
                  L24_formal_balanced_constant_interval=bounds(balanced),
                  L24_inheritance_status='NOT PROVED: signed output row permutations alone do not justify mixing source laws',
                  limiting_rare_band_coefficient_interval=bounds(rare),
                  E_certificate=str(certificate_path),
                  E_checked_boxes=certificate['checked'],
                  E_accepted_leaves=certificate['accepted'],
                  safe_same_law_order=32,
                  safe_same_law_selector_entropy_interval=bounds(epsilon32),
                  safe_same_law_balanced_cap_interval=bounds(balanced32),
                  safe_same_law_rare_coefficient_interval=bounds(mp.iv.sqrt(p32)/2),
                  safe_E_certificate=str(sylvester_path),
                  safe_E_checked_boxes=sylvester_certificate['checked'],
                  safe_E_accepted_leaves=sylvester_certificate['accepted'],
                  qualification='No full mixed-profile cap or minimax convergence claim')
    output = Path('computations/results/principle_director_stratified_selector_check_2026_09_07.json')
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
