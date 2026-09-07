"""Exact exhaustive principal-block and finite response checks; no asymptotics."""
import itertools
import json
from pathlib import Path
import numpy as np


def main():
    xs = np.array(list(itertools.product((-1, 1), repeat=4)), dtype=np.int64)
    even = xs[np.prod(xs, axis=1) == 1]
    odd = xs[np.prod(xs, axis=1) == -1]
    rr = np.ones((4, 4), dtype=np.int64) - 2*np.eye(4, dtype=np.int64)
    blocks = []
    for d in (-1, 1):
        for bits in itertools.product((-1, 1), repeat=6):
            s = d*np.eye(4, dtype=np.int64)
            for (i, j), val in zip(itertools.combinations(range(4), 2), bits):
                s[i, j] = s[j, i] = val
            if not np.all(np.prod(s, axis=0) == 1):
                continue
            transformed = rr @ s @ rr
            if not np.all(np.abs(transformed) == 4):
                continue
            v = d*s[0]
            assert np.prod(v) == 1
            assert np.array_equal(s, d*np.outer(v, v))
            blocks.append((d, s))
    assert len(blocks) == 8
    checks = 0
    # Every h below is a realizable sum of even-parity exterior columns.
    for d, s in blocks:
        for coeff in itertools.product((-1, 0, 1), repeat=4):
            h = np.array(coeff) @ even[:4]
            for e in (-5, 0, 7):
                old = e + xs @ h + (np.einsum('bi,ij,bj->b', xs, s, xs)-np.trace(s))//2
                ys2 = xs @ rr
                assert np.all(ys2 % 2 == 0)
                ys = ys2//2
                new = e + ys @ h + (np.einsum('bi,ij,bj->b', ys, s, ys)-np.trace(s))//2
                common = np.max(np.abs(e+even@h+(np.einsum('bi,ij,bj->b',even,s,even)-np.trace(s))//2))
                assert int(np.max(np.abs(old))) == max(int(common), abs(e)+int(np.max(odd@h)))
                assert int(np.max(np.abs(new))) == max(int(common), abs(e)+2*int(np.max(np.abs(h))))
                checks += 1
    rng = np.random.default_rng(202609070751)
    n = 12
    all_x = 1-2*((np.arange(1 << (n-1))[:, None] >> np.arange(n)) & 1)
    examples = {}
    for attempt in range(100):
        s = np.ones((n, n), dtype=np.int64)
        ext = even[rng.integers(0, len(even), size=n-4)].T
        s[:4, 4:] = ext
        s[4:, :4] = ext.T
        for i, j in itertools.combinations(range(4, n), 2):
            s[i, j] = s[j, i] = 1-2*int(rng.integers(2))
        new_s = s.copy()
        new_s[:4] = (rr @ new_s[:4])//2
        new_s[:, :4] = (new_s[:, :4] @ rr)//2
        assert np.all(np.abs(new_s) == 1)
        def cap(a):
            return int(np.max(np.abs((np.einsum('bi,ij,bj->b', all_x,a,all_x)-np.trace(a))//2)))
        before, after = cap(s), cap(new_s)
        if before != after:
            examples['change'] = dict(before=before, after=after, matrix=s.tolist(), traded=new_s.tolist())
            # Reverse trade supplies the opposite strict inequality.
            assert np.array_equal(rr @ new_s[:4, 4:], 2*s[:4, 4:])
            break
    assert examples, 'Finite search found no strict cap example'
    result = dict(status='PASS: EXACT INTEGER CHECKS', principal_blocks=len(blocks), response_checks=checks, examples=examples,
                  scope='Local identity only; no universal descent or convergence claim')
    Path('computations/results/flatify_director_closed_quad_response_2026_09_07.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
