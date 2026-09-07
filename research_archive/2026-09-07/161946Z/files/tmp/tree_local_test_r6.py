import itertools
import numpy as np


def spins(n):
    for bits in range(1 << max(0, n - 1)):
        x = np.ones(n, dtype=np.int64)
        for i in range(n - 1):
            if bits >> i & 1:
                x[i] = -1
        yield x


def pn(a):
    vals = [int(x @ a @ x) for x in spins(len(a))]
    return max(vals, default=0), -min(vals, default=0)


def matrices(n):
    es = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for bits in range(1 << len(es)):
        a = np.zeros((n, n), dtype=np.int64)
        for k, (i, j) in enumerate(es):
            a[i, j] = a[j, i] = 1 if bits >> k & 1 else -1
        yield a


def subsets(n):
    # modulo complementation, both shores nonempty
    for bits in range(1, 1 << (n - 1)):
        s = [i for i in range(n) if bits >> i & 1]
        t = [i for i in range(n) if not (bits >> i & 1)]
        yield s, t


def audit(n):
    worst = None
    failures = []
    count = 0
    for a0 in matrices(n):
        p0, n0 = pn(a0)
        q0 = max(p0, n0)
        # Every absolute endpoint normalization is again among matrices;
        # restrict to all-one positive endpoint gauges.
        if int(np.ones(n, dtype=np.int64) @ a0 @ np.ones(n, dtype=np.int64)) != q0:
            continue
        for s, t in subsets(n):
            d = a0[np.ix_(s, s)]
            e = a0[np.ix_(t, t)]
            b = a0[np.ix_(s, t)]
            pd, nd = pn(d)
            pe, ne = pn(e)
            qd, qe = max(pd, nd), max(pe, ne)
            hd, he = int(d.sum()), int(e.sum())
            ld = int(np.abs(b.T @ np.ones(len(s), dtype=np.int64)).sum())
            le = int(np.abs(b @ np.ones(len(t), dtype=np.int64)).sum())
            # lambda_D uses L_D = ||B^T 1_S|| and child D deficit.
            taus_d = ([1] if pd == qd else []) + ([-1] if nd == qd else [])
            taus_e = ([1] if pe == qe else []) + ([-1] if ne == qe else [])
            best = None
            for td, te in itertools.product(taus_d, taus_e):
                gd, ge = qd - td * hd, qe - te * he
                lamd, lame = 2 * ld - gd, 2 * le - ge
                cap = max(0, lamd) + max(0, lame)
                row = (cap, td, te, gd, ge, lamd, lame)
                if best is None or row < best:
                    best = row
            delta = p0 + n0 - pd - nd - pe - ne
            ratio = best[0] / delta if delta else (float('inf') if best[0] else 0)
            record = (ratio, best[0], delta, a0.copy(), s, t, (p0,n0,pd,nd,pe,ne), best)
            if worst is None or ratio > worst[0]:
                worst = record
            if best[0] > 2 * delta:
                failures.append(record)
                if len(failures) >= 3:
                    return count, worst, failures
            count += 1
    return count, worst, failures


for n in range(2, 6):
    count, worst, failures = audit(n)
    print('n', n, 'count', count, 'worst ratio/cap/delta', worst[:3], 'failures', len(failures))
    if failures:
        for z in failures:
            print('FAIL', z[0:3], 'S', z[4], 'pn', z[6], 'best', z[7])
            print(z[3].tolist())
        break
