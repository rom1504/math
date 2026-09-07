from functools import lru_cache
from itertools import product
import random


A6 = (
    (0, -1, 1, -1, -1, -1),
    (-1, 0, 1, -1, 1, 1),
    (1, 1, 0, -1, 1, -1),
    (-1, -1, -1, 0, 1, -1),
    (-1, 1, 1, 1, 0, -1),
    (-1, 1, -1, -1, -1, 0),
)


def energy(A, x):
    return sum(A[i][j] * x[i] * x[j] for i in range(len(A)) for j in range(len(A)))


def submatrix(A, inds):
    return tuple(tuple(A[i][j] for j in inds) for i in inds)


@lru_cache(None)
def extrema(A):
    xs = list(product((-1, 1), repeat=len(A)))
    es = [energy(A, x) for x in xs]
    P = max(es)
    N = -min(es)
    return P, N, tuple(x for x, e in zip(xs, es) if e == P), tuple(x for x, e in zip(xs, es) if e == -N)


def pair_data(A, p, n):
    S = tuple(i for i in range(len(A)) if p[i] * n[i] == -1)
    T = tuple(i for i in range(len(A)) if p[i] * n[i] == 1)
    out = []
    for X, Y in ((S, T), (T, S)):
        AX = submatrix(A, X)
        pX = tuple(p[i] for i in X)
        h = energy(AX, pX)
        P, N, _, _ = extrema(AX)
        Q = max(P, N)
        fields = [sum(A[i][j] * p[j] for j in X) for i in Y]
        L = sum(abs(v) for v in fields)
        muplus = max(0, 2 * L - (Q - h))
        muminus = max(0, 2 * L - (Q + h))
        out.append((X, Y, P, N, h, L, muplus, muminus))
    return out


@lru_cache(None)
def best_path(A):
    if len(A) <= 1:
        return 0, None
    P, N, ps, ns = extrema(A)
    best = (-1, None)
    for p in ps:
        for n in ns:
            # Mod global signs, a constant product gives no genuine cut.
            dat = pair_data(A, p, n)
            if any(len(row[0]) == 0 for row in dat):
                continue
            for row in dat:
                X, _, _, _, _, _, mup, mum = row
                val_child, _ = best_path(submatrix(A, X))
                val = max(mup, mum) + val_child
                if val > best[0]:
                    best = (val, (p, n, row, val_child))
    return best


def fixed_tree(A, p, n, indent=""):
    P, N, *_ = extrema(A)
    print(f"{indent}n={len(A)} P={P} N={N} R={P+N}")
    for row in pair_data(A, p, n):
        X, Y, PX, NX, h, L, mup, mum = row
        bval, choice = best_path(submatrix(A, X))
        print(
            f"{indent} retain={X}, peel={Y}, child(P,N)=({PX},{NX}), h={h}, L={L}, "
            f"mu+= {mup}, mu-= {mum}, best-desc={bval}, path={max(mup,mum)+bval}"
        )


def gauge_signing(n, bits):
    A = [[0] * n for _ in range(n)]
    # Switching fixes every edge incident to vertex zero to +1.
    for j in range(1, n):
        A[0][j] = A[j][0] = 1
    k = 0
    for i in range(1, n):
        for j in range(i + 1, n):
            v = 1 if ((bits >> k) & 1) else -1
            A[i][j] = A[j][i] = v
            k += 1
    return tuple(tuple(row) for row in A)


def search_smallest(max_n=6):
    first_bad = None
    for n in range(2, max_n + 1):
        count = 1 << ((n - 1) * (n - 2) // 2)
        bad = []
        min_ratio = None
        for bits in range(count):
            A = gauge_signing(n, bits)
            P, N, *_ = extrema(A)
            R = P + N
            V, _ = best_path(A)
            ratio = V / R if R else 1
            if min_ratio is None or ratio < min_ratio[0]:
                min_ratio = (ratio, bits, P, N, V, A)
            if V < R:
                bad.append((bits, P, N, V, A))
        print(f"n={n}: gauge signings={count}, obstructions={len(bad)}, min={min_ratio[:5]}")
        if bad:
            print("first obstruction:")
            for row in bad[0][-1]:
                print(row)
            if first_bad is None:
                first_bad = bad[0]
    return first_bad


def audit_half_lemma(max_n=6):
    bad_both = []
    bad_individual = 0
    for n in range(2, max_n + 1):
        count = 1 << ((n - 1) * (n - 2) // 2)
        for bits in range(count):
            A = gauge_signing(n, bits)
            P, N, ps, ns = extrema(A)
            for p in ps:
                for neg in ns:
                    dat = pair_data(A, p, neg)
                    if len(dat) != 2 or any(not row[0] for row in dat):
                        continue
                    ok = []
                    for row in dat:
                        _, _, PX, NX, h, *_ = row
                        this = 2 * abs(h) >= abs(PX - NX)
                        ok.append(this)
                        bad_individual += not this
                    if not any(ok):
                        bad_both.append((n, bits, P, N, p, neg, dat))
                        print("both-half-lemma fails", bad_both[-1])
                        return bad_both
    print("half lemma through", max_n, "bad individual shores", bad_individual, "bad both", len(bad_both))
    return bad_both


def random_half_lemma(trials=2000, max_n=10, seed=12345):
    rng = random.Random(seed)
    for trial in range(trials):
        n = rng.randint(2, max_n)
        A = [[0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                A[i][j] = A[j][i] = rng.choice((-1, 1))
        A = tuple(tuple(r) for r in A)
        P, N, ps, ns = extrema(A)
        p, neg = rng.choice(ps), rng.choice(ns)
        dat = pair_data(A, p, neg)
        if len(dat) != 2 or any(not row[0] for row in dat):
            continue
        oks = []
        for row in dat:
            _, _, PX, NX, h, *_ = row
            oks.append(2 * abs(h) >= abs(PX - NX))
        if not any(oks):
            print("random both half lemma fail", trial, n, P, N, p, neg, dat)
            return A, p, neg, dat
    print("random half lemma passed", trials)
    return None


def random_path_ratios(trials=1000, max_n=11, seed=6789):
    rng = random.Random(seed)
    worst = (10, None)
    for trial in range(trials):
        n = rng.randint(2, max_n)
        A = [[0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                A[i][j] = A[j][i] = rng.choice((-1, 1))
        A = tuple(tuple(r) for r in A)
        P, N, *_ = extrema(A)
        V, trace = best_path(A)
        ratio = V / (P + N)
        if ratio < worst[0]:
            worst = (ratio, n, P, N, V, A, trace)
            print("new random path worst", worst[:5])
        if ratio < 0.5:
            print("HALF FAIL", worst[:5])
    print("random path worst final", worst[:5])
    return worst


if __name__ == "__main__":
    p = (-1, -1, -1, 1, -1, 1)
    n = (-1, -1, 1, -1, 1, 1)
    fixed_tree(A6, p, n)
    print("globally best endpoint-tree path:", best_path(A6)[0])
    print("choice:", best_path(A6)[1])
    search_smallest()
    audit_half_lemma()
    random_half_lemma()
    random_path_ratios()
