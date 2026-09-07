import itertools
import random

import numpy as np


def energies(A):
    n = len(A)
    out = []
    for bits in itertools.product((-1, 1), repeat=n):
        x = np.array(bits, dtype=np.int64)
        out.append((int(x @ A @ x), x))
    return out


def qpn(A):
    es = energies(A)
    p = max(e for e, _ in es)
    mn = min(e for e, _ in es)
    return max(p, -mn), p, -mn, es


def build_tower(A, y0, rng):
    R = list(range(len(A)))
    y = y0.copy()
    blocks = []
    fresh = []
    rows = []
    while len(R) > 1:
        C = A[np.ix_(R, R)]
        q, p, neg, es = qpn(C)
        if p != q:
            break
        grounds = [x for e, x in es if e == p]
        rng.shuffle(grounds)
        chosen = None
        for x in grounds:
            # Pick the projective representative with a nonempty agreement core.
            for z in (x, -x):
                D = np.flatnonzero(z != y)
                E = np.flatnonzero(z == y)
                if len(D) and len(E):
                    chosen = (z.copy(), D, E)
                    break
            if chosen is not None:
                break
        if chosen is None:
            break
        x, Di, Ei = chosen
        Dglob = [R[i] for i in Di]
        Eglob = [R[i] for i in Ei]
        DD = C[np.ix_(Di, Di)]
        EE = C[np.ix_(Ei, Ei)]
        B = C[np.ix_(Di, Ei)]
        HD = int(x[Di] @ DD @ x[Di])
        HE = int(x[Ei] @ EE @ x[Ei])
        qD = qpn(DD)[0] if len(Di) else 0
        qE = qpn(EE)[0] if len(Ei) else 0
        inherited_energy = int(y @ C @ y)
        g = q - inherited_energy
        c = int(x[Di] @ B @ x[Ei])
        LF = int(np.abs(np.diag(x[Di]) @ B @ x[Ei]).sum())
        gammaD = qD - HD
        gammaE = qE - HE
        zeta = qD + qE - q
        f = max(2 * LF - gammaE, 0)
        h = max(g // 2 - zeta - f, 0)
        h2 = max(g - gammaD - max(gammaE, 2 * LF), 0)
        assert g == 4 * c and h == h2
        assert h <= 2 * c
        rows.append((g, c, h, HD, qD, HE, qE, LF))
        blocks.append(Dglob)
        fresh.append((Dglob, x[Di].copy()))
        R = Eglob
        y = x[Ei].copy()
    blocks.append(R)
    terminal_y = y.copy()
    return rows, blocks, fresh, terminal_y


def audit(A, y0, rng):
    q0, p0, n0, _ = qpn(A)
    if p0 != q0:
        A = -A
        q0, p0, n0, _ = qpn(A)
    rows, blocks, fresh, terminal_y = build_tower(A, y0, rng)
    if not rows:
        return
    s = np.zeros(len(A), dtype=np.int64)
    Hsum = 0
    for Dglob, sx in fresh:
        s[Dglob] = sx
        AD = A[np.ix_(Dglob, Dglob)]
        Hsum += int(sx @ AD @ sx)
    terminal = blocks[-1]
    s[terminal] = -terminal_y
    AT = A[np.ix_(terminal, terminal)]
    Hsum += int(terminal_y @ AT @ terminal_y)
    csum = sum(r[1] for r in rows)
    hsum = sum(r[2] for r in rows)
    assert int(s @ A @ s) == Hsum - 2 * csum
    assert Hsum <= p0
    assert -int(s @ A @ s) <= n0
    assert hsum <= 2 * csum <= p0 + n0 <= 2 * q0


def main():
    rng = random.Random(20260729)
    count = 0
    for n in range(3, 10):
        for _ in range(300):
            upper = [rng.choice((-1, 1)) for _ in range(n * (n - 1) // 2)]
            A = np.zeros((n, n), dtype=np.int64)
            k = 0
            for i in range(n):
                for j in range(i + 1, n):
                    A[i, j] = A[j, i] = upper[k]
                    k += 1
            y0 = np.array([rng.choice((-1, 1)) for _ in range(n)], dtype=np.int64)
            audit(A, y0, rng)
            count += 1
    print(f"audited {count} random signings/tower starts")


if __name__ == "__main__":
    main()
