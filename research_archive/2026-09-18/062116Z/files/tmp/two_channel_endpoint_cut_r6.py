import numpy as np


def spins(n):
    for bits in range(1 << max(0, n - 1)):
        x = np.ones(n, dtype=np.int64)
        for i in range(n - 1):
            if bits >> i & 1:
                x[i] = -1
        yield x


def pn_states(a):
    rows = [(int(x @ a @ x), x.copy()) for x in spins(len(a))]
    p = max(v for v, _ in rows)
    n = -min(v for v, _ in rows)
    ps = [x for v, x in rows if v == p]
    ns = [x for v, x in rows if v == -n]
    return p, n, ps, ns


def matrices(n):
    es = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for bits in range(1 << len(es)):
        a = np.zeros((n, n), dtype=np.int64)
        for k, (i, j) in enumerate(es):
            a[i, j] = a[j, i] = 1 if bits >> k & 1 else -1
        yield a


def check(a, pspin, nspin):
    n = len(a)
    # Gauge pspin to one; ratio nspin*pspin defines the endpoint cut.
    c = pspin[:,None] * a * pspin[None,:]
    ratio = pspin * nspin
    s = np.flatnonzero(ratio != ratio[-1])
    t = np.array([i for i in range(n) if i not in set(s)], dtype=int)
    if len(s) == 0 or len(t) == 0:
        return None
    d, e, b = c[np.ix_(s,s)], c[np.ix_(t,t)], c[np.ix_(s,t)]
    pd,nd,_,_=pn_states(d); pe,ne,_,_=pn_states(e)
    hd,he=int(d.sum()),int(e.sum())
    ld=int(np.abs(b.T@np.ones(len(s),dtype=int)).sum())
    le=int(np.abs(b@np.ones(len(t),dtype=int)).sum())
    one=[2*ld-(pd-hd),2*ld-(nd+hd),2*le-(pe-he),2*le-(ne+he)]
    qd,qe=max(pd,nd),max(pe,ne)
    absolute=[2*ld-(qd-hd),2*ld-(qd+hd),2*le-(qe-he),2*le-(qe+he)]
    return (s.tolist(),t.tolist(),(pd,nd,pe,ne,hd,he,ld,le),one,absolute)


for n in range(2,6):
    worst_abs=None; audited=0
    for a in matrices(n):
        p,nn,ps,ns=pn_states(a); r=p+nn
        for x in ps:
            for y in ns:
                z=check(a,x,y)
                if z is None: continue
                onecap=sum(max(0,v) for v in z[3])
                abscap=sum(max(0,v) for v in z[4])
                assert onecap>=r,(n,r,z,onecap)
                ratio=abscap/r if r else 1
                row=(ratio,abscap,r,a.copy(),x.copy(),y.copy(),z,p,nn)
                if worst_abs is None or row[0]<worst_abs[0]: worst_abs=row
                audited+=1
    print('n',n,'audited',audited,'worst absolute ratio/cap/R',worst_abs[:3])
    if worst_abs[0]<1:
        print('P,N',worst_abs[7:9], 'data',worst_abs[6]);print(worst_abs[3].tolist())
