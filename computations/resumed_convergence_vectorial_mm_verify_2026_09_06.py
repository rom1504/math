"""Exact small-field audit of the vectorial-MM transform identity.

Uses only integer arithmetic in every assertion. It checks all truth tables
on F8, and selected balanced/unbalanced profiles on F16, in quadratic field
extensions. No external finite-field library or implicit trace convention.
"""
import numpy as np


def field_tables(n, modulus):
    size=2**n
    def mul(a,b):
        value=0
        while b:
            if b&1:
                value^=a
            b>>=1
            a<<=1
            if a>=size:
                a^=modulus
        return value
    power=[]
    a=1
    for j in range(size-1):
        assert j==0 or a!=1
        power.append(a)
        a=mul(a,2)
    assert a==1
    table=np.array([[mul(a,b) for b in range(size)] for a in range(size)],dtype=np.int64)
    def powfield(a,e):
        out=1
        while e:
            if e&1:
                out=mul(out,a)
            a=mul(a,a)
            e>>=1
        return out
    return mul,powfield,table,power


def check(m, modulus, all_profiles):
    q=2**m
    n=2*m
    size=2**n
    mul,powfield,table,power=field_tables(n,modulus)
    generator=power[(size-1)//(q-1)]
    subbasis=[powfield(generator,j) for j in range(m)]
    subelements=[]
    for bits in range(q):
        a=0
        for j,basis in enumerate(subbasis):
            if (bits>>j)&1:
                a^=basis
        subelements.append(a)
    assert len(set(subelements))==q
    decode={a:j for j,a in enumerate(subelements)}
    def trace(a,degree):
        out=0
        for _ in range(degree):
            out^=a
            a=mul(a,a)
        return out
    abs_trace=np.array([trace(a,n) for a in range(size)],dtype=np.int64)
    assert set(abs_trace)=={0,1}
    rel_trace=np.array([a^powfield(a,q) for a in range(size)],dtype=np.int64)
    assert set(rel_trace)==set(subelements)
    rel_label=np.array([decode[a] for a in rel_trace],dtype=np.int64)
    small_trace=np.array([trace(a,m) for a in subelements],dtype=np.int64)
    assert set(small_trace)=={0,1}
    H=1-2*abs_trace[table]
    Hq=np.array([[1-2*trace(mul(a,b),m) for b in subelements] for a in subelements],dtype=np.int64)
    assert np.array_equal(H@H.T,size*np.eye(size,dtype=np.int64))
    assert np.array_equal(Hq@Hq.T,q*np.eye(q,dtype=np.int64))
    e=2
    d=(size//2)  # 2*d == 1 modulo size-1.
    assert (e*d)%(size-1)==1
    yd=np.array([powfield(y,d) for y in range(size)])
    ae=np.array([powfield(a,e) for a in range(size)])
    phi=rel_label[table[:,yd]]
    dual=rel_label[table[ae,:]]
    expected_counts=np.full(q,(size*size-size)//q,dtype=np.int64)
    expected_counts[0]+=size
    assert np.array_equal(np.bincount(phi.ravel(),minlength=q),expected_counts)
    assert np.array_equal(np.bincount(dual.ravel(),minlength=q),expected_counts)
    p=(-e)%(q-1)
    perm=np.array([0]+[decode[powfield(a,p)] for a in subelements[1:]])
    kernel=Hq[:,perm]@Hq
    assert np.array_equal(kernel@kernel.T,q*q*np.eye(q,dtype=np.int64))
    if all_profiles:
        codes=np.arange(2**q,dtype=np.uint64)
    else:
        rng=np.random.default_rng(749)
        codes=np.unique(np.concatenate(([0,2**q-1],rng.integers(0,2**q,32))))
    for code in codes:
        f=np.array([1-2*((int(code)>>j)&1) for j in range(q)],dtype=np.int64)
        total=f.sum()
        lhs=H@f[phi]@H.T
        rhs=(size//q)*((kernel@f)[dual]-total)
        rhs[0,0]+=(size*size//q)*total
        assert np.array_equal(lhs,rhs),(m,int(code),np.max(np.abs(lhs-rhs)))
    print('PASS q=%d N=%d p=%d: primitive orbit, both traces, distributions, '%(q,size,p)
          +'and %d exact transforms'%len(codes),flush=True)


if __name__=='__main__':
    check(3,0b1000011,True)
    check(4,0b100011101,False)
