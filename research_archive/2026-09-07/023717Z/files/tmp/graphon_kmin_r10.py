#!/usr/bin/env python3
"""Limiting allocation/path minimax for the strict nine-type witness."""

from fractions import Fraction
import sys

sys.path.insert(0, "/home/math/quadra/tmp")
import minimax_allocation_r8 as mm


M = (
    (0,54,54,54,28,28,28,55,29),
    (54,0,54,54,28,28,28,55,29),
    (54,54,0,54,28,28,28,55,29),
    (54,54,54,0,28,28,28,55,29),
    (28,28,28,28,0,83,83,-55,-55),
    (28,28,28,28,83,0,83,-55,-55),
    (28,28,28,28,83,83,0,-55,-55),
    (55,55,55,55,-55,-55,-55,0,113),
    (29,29,29,29,-55,-55,-55,113,0),
)


def alternate_matrix():
    """6000 times the alternate 49/50 strictification."""
    a=[[0]*9 for _ in range(9)]
    blocks=(range(3),range(3,5),range(5,8),range(8,9))
    def put(I,J,z,same=False):
        for i in I:
            for j in J:
                if i==j: continue
                if same and i>j: continue
                a[i][j]=a[j][i]=z
    U,V,C,D=blocks
    put(U,U,735,True);put(U,V,-490);put(V,V,980,True)
    put(C,C,735,True);put(C,D,245)
    put(U,C,251);put(U,D,251);put(V,C,496);put(V,D,6)
    return tuple(tuple(r) for r in a)

ALT=alternate_matrix()


def main():
    P,N,ps,ns=mm.extrema(M)
    print("root extrema",P,N,len(ps),len(ns),"R",P+N,"d",abs(P-N))
    print("ps",ps,"ns",ns)
    sigs=mm.tree_signatures(M)
    print("signatures",len(sigs))
    rows=[]
    for sig in sigs:
        root=mm.node_from_signature(sig)
        root.P=P; root.N=N
        val=mm.exact_primal_certificate(root,P+N)[0]
        dual=mm.exact_dual_certificate(root,P+N)[0]
        assert val==dual
        rows.append((val,sig,root))
    rows.sort(key=lambda x:(x[0],repr(x[1])))
    print("value histogram")
    for value in sorted({r[0] for r in rows}):
        print(value, sum(r[0]==value for r in rows))
    for j,(value,sig,root) in enumerate(rows[:20]):
        print("row",j,"K",value,"sig",sig)
    value,sig,root=rows[0]
    print("BEST",value)
    p,xx=mm.exact_primal_certificate(root,P+N)
    d,yy=mm.exact_dual_certificate(root,P+N)
    print("PRIMAL")
    for key,v in xx.items():
        if v: print(key,v)
    print("DUAL")
    for key,v in yy.items():
        if v: print(key,v)


def leaf(): return (0, ())


def edge_sig(c):
    """Two-type edge of magnitude c/2; each leaf capacity is c."""
    return (0, (((c, c), leaf()), ((c, c), leaf())))


def f_sig(p, q):
    """Signature of F(t), t=p/q, after multiplying all weights by q."""
    big = 110*q + 166*p
    low = 110*q - 54*p
    if p == 0:
        return (0, (((110*q, 110*q), leaf()), ((110*q, 110*q), leaf())))
    return (220*p, (((big, big), leaf()), ((low, big), edge_sig(110*p))))


def continuum_root(p, q):
    """Macro tree with a flat T-endpoint parameter t=p/q."""
    fs0 = f_sig(p, q)
    fs1 = f_sig(q-p, q)
    tsig = (992*q, (
        ((668*q-552*p, 888*q), fs0),
        ((668*q-552*(q-p), 888*q), fs1),
    ))
    pair = edge_sig(108*q)
    ssig = (432*q, (((216*q,432*q), pair), ((216*q,432*q), pair)))
    sig = (1424*q, (((0,24*q),tsig), ((48*q,1344*q),ssig)))
    root=mm.node_from_signature(sig); root.P=2056*q;root.N=632*q
    return root


def alt_f_sig(p,q):
    big=980*q+1470*p
    low=980*q-490*p
    if p==0:
        return (0, (((980*q,980*q),leaf()),((980*q,980*q),leaf())))
    return (1960*p, (((big,big),leaf()),((low,big),edge_sig(980*p))))


def alt_continuum_root(p,q):
    fs0=alt_f_sig(p,q);fs1=alt_f_sig(q-p,q)
    flat=(8820*q,(
        ((5880*q-4900*p,7840*q),fs0),
        ((5880*q-4900*(q-p),7840*q),fs1),
    ))
    other=(3944*q,(
        ((2940*q,3920*q),edge_sig(490*q)),
        ((980*q,3920*q),edge_sig(1470*q)),
    ))
    other=(3920*q,other[1])
    sig=(12740*q,(((0,240*q),flat),((240*q,12000*q),other)))
    root=mm.node_from_signature(sig);root.P=18370*q;root.N=5630*q
    return root


BASE_SIG = (1424, (
    ((0,24), (992, (
        ((116,888), (220, (
            ((276,276), (0,())),
            ((56,276), (0, (((110,110),(0,())),((110,110),(0,()))))),
        ))),
        ((668,888), (0, (((110,110),(0,())),((110,110),(0,()))))),
    ))),
    ((48,1344), (432, (
        ((216,432),(0,(((108,108),(0,())),((108,108),(0,()))))),
        ((216,432),(0,(((108,108),(0,())),((108,108),(0,()))))),
    ))),
))


def scale_sig(sig,k):
    d,branches=sig
    return (d*k,tuple((tuple(c*k for c in caps),scale_sig(ch,k)) for caps,ch in branches))


def graft(sig,repl):
    d,branches=sig
    if not branches:return repl
    return (d,tuple((caps,graft(ch,repl)) for caps,ch in branches))


def hierarchical_sig(depth):
    """Exploratory leaf graft, top scale 81^(depth-1), one 9-type copy/cell."""
    assert depth>=1
    if depth==1:return BASE_SIG
    top=scale_sig(BASE_SIG,81**(depth-1))
    return graft(top,hierarchical_sig(depth-1))


def hierarchy_scan():
    print("HIERARCHICAL LEAF-GRAFT TEST")
    for depth in range(1,4):
        sig=hierarchical_sig(depth);root=mm.node_from_signature(sig)
        top=81**(depth-1)
        # Identical internal positive endpoints contribute equally to both
        # parent endpoint energies, so they cancel from the parent range.
        W=2688*top
        v=mm.primal_lp(root,W)[0].fun
        print("depth",depth,"W",W,"K",Fraction(float(v)).limit_denominator(10**8),v,"nodes",len(mm.relabel(root)))


def faithful_hierarchy(depth):
    """Hierarchical M substitution, assuming dominant block-uniform extrema."""
    counter=[0]
    def root_at(d):
        if d==1:
            return node_at(M,tuple(range(9)),1,None)
        return node_at(M,tuple(range(9)),d,root_at(d-1).P)
    def node_at(a,ids,d,innerP):
        # A singleton top block exposes the complete next-level copy.
        if len(a)==1:
            if d==1:return mm.Node(a=a,ids=ids,P=0,N=0,d=0)
            return root_at(d-1)
        P0,N0,ps,ns=mm.extrema(a)
        p,n=ps[0],ns[0]
        if d==1:
            scale=1; pin=0
        else:
            scale=9**(2*(d-1));pin=innerP
        k=len(a)
        P=P0*scale+k*pin;N=N0*scale-k*pin
        assert N>=0,(d,ids,N0,scale,k,pin)
        v=mm.Node(a=a,ids=ids,p=p,n=n,P=P,N=N,d=abs(P-N))
        ratio=tuple(x*y for x,y in zip(p,n))
        for sg in (-1,1):
            side=tuple(i for i,z in enumerate(ratio) if z==sg)
            other=tuple(i for i in range(len(a)) if i not in side)
            ca=mm.induced(a,side);cids=tuple(ids[i] for i in side)
            if len(side)==1 and d>1:
                ch=root_at(d-1)
            else:
                ch=node_at(ca,cids,d,pin)
            h0=mm.energy(ca,tuple(p[i] for i in side))
            h=h0*scale+len(side)*pin
            L0=sum(abs(sum(a[j][i]*p[i] for i in side)) for j in other)
            L=L0*scale
            Q=max(ch.P,ch.N)
            caps=(max(0,2*L-(Q-h)),max(0,2*L-(Q+h)))
            v.children.append(ch);v.hL.append((h,L));v.caps.append(caps)
        return v
    return root_at(depth)


def faithful_scan():
    print("FAITHFUL DOMINANT HIERARCHY TEST")
    for d in range(1,4):
        root=faithful_hierarchy(d)
        v=mm.primal_lp(root,root.R)[0].fun
        print("depth",d,"P,N,R",root.P,root.N,root.R,"K",Fraction(float(v)).limit_denominator(10**8),v,"nodes",len(mm.relabel(root)))


def two_level_scaled(p,q):
    """Two-level formal hierarchy: outer aggregate scale 81q, inner M scale p."""
    outer=81*q;pin=2056*p
    inner_sig=scale_sig(BASE_SIG,p)
    def rec(a,ids):
        if len(a)==1:return mm.node_from_signature(inner_sig)
        P0,N0,ps,ns=mm.extrema(a);P=P0*outer+len(a)*pin;N=N0*outer-len(a)*pin
        if N<0:raise ValueError((ids,N))
        pp,nn=ps[0],ns[0]
        v=mm.Node(a=a,ids=ids,p=pp,n=nn,P=P,N=N,d=abs(P-N))
        ratio=tuple(x*y for x,y in zip(pp,nn))
        for sg in (-1,1):
            side=tuple(i for i,z in enumerate(ratio) if z==sg);other=tuple(i for i in range(len(a)) if i not in side)
            ca=mm.induced(a,side);ch=rec(ca,tuple(ids[i] for i in side))
            h0=mm.energy(ca,tuple(pp[i] for i in side));h=h0*outer+len(side)*pin
            L0=sum(abs(sum(a[j][i]*pp[i] for i in side)) for j in other);L=L0*outer
            Q=max(ch.P,ch.N);caps=(max(0,2*L-(Q-h)),max(0,2*L-(Q+h)))
            v.children.append(ch);v.hL.append((h,L));v.caps.append(caps)
        return v
    return rec(M,tuple(range(9)))


def scaled_hierarchy_scan():
    print("TWO-LEVEL SCALED FORMAL TEST")
    for lam in (Fraction(0),Fraction(1,100),Fraction(1,20),Fraction(1,10),Fraction(1,4),Fraction(1,2),Fraction(1)):
        if lam==0:continue
        r=two_level_scaled(lam.numerator,lam.denominator);v=mm.primal_lp(r,r.R)[0].fun
        print("lambda",lam,"K",Fraction(float(v)).limit_denominator(10**8),v)


def scan_continuum():
    vals=[]
    for q in range(1,101):
        for p in range(q+1):
            if 2*p>q: continue
            root=continuum_root(p,q)
            val=mm.primal_lp(root,2688*q)[0].fun
            vals.append((val,Fraction(p,q),root))
    vals.sort(key=lambda z:(z[0],z[1]))
    print("CONTINUUM BEST")
    for v,t,_ in vals[:20]:print(t,Fraction(float(v)).limit_denominator(10**6),float(v))
    for t in (Fraction(0),Fraction(1,10),Fraction(1,4),Fraction(1,3),Fraction(2,5),Fraction(1,2)):
        root=continuum_root(t.numerator,t.denominator)
        v=mm.primal_lp(root,2688*t.denominator)[0].fun
        print("t",t,"K",Fraction(float(v)).limit_denominator(10**6),v)


def alternate():
    P,N,ps,ns=mm.extrema(ALT)
    print("ALTERNATE root",P,N,P+N,abs(P-N),len(ps),len(ns),ps,ns)
    sigs=mm.tree_signatures(ALT)
    print("ALTERNATE signatures",len(sigs))
    rows=[]
    for sig in sigs:
        root=mm.node_from_signature(sig);root.P=P;root.N=N
        val=mm.primal_lp(root,P+N)[0].fun
        rows.append((val,sig,root))
    rows.sort(key=lambda z:z[0])
    print("ALTERNATE values")
    for v,sig,_ in rows[:10]:print(Fraction(float(v)).limit_denominator(10**7),v,sig)
    root=rows[0][2]
    print("alt exact primal",mm.exact_primal_certificate(root,P+N)[0])
    try:print("alt exact dual",mm.exact_dual_certificate(root,P+N,den=12636000)[0])
    except AssertionError as exc:print("alt exact dual rationalizer failed",exc)
    vals=[]
    for q in range(1,101):
        for p in range(q//2+1):
            v=mm.primal_lp(alt_continuum_root(p,q),24000*q)[0].fun
            vals.append((v,Fraction(p,q)))
    print("alt continuum min/max",min(vals),max(vals))
    for t in (Fraction(0),Fraction(1,10),Fraction(1,4),Fraction(1,3),Fraction(2,5),Fraction(1,2)):
        v=mm.primal_lp(alt_continuum_root(t.numerator,t.denominator),24000*t.denominator)[0].fun
        print("alt t",t,"K",Fraction(float(v)).limit_denominator(10**8),v)


if __name__ == "__main__":
    main()
    scan_continuum()
    alternate()
    hierarchy_scan()
    faithful_scan()
    scaled_hierarchy_scan()
