#!/usr/bin/env python3
"""Joint conserved-allocation / path-Carleson minimax experiments.

All arithmetic used to build trees is integral.  The LP is solved by scipy
HiGHS; solutions can then be rationalized and checked independently.
Scratch file for Wave 8; never commit.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from itertools import product
from functools import lru_cache
import numpy as np
from scipy.optimize import linprog


def spins(n):
    if n == 0:
        return ((),)
    return tuple(x + (1,) for x in product((-1, 1), repeat=n - 1))


def energy(a, x):
    return sum(a[i][j] * x[i] * x[j] for i in range(len(a)) for j in range(len(a)))


@lru_cache(None)
def extrema(a):
    vals = tuple((energy(a, x), x) for x in spins(len(a)))
    P = max((z for z, _ in vals), default=0)
    lo = min((z for z, _ in vals), default=0)
    return P, -lo, tuple(x for z, x in vals if z == P), tuple(x for z, x in vals if z == lo)


def induced(a, ids):
    return tuple(tuple(a[i][j] for j in ids) for i in ids)


@dataclass
class Node:
    a: tuple
    ids: tuple
    p: tuple | None = None
    n: tuple | None = None
    P: int = 0
    N: int = 0
    d: int = 0
    children: list["Node"] = field(default_factory=list)
    # For each child: capacities (plus, minus), inherited h, cross L.
    caps: list[tuple[int, int]] = field(default_factory=list)
    hL: list[tuple[int, int]] = field(default_factory=list)
    label: int = -1

    @property
    def R(self):
        return self.P + self.N


def make_tree(a, ids=None, chooser=None, depth=0):
    """Build a fixed endpoint tree.

    chooser(a,P,N,positive_endpoints,negative_endpoints,depth) -> (p,n).
    Default is lexicographically first pair.
    """
    if ids is None:
        ids = tuple(range(len(a)))
    P, N, ps, ns = extrema(a)
    v = Node(a=a, ids=ids, P=P, N=N, d=abs(P-N))
    if len(a) <= 1 or P + N == 0:
        return v
    p, n = chooser(a, P, N, ps, ns, depth) if chooser else (ps[0], ns[0])
    ratio = tuple(x*y for x, y in zip(p, n))
    sides = [tuple(i for i,z in enumerate(ratio) if z == s) for s in (-1,1)]
    assert all(sides), (len(a), P, N, p, n)
    v.p, v.n = p, n
    for side in sides:
        other = tuple(i for i in range(len(a)) if i not in side)
        child_a = induced(a, side)
        child_ids = tuple(ids[i] for i in side)
        child = make_tree(child_a, child_ids, chooser, depth+1)
        h = energy(child_a, tuple(p[i] for i in side))
        L = sum(abs(sum(a[j][i]*p[i] for i in side)) for j in other)
        Q = max(child.P, child.N)
        caps = (max(0, 2*L-(Q-h)), max(0, 2*L-(Q+h)))
        v.children.append(child)
        v.caps.append(caps)
        v.hL.append((h,L))
    return v


def relabel(root):
    nodes=[]
    def go(v):
        v.label=len(nodes); nodes.append(v)
        for x in v.children: go(x)
    go(root)
    return nodes


def primal_lp(root, W=None):
    """Minimize theta_root using the recursive formulation (10.412)."""
    nodes=relabel(root)
    if W is None: W=root.R
    # Variables: a(v,j,s) for positive cap buckets; w(nonroot); theta(v);
    # z(v,j) for each tree edge.  Zero caps have no a variable.
    names=[]; idx={}
    def add(name): idx[name]=len(names); names.append(name); return idx[name]
    for v in nodes:
        for j,caps in enumerate(v.caps):
            for s,c in enumerate(caps):
                if c>0: add(("a",v.label,j,s))
    for v in nodes[1:]: add(("w",v.label))
    for v in nodes: add(("th",v.label))
    for v in nodes:
        for j in range(len(v.children)): add(("z",v.label,j))
    nv=len(names)
    cobj=np.zeros(nv); cobj[idx[("th",root.label)]]=1
    Aub=[]; bub=[]; Aeq=[]; beq=[]
    # Conservation.
    for v in nodes:
        row=np.zeros(nv)
        for j,caps in enumerate(v.caps):
            for s,c in enumerate(caps):
                if c>0: row[idx[("a",v.label,j,s)]] += 1
            row[idx[("w",v.children[j].label)]] += 1
        if v is root:
            Aeq.append(row); beq.append(W)
        elif v.children:
            row[idx[("w",v.label)]] -= 1
            Aeq.append(row); beq.append(0)
        else:
            # leaf conservation forces incoming obligation to zero
            row[idx[("w",v.label)]] -= 1
            Aeq.append(row); beq.append(0)
    # z >= load and z >= theta_child; theta_v >= sum z.
    for v in nodes:
        if not v.children:
            row=np.zeros(nv); row[idx[("th",v.label)]]=1
            Aeq.append(row); beq.append(0)
            continue
        for j,ch in enumerate(v.children):
            row=np.zeros(nv)
            ze=idx[("z",v.label,j)]
            row[ze]-=1
            for s,cap in enumerate(v.caps[j]):
                if cap>0: row[idx[("a",v.label,j,s)]] += 1/cap
            Aub.append(row); bub.append(0)
            row=np.zeros(nv); row[idx[("th",ch.label)]] += 1; row[ze]-=1
            Aub.append(row); bub.append(0)
        row=np.zeros(nv); row[idx[("th",v.label)]]=-1
        for j in range(len(v.children)): row[idx[("z",v.label,j)]] += 1
        Aub.append(row); bub.append(0)
    bounds=[]
    for name in names:
        if name[0]=="a":
            _,vi,j,s=name; bounds.append((0,nodes[vi].caps[j][s]))
        elif name[0]=="w":
            bounds.append((0,nodes[name[1]].d))
        else: bounds.append((0,None))
    ans=linprog(cobj,A_ub=np.array(Aub),b_ub=np.array(bub),A_eq=np.array(Aeq),b_eq=np.array(beq),bounds=bounds,method="highs")
    if not ans.success: raise RuntimeError(ans.message)
    return ans, names, nodes


def minimax_dual_lp(root,W=None):
    """Dual/minimax LP: antichain distribution plus min-cost-flow dual."""
    nodes=relabel(root)
    if W is None:W=root.R
    # Identify tree edges by (parent label, child index).
    edges=[]; parent_edge={}
    for v in nodes:
        for j,ch in enumerate(v.children):
            e=(v.label,j);edges.append(e);parent_edge[ch.label]=e
    names=[];idx={}
    def add(z):idx[z]=len(names);names.append(z)
    for e in edges:add(("q",)+e)
    for v in nodes:add(("y",v.label))
    for v in nodes:
        for j,caps in enumerate(v.caps):
            for s,c in enumerate(caps):
                if c>0:add(("z",v.label,j,s))
    for v in nodes[1:]:add(("s",v.label))
    nv=len(names);cobj=np.zeros(nv);cobj[idx[("y",root.label)]]=-W
    for v in nodes:
        for j,caps in enumerate(v.caps):
            for s,c in enumerate(caps):
                if c>0:cobj[idx[("z",v.label,j,s)]]=c
    for v in nodes[1:]:cobj[idx[("s",v.label)]]=v.d
    Aub=[];bub=[]
    # q is in the edge-antichain (chain) polytope: every root-leaf path <=1.
    def paths(v,pref):
        if not v.children:return [pref]
        out=[]
        for j,ch in enumerate(v.children):out += paths(ch,pref+[(v.label,j)])
        return out
    for path in paths(root,[]):
        row=np.zeros(nv)
        for e in path:row[idx[("q",)+e]]=1
        Aub.append(row);bub.append(1)
    # Bucket reduced-cost constraints y_v-z_b <= q_e/c_b.
    for v in nodes:
        for j,caps in enumerate(v.caps):
            for s,c in enumerate(caps):
                if c<=0:continue
                row=np.zeros(nv);row[idx[("y",v.label)]]=1;row[idx[("z",v.label,j,s)]]=-1;row[idx[("q",v.label,j)]]=-1/c
                Aub.append(row);bub.append(0)
    # Child-flow reduced-cost constraints y_parent-y_child-s_child <=0.
    for v in nodes[1:]:
        pe=parent_edge[v.label];par=nodes[pe[0]]
        row=np.zeros(nv);row[idx[("y",par.label)]]=1;row[idx[("y",v.label)]]=-1;row[idx[("s",v.label)]]=-1
        Aub.append(row);bub.append(0)
    bounds=[]
    for name in names:
        bounds.append((None,None) if name[0]=="y" else (0,None))
    ans=linprog(cobj,A_ub=np.array(Aub),b_ub=np.array(bub),bounds=bounds,method="highs")
    if not ans.success:raise RuntimeError(ans.message)
    return ans,names,nodes


def audit_primal_dual(root,W=None):
    p,pnames,nodes=primal_lp(root,W);d,dnames,_=minimax_dual_lp(root,W)
    print("PD",frac(p.fun,10**9),frac(-d.fun,10**9),"gap",p.fun+d.fun)
    assert abs(p.fun+d.fun)<1e-7
    print("DUAL NONZERO")
    for x,nm in zip(d.x,dnames):
        if abs(x)>1e-8:print(nm,frac(x,10**7))
    return p,d


def exact_primal_certificate(root,W=None,den=10**7):
    """Rationalize a HiGHS primal and verify conservation/K exactly."""
    if W is None:W=root.R
    ans,names,nodes=primal_lp(root,W);xx={nm:Fraction(float(x)).limit_denominator(den) for nm,x in zip(names,ans.x)}
    for v in nodes:
        incoming=Fraction(W) if v is root else xx[("w",v.label)]
        outgoing=Fraction()
        for j,ch in enumerate(v.children):
            outgoing += xx[("w",ch.label)]
            for s,c in enumerate(v.caps[j]):
                if c:
                    a=xx[("a",v.label,j,s)];assert 0<=a<=c
                    outgoing += a
        assert incoming==outgoing,(v.label,incoming,outgoing)
        if v is not root:assert incoming<=v.d
    def theta(v):
        if not v.children:return Fraction()
        total=Fraction()
        for j,ch in enumerate(v.children):
            load=Fraction()
            for s,c in enumerate(v.caps[j]):
                if c:load += xx[("a",v.label,j,s)]/c
            total += max(load,theta(ch))
        return total
    kval=theta(root);assert kval==Fraction(float(ans.fun)).limit_denominator(den),(kval,ans.fun)
    return kval,xx


def exact_dual_certificate(root,W=None,den=10**7):
    """Rationalize a HiGHS minimax dual and check every inequality exactly."""
    if W is None:W=root.R
    ans,names,nodes=minimax_dual_lp(root,W);xx={nm:Fraction(float(x)).limit_denominator(den) for nm,x in zip(names,ans.x)}
    def paths(v,pref):
        if not v.children:return [pref]
        return sum((paths(ch,pref+[(v.label,j)]) for j,ch in enumerate(v.children)),[])
    for path in paths(root,[]):assert sum(xx[("q",)+e] for e in path)<=1
    parent={ch.label:v for v in nodes for ch in v.children}
    for v in nodes:
        for j,caps in enumerate(v.caps):
            q=xx[("q",v.label,j)]
            for s,c in enumerate(caps):
                if c:assert xx[("y",v.label)]-xx[("z",v.label,j,s)]<=q/c
    for v in nodes[1:]:assert xx[("y",parent[v.label].label)]-xx[("y",v.label)]-xx[("s",v.label)]<=0
    obj=Fraction(W)*xx[("y",root.label)]
    for v in nodes:
        for j,caps in enumerate(v.caps):
            for s,c in enumerate(caps):
                if c:obj-=c*xx[("z",v.label,j,s)]
    for v in nodes[1:]:obj-=v.d*xx[("s",v.label)]
    assert obj==Fraction(float(-ans.fun)).limit_denominator(den),(obj,ans.fun)
    return obj,xx


def frac(x, lim=100000): return Fraction(float(x)).limit_denominator(lim)


def summarize(root,W=None):
    ans,names,nodes=primal_lp(root,W)
    print("LP", "n",len(root.a),"P,N,R,d",root.P,root.N,root.R,root.d,"W",root.R if W is None else W,"K",frac(ans.fun),float(ans.fun),"nodes",len(nodes))
    for val,name in zip(ans.x,names):
        if val>1e-8 and name[0] in ("a","w","th","z"):
            print(name,frac(val))
    print("TREE")
    for v in nodes:
        print(v.label,"ids",v.ids,"P,N,d",v.P,v.N,v.d,"caps",v.caps,"hL",v.hL,"p",v.p,"n",v.n)
    return ans,names,nodes


def optimize_root_pairs(a, W=None):
    """Exact enumeration of root endpoint ties; descendants lexicographic."""
    P,N,ps,ns=extrema(a)
    rows=[]
    for ip,p in enumerate(ps):
        for inn,n in enumerate(ns):
            def choose(_a,_P,_N,_ps,_ns,depth,p=p,n=n):
                return (p,n) if depth==0 else (_ps[0],_ns[0])
            root=make_tree(a,chooser=choose)
            ans,_,_=primal_lp(root,W)
            rows.append((Fraction(float(ans.fun)).limit_denominator(100000),ip,inn,p,n,root))
    rows.sort(key=lambda z:(z[0],z[1],z[2]))
    for z in rows[:10]:
        print("ROOTOPT",z[0],"indices",z[1:3],"p",z[3],"n",z[4],"split",[x.ids for x in z[5].children],"caps",z[5].caps,"d",[x.d for x in z[5].children])
    return rows


def all_sign(n,s=1):
    return tuple(tuple(0 if i==j else s for j in range(n)) for i in range(n))


def matrix_code(n,code):
    a=[[0]*n for _ in range(n)]; k=0
    for i in range(n):
        for j in range(i+1,n):
            a[i][j]=a[j][i]=1 if (code>>k)&1 else -1; k+=1
    return tuple(tuple(r) for r in a)


@lru_cache(None)
def tree_signatures(a):
    """All LP-distinct fixed endpoint-tree signatures, with tie choices."""
    P,N,ps,ns=extrema(a); d=abs(P-N)
    if len(a)<=1 or P+N==0:
        return ((d,()),)
    out=set()
    for p in ps:
        for n in ns:
            ratio=tuple(x*y for x,y in zip(p,n))
            sides=[tuple(i for i,z in enumerate(ratio) if z==s) for s in (-1,1)]
            if not all(sides): continue
            choices=[]
            for side in sides:
                other=tuple(i for i in range(len(a)) if i not in side)
                ca=induced(a,side); Pc,Nc,*_=extrema(ca); Q=max(Pc,Nc)
                h=energy(ca,tuple(p[i] for i in side))
                L=sum(abs(sum(a[j][i]*p[i] for i in side)) for j in other)
                caps=tuple(sorted((max(0,2*L-(Q-h)),max(0,2*L-(Q+h)))))
                choices.append(tuple((caps,sig) for sig in tree_signatures(ca)))
            for b0 in choices[0]:
                for b1 in choices[1]:
                    out.add((d,tuple(sorted((b0,b1),key=repr))))
    return tuple(sorted(out,key=repr))


def node_from_signature(sig):
    d,branches=sig
    v=Node(a=tuple(),ids=tuple(),P=d,N=0,d=d)
    for caps,csig in branches:
        v.caps.append(caps); v.hL.append((0,0)); v.children.append(node_from_signature(csig))
    return v


def optimize_all_trees(a,W=None,limit_print=10):
    P,N,*_=extrema(a); W=P+N if W is None else W
    sigs=tree_signatures(a)
    print("SIGNATURES",len(a),len(sigs))
    rows=[]
    for sig in sigs:
        root=node_from_signature(sig); root.P=P; root.N=N; root.d=abs(P-N)
        ans,_,_=primal_lp(root,W)
        rows.append((frac(ans.fun),sig,root,ans))
    rows.sort(key=lambda x:(x[0],repr(x[1])))
    for k,(val,sig,root,ans) in enumerate(rows[:limit_print]):
        print("ALLOPT",k,val,"rootsig",sig[:1],"branches",sig[1])
    return rows


A37=matrix_code(6,37)

RESET12=(
 (0,-1,-1,-1,-1,-1,-1,1,-1,1,-1,-1),
 (-1,0,-1,-1,-1,-1,1,-1,1,-1,-1,-1),
 (-1,-1,0,1,-1,1,1,-1,1,1,1,-1),
 (-1,-1,1,0,-1,1,1,1,1,1,1,1),
 (-1,-1,-1,-1,0,1,1,1,1,-1,1,1),
 (-1,-1,1,1,1,0,1,-1,-1,-1,-1,-1),
 (-1,1,1,1,1,1,0,-1,1,-1,-1,-1),
 (1,-1,-1,1,1,-1,-1,0,-1,1,1,1),
 (-1,1,1,1,1,-1,1,-1,0,1,1,1),
 (1,-1,1,1,-1,-1,-1,1,1,0,1,-1),
 (-1,-1,1,1,1,-1,-1,1,1,1,0,-1),
 (-1,-1,-1,1,1,-1,-1,1,1,-1,-1,0),
)


def family_A(m):
    k=int(round(m**.5)); assert k*k==m
    s=(1,)*((m+k)//2)+(-1,)*((m-k)//2)
    C=tuple(tuple(0 if i==j else s[i]*s[j] for j in range(m)) for i in range(m))
    return tuple(tuple(C[i][j] if i<m and j<m else -C[i-m][j-m] if i>=m and j>=m else 1 for j in range(2*m)) for i in range(2*m))


if __name__=="__main__":
    import argparse
    ap=argparse.ArgumentParser(); ap.add_argument("which",choices=("triangle","a37","am4","reset")); ap.add_argument("--W",type=int); ap.add_argument("--root-opt",action="store_true"); ap.add_argument("--all-trees",action="store_true")
    z=ap.parse_args()
    a={"triangle":all_sign(3),"a37":A37,"am4":family_A(4),"reset":RESET12}[z.which]
    if z.all_trees: optimize_all_trees(a,z.W)
    elif z.root_opt: optimize_root_pairs(a,z.W)
    else: summarize(make_tree(a),z.W)
