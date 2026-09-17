from itertools import combinations

A = [
 [0,1,-1,1,1,1,1,-1,-1],
 [1,0,1,-1,-1,-1,1,-1,-1],
 [-1,1,0,1,1,1,1,1,-1],
 [1,-1,1,0,1,1,1,-1,1],
 [1,-1,1,1,0,1,-1,1,-1],
 [1,-1,1,1,1,0,-1,-1,-1],
 [1,1,1,1,-1,-1,0,1,1],
 [-1,-1,1,-1,1,-1,1,0,-1],
 [-1,-1,-1,1,-1,-1,1,-1,0],
]
n=len(A)

def decode(v):
    # coordinates: bit 0 is t=-1; bits 1..n-1 are x_j=-1, x_0=1.
    t = -1 if (v & 1) else 1
    x=[1]+[-1 if ((v >> j) & 1) else 1 for j in range(1,n)]
    return t,x

def score(v):
    t,x=decode(v)
    return sum(A[i][j]*t*x[i]*x[j] for i in range(n) for j in range(i+1,n))

E=[v for v in range(1<<n) if score(v)==12]
assert len(E)==25

def fields(v):
    t,x=decode(v)
    return [sum(A[i][j]*t*x[i]*x[j] for j in range(n) if j != i) for i in range(n)]

print("endpoint count",len(E))
print("tight cover", [sum(fields(v)[i]==0 for v in E) for i in range(n)])
print("field profiles")
for v in E:
    print(v, fields(v))

# Enumerate every affine plane and greedily extend it to all affine subspaces.
spaces=set()
for a in E:
    D={a^v for v in E}
    # Store (base, frozenset vector-subspace), canonicalized by point set later.
    def rec(U, candidates):
        pts=frozenset(a^u for u in U)
        spaces.add(pts)
        for w in sorted(candidates):
            if w in U:
                continue
            U2=U | {u^w for u in U}
            if U2 <= D:
                rec(U2, {z for z in D if z>w})
    rec({0},D-{0})

by_size={}
for F in spaces:
    by_size.setdefault(len(F),[]).append(F)
print("affine space counts", {k:len(v) for k,v in sorted(by_size.items())})
for size in sorted(by_size, reverse=True):
    best=[]
    maxcov=-1
    for F in by_size[size]:
        cov={i for i in range(n) if any(fields(v)[i]==0 for v in F)}
        if len(cov)>maxcov:
            maxcov=len(cov); best=[(F,cov)]
        elif len(cov)==maxcov:
            best.append((F,cov))
    print("size",size,"max tight-coordinate cover",maxcov,"examples",len(best))
    for F,cov in best[:3]: print(" ",sorted(F),sorted(cov))
