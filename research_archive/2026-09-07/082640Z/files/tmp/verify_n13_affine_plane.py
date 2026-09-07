"""Exact verifier for the prime-order affine-plane counterexample."""
A = [
[0,1,-1,1,1,-1,-1,1,1,-1,1,1,-1],
[1,0,1,1,1,1,-1,1,1,1,1,1,1],
[-1,1,0,1,1,1,1,1,1,1,1,-1,1],
[1,1,1,0,-1,1,-1,1,1,-1,1,1,1],
[1,1,1,-1,0,1,1,1,1,1,1,1,-1],
[-1,1,1,1,1,0,1,1,1,1,-1,1,1],
[-1,-1,1,-1,1,1,0,1,-1,1,1,-1,-1],
[1,1,1,1,1,1,1,0,-1,-1,-1,-1,-1],
[1,1,1,1,1,1,-1,-1,0,-1,-1,-1,-1],
[-1,1,1,-1,1,1,1,-1,-1,0,-1,1,1],
[1,1,1,1,1,-1,1,-1,-1,-1,0,-1,1],
[1,1,-1,1,1,1,-1,-1,-1,1,-1,0,-1],
[-1,1,1,1,-1,1,-1,-1,-1,1,1,-1,0],
]
n=len(A)
types=[0]*6+[1]+[2]*6
w=2
def dot2(a,b): return (a&1)*(b&1)+((a>>1)&1)*((b>>1)&1)
def q(h,i,j): return -1 if dot2(types[i]^types[j]^w,h)%2 else 1
scores=[]; fields=[]
for h in range(4):
    scores.append(sum(A[i][j]*q(h,i,j) for i in range(n) for j in range(i+1,n)))
    fields.append([sum(A[i][j]*q(h,i,j) for j in range(n) if j!=i) for i in range(n)])
energies=[]
child_norms=[]
for mask in range(1<<(n-1)):
    x=[1]+[-1 if mask>>(i-1)&1 else 1 for i in range(1,n)]
    energies.append(sum(A[i][j]*x[i]*x[j] for i in range(n) for j in range(i+1,n)))
for omit in range(n):
    verts=[i for i in range(n) if i!=omit]
    vals=[]
    for mask in range(1<<(len(verts)-1)):
        x=[1]+[-1 if mask>>(k-1)&1 else 1 for k in range(1,len(verts))]
        vals.append(sum(A[verts[a]][verts[b]]*x[a]*x[b]
                        for a in range(len(verts)) for b in range(a+1,len(verts))))
    child_norms.append(max(abs(v) for v in vals))
print('scores',scores)
print('fields')
for row in fields: print(row)
print('zero cover',[any(fields[h][i]==0 for h in range(4)) for i in range(n)])
print('energy range',min(energies),max(energies),'absolute norm',max(map(abs,energies)))
print('child norms',child_norms)
assert scores==[24]*4
assert all(min(row)>=0 for row in fields)
assert all(any(fields[h][i]==0 for h in range(4)) for i in range(n))
assert max(map(abs,energies))==24
assert child_norms==[24]*n
