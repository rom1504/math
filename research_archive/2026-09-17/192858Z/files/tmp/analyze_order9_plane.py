from itertools import combinations
A=[
[0,1,-1,1,1,1,1,-1,-1],
[1,0,1,-1,-1,-1,1,-1,-1],
[-1,1,0,1,1,1,1,1,-1],
[1,-1,1,0,1,1,1,-1,1],
[1,-1,1,1,0,1,-1,1,-1],
[1,-1,1,1,1,0,-1,-1,-1],
[1,1,1,1,-1,-1,0,1,1],
[-1,-1,1,-1,1,-1,1,0,-1],
[-1,-1,-1,1,-1,-1,1,-1,0]]
n=len(A)
vals={}
for mask in range(1<<(n-1)):
 x=[1]+[(-1 if mask>>(i-1)&1 else 1) for i in range(1,n)]
 vals[mask]=sum(A[i][j]*x[i]*x[j] for i in range(n) for j in range(i+1,n))
M=max(map(abs,vals.values()))
ends=[((1 if e==M else -1),m) for m,e in vals.items() if abs(e)==M]
es=set(ends)
def mul(a,b):return(a[0]*b[0],a[1]^b[1])
def field(q):
 t,m=q;x=[1]+[(-1 if m>>(i-1)&1 else 1) for i in range(1,n)]
 return [sum(A[i][j]*t*x[i]*x[j] for j in range(n) if i!=j) for i in range(n)]
fs={q:field(q) for q in ends}
planes=set()
for q0 in ends:
 rel=[mul(q0,q) for q in ends];rs=set(rel)
 for a,b in combinations([z for z in rel if z!=(1,0)],2):
  ab=mul(a,b)
  if ab in rs:
   P=tuple(sorted((q0,mul(q0,a),mul(q0,b),mul(q0,ab))))
   if all(any(fs[q][i]==0 for q in P) for i in range(n)):planes.add(P)
print(M,len(ends),len(planes))
for P in list(planes)[:10]:
 print('PLANE',P)
 for q in P:print(q,fs[q])
 # take q0 P[0], rel gens next independent; derive orientation w bits and vertex types
 q0=P[0]; rels=[mul(q0,q) for q in P]
 gens=[z for z in rels if z!=(1,0)][:2]
 if mul(gens[0],gens[1]) not in rels:raise ValueError
 w=(0 if gens[0][0]==1 else 1)|(0 if gens[1][0]==1 else 2)
 types=[]
 for i in range(n):
  b1=0 if i==0 or not(gens[0][1]>>(i-1)&1) else 1
  b2=0 if i==0 or not(gens[1][1]>>(i-1)&1) else 1
  types.append(b1|(b2<<1))
 print('w',w,'types',types,'counts',[types.count(a) for a in range(4)])
 # switch A by q0 and display in type order
 tq,mq=q0
 xq=[1]+[(-1 if mq>>(i-1)&1 else 1) for i in range(1,n)]
 As=[[0 if i==j else A[i][j]*tq*xq[i]*xq[j] for j in range(n)] for i in range(n)]
 order=sorted(range(n),key=lambda i:types[i])
 print('order',order)
 for i in order:print([As[i][j] for j in order])
 break
