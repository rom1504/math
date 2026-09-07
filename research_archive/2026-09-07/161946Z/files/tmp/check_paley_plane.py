from itertools import combinations

def legendre(a,p):
    a%=p
    if a==0:return 0
    return 1 if pow(a,(p-1)//2,p)==1 else -1

def analyze(p):
    A=[[0 if i==j else legendre(i-j,p) for j in range(p)] for i in range(p)]
    # x represented with x0=1; edge word augmented represented (t, mask)
    words=[]
    vals={}
    for mask in range(1<<(p-1)):
        x=[1]+[1 if ((mask>>(i-1))&1)==0 else -1 for i in range(1,p)]
        E=sum(A[i][j]*x[i]*x[j] for i in range(p) for j in range(i+1,p))
        vals[mask]=E
    M=max(abs(e) for e in vals.values())
    endpoints=[]
    for mask,e in vals.items():
        if abs(e)==M:
            # oriented q=t xx with score M
            t=1 if e==M else -1
            endpoints.append((t,mask))
    eset=set(endpoints)
    # group product: (t,m)*(s,n)=(ts,m xor n)
    def mul(a,b):return (a[0]*b[0],a[1]^b[1])
    def field(q):
        t,mask=q
        x=[1]+[1 if ((mask>>(i-1))&1)==0 else -1 for i in range(1,p)]
        return [sum(A[i][j]*t*x[i]*x[j] for j in range(p) if j!=i) for i in range(p)]
    fs={q:field(q) for q in endpoints}
    planes=[]
    for q0 in endpoints:
      rel=[mul(q0,q) for q in endpoints]
      relset=set(rel)
      one=(1,0)
      for a,b in combinations([z for z in rel if z!=one],2):
        if a==b:continue
        ab=mul(a,b)
        if ab in relset:
          P=tuple(sorted((q0,mul(q0,a),mul(q0,b),mul(q0,ab))))
          if all(any(fs[q][i]==0 for q in P) for i in range(p)):
            planes.append(P)
    print('p',p,'M',M,'energies',sorted(set(vals.values())),'endpoints',len(endpoints),'cover planes',len(set(planes)))
    if planes:
      P=planes[0]
      print(P)
      for q in P:print(q,fs[q])

for p in [5,13,17,29]:
    if p<=17: analyze(p)
