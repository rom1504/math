import itertools

A9 = [
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

def energy(A,x):
    return sum(A[i][j]*x[i]*x[j] for i in range(len(A)) for j in range(i+1,len(A)))

def grounds(A, omit=None):
    n=len(A); vs=[i for i in range(n) if i != omit]
    vals=[]; m=-1
    for bits in itertools.product((-1,1), repeat=len(vs)-1):
        x={vs[0]:1, **{v:b for v,b in zip(vs[1:],bits)}}
        h=sum(A[u][v]*x[u]*x[v] for ii,u in enumerate(vs) for v in vs[ii+1:])
        ah=abs(h)
        if ah>m: m=ah; vals=[]
        if ah==m:
            t=1 if h>=0 else -1
            vals.append((t,x))
    return m, vals

def extensions(A,i):
    m,gs=grounds(A,i); out=[]
    for t,x0 in gs:
        h=t*sum(A[i][j]*x0[j] for j in x0)
        xi=1 if h>=0 else -1
        x=dict(x0); x[i]=xi
        q=[[0]*len(A) for _ in A]
        for u in range(len(A)):
            for v in range(u+1,len(A)):
                q[u][v]=q[v][u]=t*x[u]*x[v]
        score=sum(A[u][v]*q[u][v] for u in range(len(A)) for v in range(u+1,len(A)))
        fields=[sum(A[j][k]*q[j][k] for k in range(len(A)) if k!=j) for j in range(len(A))]
        out.append((score,fields,q))
    return m,out

def audit(A, name, expected_M, expected_child_M):
    n=len(A); M,_=grounds(A); ex=[]
    for i in range(n):
        mi, xs=extensions(A,i); ex.append((mi,xs))
    assert M == expected_M
    assert all(mi == expected_child_M for mi,_ in ex)
    ds=[M-mi for mi,_ in ex]
    # For every possible optimal extension q_i, verify u_ii=0,
    # u_ij>=0, and the row conservation identity.
    all_u=set(); row_sums=set(); slacks=set()
    for i,(mi,xs) in enumerate(ex):
      for score,fields,q in xs:
        s=M-score
        assert 0 <= s <= ds[i]
        assert fields[i] == ds[i]-s
        us=[s+fields[j]-ds[j] for j in range(n)]
        assert us[i] == 0 and min(us) >= 0
        assert sum(us) == 2*M-sum(ds)+(n-2)*s
        all_u.update(us); row_sums.add(sum(us)); slacks.add(s)
    mx=(-1,None)
    for i in range(n):
      for j in range(i+1,n):
       for si,fi,qi in ex[i][1]:
        for sj,fj,qj in ex[j][1]:
          sli=M-si; slj=M-sj
          uij=sli+fi[j]-ds[j]
          uji=slj+fj[i]-ds[i]
          val=uij+uji
          if val>mx[0]:mx=(val,(i,j,ds[i],ds[j],sli,slj,uij,uji))
    print(name, 'n=',n, 'M=',M, 'child_M=',[z for z,_ in ex],
          'ground_extensions=',[len(z) for _,z in ex], 'd=',ds,
          's=',sorted(slacks), 'u_values=',sorted(all_u),
          'row_sums=',sorted(row_sums), 'max_pair=',mx)

A8 = [
[0,1,1,1,1,1,1,1],
[1,0,1,-1,1,1,-1,-1],
[1,1,0,1,-1,1,-1,-1],
[1,-1,1,0,-1,-1,-1,1],
[1,1,-1,-1,0,-1,1,-1],
[1,1,1,-1,-1,0,1,1],
[1,-1,-1,-1,1,1,0,1],
[1,-1,-1,1,-1,1,1,0],
]

audit(A8, 'A8', 10, 9)
audit(A9, 'A9', 12, 12)
