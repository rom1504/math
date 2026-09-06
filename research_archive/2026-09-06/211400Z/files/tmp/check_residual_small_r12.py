from itertools import product


A5 = [
    [0,-1,1,-1,1],
    [-1,0,-1,1,1],
    [1,-1,0,1,1],
    [-1,1,1,0,1],
    [1,1,1,1,0],
]

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


def energy(A, x):
    return sum(A[i][j]*x[i]*x[j] for i in range(len(A)) for j in range(len(A)))


def stats(A):
    xs = [(1,)+tail for tail in product((-1,1), repeat=len(A)-1)]
    es = [(energy(A,x),x) for x in xs]
    P=max(e for e,x in es); N=-min(e for e,x in es); Q=max(P,N)
    pos=[x for e,x in es if e==P]; neg=[x for e,x in es if e==-N]
    return P,N,Q,pos,neg


def principal(A, inds):
    return [[A[i][j] for j in inds] for i in inds]


def root_data(A):
    P,N,Q,ps,ns=stats(A)
    print('parent',len(A),'P,N,Q',P,N,Q,'endpoints',len(ps),len(ns))
    out=[]
    for p in ps:
      for n in ns:
        shores=[]
        for val in (-1,1):
          shores.append([i for i in range(len(A)) if p[i]*n[i]==val])
        if not all(shores):
          continue
        row=[]
        for X in shores:
          Y=[i for i in range(len(A)) if i not in X]
          AX=principal(A,X)
          PX,NX,QX,_,_=stats(AX)
          h=sum(A[i][j]*p[i]*p[j] for i in X for j in X)
          # Fix p on X; expose cross by signs on Y, equivalently l1 of fields into Y.
          fields=[sum(A[i][j]*p[i] for i in X) for j in Y]
          L=sum(abs(v) for v in fields)
          caps=tuple(max(0,2*L-(QX-sig*h)) for sig in (1,-1))
          d=Q-QX
          residual=tuple(max(0,c-d) for c in caps)
          row.append((tuple(X),QX,h,L,caps,d,residual))
        out.append((p,n,row))
    for k,(p,n,row) in enumerate(out):
      if len(A)==5 and any(len(item[0])==1 for item in row):
        fields=tuple(sum(A[i][j]*p[j] for j in range(len(A))) * p[i] for i in range(len(A)))
        print('pair',k,'p',p,'n',n,'gauged fields',fields,row)
      elif len(A)!=5:
        print('pair',k,row)
    return P,N,Q,out


p5,n5,q5,out5=root_data(A5)
assert (p5,n5,q5)==(8,8,8)
target_p=(1,1,1,1,1)
target_n=(1,1,1,1,-1)
target=[row for p,n,row in out5 if p==target_p and n==target_n]
assert len(target)==1
assert target[0]==[
    ((4,),0,0,4,(8,8),8,(0,0)),
    ((0,1,2,3),8,0,4,(0,0),0,(0,0)),
]
residual_totals=[sum(sum(item[-1]) for item in row) for _,_,row in out5]
assert residual_totals.count(0)==5
assert residual_totals.count(4)==20

core=principal(A5,[0,1,2,3])
core_energy=energy(core,(1,1,1,1))
assert stats(core)[2]==8 and core_energy==0
assert (q5-stats(core)[2],stats(core)[2]-core_energy)==(0,8)
assert 8-64*(5**0.5)/25>0

p8,n8,q8,out8=root_data(A8)
assert (p8,n8,q8)==(20,20,20)
assert all(sum(sum(item[-1]) for item in row)==0 for _,_,row in out8)
print('all residual certificates: PASS')
