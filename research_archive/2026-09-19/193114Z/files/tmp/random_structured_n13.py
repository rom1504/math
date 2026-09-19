import itertools, random, time
import numpy as np

m=6;n=13
edges=[(i,j) for i in range(n) for j in range(i+1,n)]
X=list(range(m)); z=m; Y=list(range(m+1,2*m+1))
spins=np.ones((1<<(n-1),n),dtype=np.int8)
for mask in range(1<<(n-1)):
 for i in range(1,n):
  if mask>>(i-1)&1:spins[mask,i]=-1
edgechars=np.array([spins[:,i]*spins[:,j] for i,j in edges],dtype=np.int8).T
allpairs=list(itertools.combinations(range(m),2))
balanced=[]
for neg in itertools.combinations(range(m),m//2):
 s=np.ones(m,dtype=np.int8);s[list(neg)]=-1;balanced.append(s)

ECOUNT=int(__import__('os').environ.get('ECOUNT','4'))
def rand_graph():
 while True:
  es=random.sample(allpairs,ECOUNT)
  deg=[0]*m
  for i,j in es:deg[i]+=1;deg[j]+=1
  if max(deg)<=2:return es,deg

def rand_bip(degr,degc):
 rs=[i for i,d in enumerate(degr) for _ in range(d)]
 cs0=[j for j,d in enumerate(degc) for _ in range(d)]
 for _ in range(30):
  cs=cs0[:];random.shuffle(cs)
  ps=list(zip(rs,cs))
  if len(set(ps))==len(ps):return ps

def make():
 while True:
  eg,dg=rand_graph(); eh,dh=rand_graph(); hb=rand_bip(dg,dh)
  if hb is not None:break
 A=np.ones((n,n),dtype=np.int8);np.fill_diagonal(A,0)
 # B_X positive except G; B_Y is negative except positives E_h = -B(graph)
 for i,j in eg:A[i,j]=A[j,i]=-1
 for i,j in itertools.combinations(range(m),2):A[Y[i],Y[j]]=A[Y[j],Y[i]]=-1
 for i,j in eh:A[Y[i],Y[j]]=A[Y[j],Y[i]]=1
 for i,j in hb:A[X[i],Y[j]]=A[Y[j],X[i]]=-1
 sx=random.choice(balanced);sy=random.choice(balanced)
 for i in range(m):A[X[i],z]=A[z,X[i]]=sx[i];A[Y[i],z]=A[z,Y[i]]=sy[i]
 return A,(eg,eh,hb,sx,sy)

def avec(A):return np.array([A[i,j] for i,j in edges],dtype=np.int16)
TARGET=m*m-4*ECOUNT
best=999
start=time.time()
for it in range(200000):
 A,data=make(); en=edgechars@avec(A); mm=int(np.abs(en).max())
 if mm<best:
  best=mm;print('best',best,'it',it,'elapsed',time.time()-start,'range',en.min(),en.max(),flush=True)
 if mm<=TARGET:
  print('FOUND',it)
  print(A.tolist())
  print(data[:3],data[3].tolist(),data[4].tolist())
  vals,cts=np.unique(en,return_counts=True);print(list(zip(vals.tolist(),cts.tolist())))
  break
else:print('none best',best)
