import numpy as np,math
z=np.load('/home/math/quadra/tmp/rowcert/cand.npz');L=z['L'];lam=5.382104195764755;m=4;n=4;q0=16
w=np.exp(-lam*(L-L.min()));q=w/w.sum(); idx=np.arange(len(L),dtype=np.uint64); rows=np.stack([((idx>>(i*n))&15).astype(int) for i in range(m)])
chi=np.empty((16,16))
for s in range(16):chi[s]=1-2*(np.bitwise_count((np.arange(16)&s).astype(np.uint8))%2).astype(int)
means=np.array([[np.dot(q,chi[s,rows[i]]) for s in range(1,16)] for i in range(m)])
for i in range(m):
 for j in range(i):
  best=[]
  for s in range(1,16):
   for t in range(1,16):
    c=np.dot(q,chi[s,rows[i]]*chi[t,rows[j]])-means[i,s-1]*means[j,t-1]
    best.append((abs(c),c,s,t))
  print(i,j,sorted(best,reverse=True)[:10])
# brute combos based top pair candidates perhaps print mean
