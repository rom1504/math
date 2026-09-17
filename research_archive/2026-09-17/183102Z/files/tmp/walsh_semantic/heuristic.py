import numpy as np
from probe import mats


def local(m, labels, edges, roots=(), starts=500, seed=0):
    rng = np.random.default_rng(seed)
    cs = [mats(m, a)[0] for a in labels]
    w = mats(m, labels[0])[1]
    n = len(w)
    total = n * len(labels)
    K = np.zeros((total, total), dtype=np.int64)
    for v, c in enumerate(cs):
        K[v*n:(v+1)*n, v*n:(v+1)*n] += c
    for u, v in edges:
        K[u*n:(u+1)*n, v*n:(v+1)*n] += w
        K[v*n:(v+1)*n, u*n:(u+1)*n] += w
    omega = (1 << m) - 1
    y = np.array([1 if ((bin(u&v).count('1')+bin(omega&u).count('1'))&1)==0 else -1
                  for u in range(1<<m) for v in range(1<<m)], dtype=np.int64)
    h=np.zeros(total,dtype=np.int64)
    for v,mult in roots: h[v*n:(v+1)*n]+=mult*(w@y)
    best=-10**18
    for _ in range(starts):
        x=rng.choice([-1,1],size=total).astype(np.int64)
        for sweep in range(100):
            order=rng.permutation(total); changed=False
            field=K@x+h
            for i in order:
                # objective 1/2 xKx + hx; flip delta=-2*x_i*(field_i-Kii*x_i)
                delta=-2*x[i]*(field[i]-K[i,i]*x[i])
                if delta>0:
                    old=x[i]; x[i]=-old; field += -2*old*K[:,i]; changed=True
            if not changed: break
        val=int(x@(K@x)//2+h@x)
        best=max(best,val)
    return best


if __name__=='__main__':
  for labels in ([4,2],[4,7]):
    for edges,roots in [([(0,1)],()), ([(0,1)],((0,1),)), ([],((0,1),(1,1))), ([(0,1)],((0,4),))]:
      print(labels,edges,roots,local(3,labels,edges,roots))
