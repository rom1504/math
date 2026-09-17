import sys,math
sys.path.insert(0,'/home/math/quadra/extremal_information/experiments')
import actual_child_bridge_law_exact as ex
import actual_child_row_product_shadow as sh
import numpy as np,mpmath as mp
mp.mp.dps=50

def pressure(N,beta,eps=-1):
 m=N//2;n=N-m
 spaces={z:ex.build_signing_space(z) for z in {m,n}}
 mats=[]
 for z in [m,n]:
  c,_=ex.thermal_minimizer_classes(spaces[z],str(beta),N);mats.append(np.asarray(c[0]['representative_matrix'],dtype=np.int8))
 L,_=ex.bridge_pressures(mats[0],mats[1],beta,N,eps)
 return m,n,L,sh.pressure_tensor(L,m,1<<n)

def rect(M):
 # axes first two q,q, rest flattened
 q=M.shape[0]; rest=int(np.prod(M.shape[2:]));X=M.reshape(q,q,rest)
 best=0
 for a in range(q):
  for ap in range(a):
   D=X[a]-X[ap]
   v=np.max(D,axis=0)-np.min(D,axis=0)
   best=max(best,float(np.max(v)))
 return best

def cmat(T):
 m=T.ndim; C=np.zeros((m,m))
 for i in range(m):
  for j in range(i):
   axes=[i,j]+[k for k in range(m) if k not in [i,j]]
   C[i,j]=C[j,i]=rect(np.transpose(T,axes))
 return C
for N in [6,7,8,9]:
 for beta in [2.,4.]:
  m,n,L,T=pressure(N,beta)
  C=cmat(T);print('N beta',N,beta,'m q',m,1<<n,'C',C,'rho',max(np.linalg.eigvalsh(C)),'lambda1 criterion',max(np.linalg.eigvalsh(C))/4)
