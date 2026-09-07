import sys,math
sys.path.insert(0,'/home/math/quadra/extremal_information/experiments')
import actual_child_bridge_law_exact as ex
import actual_child_row_product_shadow as sh
import numpy as np, mpmath as mp
mp.mp.dps=50
N=8;m=n=4;beta=4.;lam=5.382104195764755
sp=ex.build_signing_space(4);cs,_=ex.thermal_minimizer_classes(sp,str(beta),N);A=np.asarray(cs[0]['representative_matrix'],np.int8);L,_=ex.bridge_pressures(A,A,beta,N,-1);T=sh.pressure_tensor(L,m,1<<n)
p=[np.full(16,1/16) for _ in range(m)]
for sweep in range(1000):
 old=[x.copy() for x in p]
 for i in range(m):p[i]=sh.softmax(-lam*sh.contract_except(T,p,i))
 if max(np.sum(abs(p[i]-old[i])) for i in range(m))<1e-14:break
print('sweep',sweep,'obj',sh.product_objective(T,p,lam))
np.set_printoptions(precision=12,suppress=True)
for i,x in enumerate(p):print(i,np.argsort(x)[::-1],x[np.argsort(x)[::-1]])
np.savez('/home/math/quadra/tmp/rowcert/cand.npz',L=L,p=np.stack(p),T=T)
