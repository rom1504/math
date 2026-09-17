import sys,math,time
from pathlib import Path
sys.path.insert(0,'/home/math/quadra/extremal_information/experiments')
import actual_child_bridge_law_exact as ex
import numpy as np, mpmath as mp
import cvxpy as cp
mp.mp.dps=50
N=8;m=n=4;beta=4.;lam=5.382104195764755
sp=ex.build_signing_space(4)
cls,cert=ex.thermal_minimizer_classes(sp,str(beta),N)
A=np.asarray(cls[0]['representative_matrix'],dtype=np.int8)
L,aud=ex.bridge_pressures(A,A,beta,N,-1)
q=1<<n; d=m*n
idx=np.arange(1<<d,dtype=np.uint64)
rowidx=np.stack([((idx>>(i*n))&(q-1)).astype(np.int64) for i in range(m)])
phi=[cp.Variable(q) for i in range(m)]
c=cp.Variable()
expr=c
for i in range(m): expr += phi[i][rowidx[i]]
cons=[expr <= L]
# gauges
cons += [phi[i][0]==0 for i in range(m)]
obj=c
for i in range(m): obj += -(cp.log_sum_exp(-lam*phi[i])-math.log(q))/lam
prob=cp.Problem(cp.Maximize(obj),cons)
print('DCP',prob.is_dcp(),'shape',len(L))
t=time.time();val=prob.solve(solver='CLARABEL',verbose=True,tol_gap_abs=1e-9,tol_feas=1e-9,tol_gap_rel=1e-9,max_iter=500);print('status',prob.status,'val',val,'time',time.time()-t)
phis=np.stack([x.value for x in phi]); residual=L-cp.Constant(0).value if False else L-np.sum([phis[i,rowidx[i]] for i in range(m)],axis=0)
cvalid=float(np.min(residual))
lb=cvalid+sum(-(np.log(np.mean(np.exp(-lam*phis[i]))))/lam for i in range(m))
print('c solve',c.value,'cvalid',cvalid,'lb',lb,'minL',L.min(),'target V',ex.negative_moment_soft_pressure(L,lam),'Lmean',L.mean())
print('violation',c.value-cvalid,'active',np.sum(residual<cvalid+1e-7))
np.savez('/home/math/quadra/tmp/rowcert/dual.npz',L=L,phis=phis,cvalid=cvalid,rowidx=rowidx)
