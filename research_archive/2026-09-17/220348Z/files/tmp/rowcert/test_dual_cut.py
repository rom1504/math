import sys,math,time
sys.path.insert(0,'/home/math/quadra/extremal_information/experiments')
import actual_child_bridge_law_exact as ex
import numpy as np, mpmath as mp
import cvxpy as cp
mp.mp.dps=50
N=8;m=n=4;beta=4.;lam=5.382104195764755
sp=ex.build_signing_space(4); cls,_=ex.thermal_minimizer_classes(sp,str(beta),N)
A=np.asarray(cls[0]['representative_matrix'],dtype=np.int8); L,_=ex.bridge_pressures(A,A,beta,N,-1)
q=1<<n; idx=np.arange(len(L),dtype=np.uint64); rowidx=np.stack([((idx>>(i*n))&(q-1)).astype(np.int64) for i in range(m)])
phi=[cp.Variable(q) for i in range(m)]; c=cp.Variable();
obj=c+sum(-(cp.log_sum_exp(-lam*x)-math.log(q))/lam for x in phi)
active=set(np.argsort(L)[:256].tolist())
phis=np.zeros((m,q));cv=0
for it in range(100):
 ai=np.array(sorted(active),dtype=np.int64)
 cons=[c+sum(phi[i][rowidx[i,ai]] for i in range(m)) <= L[ai]]+[phi[i][0]==0 for i in range(m)]
 prob=cp.Problem(cp.Maximize(obj),cons)
 val=prob.solve(solver='CLARABEL',tol_gap_abs=1e-10,tol_feas=1e-10,tol_gap_rel=1e-10,max_iter=500,verbose=False)
 phis=np.stack([x.value for x in phi]); cv=float(c.value)
 resid=L-np.sum([phis[i,rowidx[i]] for i in range(m)],axis=0)-cv
 order=np.argsort(resid)
 neg=order[resid[order]<-1e-8]
 print(it,prob.status,val,'active',len(active),'minres',resid[order[0]],'neg',len(neg))
 if len(neg)==0:break
 # add worst 256 plus diverse
 active.update(order[:min(1024,len(order))].tolist())
res0=L-np.sum([phis[i,rowidx[i]] for i in range(m)],axis=0)
cvalid=float(np.min(res0)); lb=cvalid+sum(-np.log(np.mean(np.exp(-lam*phis[i])))/lam for i in range(m))
print('c',cv,'valid',cvalid,'LB',lb,'V',ex.negative_moment_soft_pressure(L,lam),'minL',L.min(),'mean',L.mean(),'active',len(active))
np.savez('/home/math/quadra/tmp/rowcert/dual.npz',L=L,phis=phis,cvalid=cvalid,rowidx=rowidx)
