import sys,math
sys.path.insert(0,'/home/math/quadra/extremal_information/experiments')
import actual_child_bridge_law_exact as ex
import actual_child_row_product_shadow as sh
import numpy as np,mpmath as mp
from scipy.optimize import minimize
mp.mp.dps=50

def setup(N,beta,eps=-1):
 m=N//2;n=N-m; mats=[]
 for z in [m,n]:
  sp=ex.build_signing_space(z);c,_=ex.thermal_minimizer_classes(sp,str(beta),N);mats.append(np.asarray(c[0]['representative_matrix'],np.int8))
 A,D=mats;t=beta/math.sqrt(N);x=ex.projective_spins(m).astype(float);y=ex.projective_spins(n).astype(float);ea=ex.energies_for_matrix(A,x).astype(float);ed=ex.energies_for_matrix(D,y).astype(float)
 internal=ea[:,None]+eps*ed[None,:]
 def fg(z):
  M=z.reshape(m,n)
  cross=x@M@y.T
  u=t*(internal+cross)
  umax=np.max(np.abs(u)); # log mean cosh stable
  # direct fine u small
  c=np.cosh(u);Z=np.mean(c);L=math.log(Z)
  # derivative logZ: mean t*sinh(u)*xiyj /Z
  w=np.sinh(u)/np.sum(c)
  gradL=t*np.einsum('ab,ai,bj->ij',w,x,y)
  zz=np.clip(z,-1+1e-15,1-1e-15)
  I=.5*(1+zz)*np.log1p(zz)+.5*(1-zz)*np.log1p(-zz)
  gradI=np.arctanh(zz)
  return L+np.sum(I)/lam,(gradL.reshape(-1)+gradI/lam)
 return m,n,fg
for N,beta,lam0 in [(8,4.,5.382104195764755),(8,2.,1.),(8,4.,1.),(9,2.,1.),(9,4.,1.)]:
 lam=lam0;m,n,fg=setup(N,beta)
 best=None
 for seed in range(4):
  rng=np.random.default_rng(seed);z0=np.zeros(m*n) if seed==0 else rng.uniform(-.5,.5,m*n)
  r=minimize(lambda z:fg(z),z0,jac=True,method='L-BFGS-B',bounds=[(-.999999999,.999999999)]*(m*n),options={'ftol':1e-14,'gtol':1e-11,'maxiter':10000,'maxls':100})
  if best is None or r.fun<best.fun:best=r
 print('case',N,beta,lam,'val',best.fun,'grad',np.max(abs(best.jac)),'nit',best.nit,'success',best.success,'M range',best.x.min(),best.x.max(),'distinct',np.unique(np.round(abs(best.x),6)))
