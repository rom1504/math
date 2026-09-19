import numpy as np,itertools,math
from scipy.optimize import minimize_scalar,differential_evolution,shgo
z=np.load('/home/math/quadra/tmp/rowcert/cand.npz');p=z['p'];lam=5.382104195764755
m=4;q=16
Ks={
(0,1):12.3700902097945,(0,2):8.892458332610357,(0,3):9.901839036283434,(1,2):9.169817644674058,(1,3):9.566248371388621,(2,3):11.808414730006925,
(0,1,2):12.367470266132912,(0,1,3):13.735790857826085,(0,2,3):12.105606008335538,(1,2,3):12.366644821952725,(0,1,2,3):18.457243323930605}
ss=[]
for x in p:
 v=np.array([0.0])
 for a in x:v=np.concatenate([v,v+a])
 ss.append(np.sort(v[(v>1e-15)&(v<1-1e-15)]))
def bin_kl(d,a):
 b=a+d
 if a<=0 or b>=1:return math.inf
 return b*math.log(b/a)+(1-b)*math.log((1-b)/(1-a))
def psi(i,d):
 if d<1e-14:return 0.
 arr=ss[i]; hi=np.searchsorted(arr,1-d,side='right')
 if not hi:return math.inf
 # continuous optimum
 r=minimize_scalar(lambda a:bin_kl(d,a),bounds=(1e-14,1-d-1e-14),method='bounded',options={'xatol':1e-13})
 ix=np.searchsorted(arr[:hi],r.x)
 cand=arr[max(0,ix-4):min(hi,ix+5)]
 return min(bin_kl(d,a) for a in cand)
def fun(d):
 return sum(psi(i,float(d[i])) for i in range(m))/lam-sum(K*np.prod([d[i] for i in S]) for S,K in Ks.items())
for i in range(4):
 print('row',i,'sub around half',ss[i][np.argmin(abs(ss[i]-.5))],[(d,psi(i,d)) for d in [.001,.01,.1,.25,.5,.9]])
res=differential_evolution(fun,[(0,1-min(x)) for x in p],tol=1e-10,popsize=30,maxiter=1000,polish=True,seed=3,workers=1);print(res.fun,res.x)
print('corners')
for b in itertools.product([0,1],repeat=4):
 d=np.array([b[i]*(1-p[i].min()) for i in range(4)]);print(b,fun(d))
