import json,numpy as np,itertools,math
for beta in [2,4]:
 d=json.load(open(f'/home/math/quadra/tmp/rowcert/n9b{beta}.json'));co=np.array(d['coef']);f=-co
 C=np.zeros((4,4))
 for i in range(4):
  for j in range(i):
   rest=[k for k in range(4) if k not in [i,j]];best=0
   for z in itertools.product([-1,1],repeat=2):
    s=0
    for mask in range(16):
     if (mask>>i&1) and (mask>>j&1):s+=f[mask]*math.prod(z[rest.index(k)] for k in rest if mask>>k&1)
    best=max(best,abs(4*s))
   C[i,j]=C[j,i]=best
 print(beta,'C',C,'rho',np.linalg.eigvalsh(C)[-1],'Iuniform',f[0]-4*math.log(2))
