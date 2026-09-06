import mpmath as mp,numpy as np,math,json
mp.iv.dps=80
z=np.load('/home/math/quadra/tmp/rowcert/signatures.npz');sigs=z['sigs'];counts=z['counts'];
t=mp.iv.mpf('4')/mp.iv.sqrt(mp.iv.mpf('8')); lam=mp.iv.mpf('5.382104195764755')
def cosh(x):return (mp.iv.exp(x)+mp.iv.exp(-x))/2
cv=[cosh(t*c) for c in range(17)];ivv=[cosh(t*a) for a in range(9)]
weights=[]
for sig in sigs:
 Z=mp.iv.mpf('0')
 for ix,cnt in enumerate(sig):
  if cnt: Z += int(cnt)*ivv[ix//17]*cv[ix%17]
 Z/=64
 weights.append(Z**(-lam))
nums=[]
for code in range(16):nums.append(sum((int(counts[k,code])*weights[k] for k in range(len(weights))),mp.iv.mpf('0')))
den=sum(nums,mp.iv.mpf('0'));Q=[x/den for x in nums];logs=[mp.iv.log(x) for x in Q]
co=[]
for mask in range(16):
 v=sum((((-1)**(bin(s&mask).count('1')))*logs[s] for s in range(16)),mp.iv.mpf('0'))/16
 co.append(v)
print('Q');[print(i,x) for i,x in enumerate(Q)];print('co');[print(i,x) for i,x in enumerate(co) if not (float(mp.libmp.to_float(x._mpi_[0]))<=0<=float(mp.libmp.to_float(x._mpi_[1])) and abs(float(mp.libmp.to_float(x._mpi_[1]))-float(mp.libmp.to_float(x._mpi_[0])))<1e-50)]
# outward float bounds
def bounds(x):
 lo=mp.libmp.to_float(x._mpi_[0],rnd='f');hi=mp.libmp.to_float(x._mpi_[1],rnd='c');return [math.nextafter(lo,-math.inf),math.nextafter(hi,math.inf)]
out={'Q_bounds':[bounds(x) for x in Q],'logQ_walsh_bounds':[bounds(x) for x in co]}
open('/home/math/quadra/tmp/rowcert/ivq.json','w').write(json.dumps(out,indent=2)+'\n')
