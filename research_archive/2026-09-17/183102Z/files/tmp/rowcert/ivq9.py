import mpmath as mp,numpy as np,math,json
mp.iv.dps=60;z=np.load('/home/math/quadra/tmp/rowcert/signatures9.npz');sigs=z['sigs'];counts=z['counts'];width=21

def cosh(x):return (mp.iv.exp(x)+mp.iv.exp(-x))/2
def bounds(x):
 lo=mp.libmp.to_float(x._mpi_[0],rnd='f');hi=mp.libmp.to_float(x._mpi_[1],rnd='c');return [math.nextafter(lo,-math.inf),math.nextafter(hi,math.inf)]
out={}
for bi,beta in enumerate(['2','4']):
 t=mp.iv.mpf(beta)/mp.iv.sqrt(mp.iv.mpf('9'));cv=[cosh(t*c) for c in range(21)];ivv=[cosh(t*a) for a in range(sigs.shape[1]//width)];weights=[]
 for sig in sigs:
  Z=mp.iv.mpf('0')
  for ix,cnt in enumerate(sig):
   if cnt:Z+=int(cnt)*ivv[ix//width]*cv[ix%width]
  Z/=128;weights.append(Z**(-1))
 nums=[sum((int(counts[k,bi,code])*weights[k] for k in range(len(weights))),mp.iv.mpf('0')) for code in range(16)];den=sum(nums,mp.iv.mpf('0'));Q=[x/den for x in nums];logs=[mp.iv.log(x) for x in Q];co=[]
 for mask in range(16):co.append(sum((((-1)**bin(s&mask).count('1'))*logs[s] for s in range(16)),mp.iv.mpf('0'))/16)
 out[beta]={'Q_bounds':[bounds(x) for x in Q],'logQ_walsh_bounds':[bounds(x) for x in co]}
 print(beta,'Q0',Q[0],'coef0',co[0],'nonzero',[(i,co[i]) for i in range(16) if i in [0,3,5,6,9,10,12,15]])
open('/home/math/quadra/tmp/rowcert/ivq9.json','w').write(json.dumps(out,indent=2)+'\n')
