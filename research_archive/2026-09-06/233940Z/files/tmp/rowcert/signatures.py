import sys,numpy as np,math,time,collections
sys.path.insert(0,'/home/math/quadra/extremal_information/experiments')
import actual_child_bridge_law_exact as ex
import mpmath as mp
mp.mp.dps=80
N=8;m=n=4;beta=4.;eps=-1
sp=ex.build_signing_space(4);cs,_=ex.thermal_minimizer_classes(sp,str(beta),N);A=np.asarray(cs[0]['representative_matrix'],np.int8)
x=ex.projective_spins(m).astype(np.int16);y=ex.projective_spins(n).astype(np.int16);ea=ex.energies_for_matrix(A,x);ed=ex.energies_for_matrix(A,y);internal=np.abs(ea[:,None]+eps*ed[None,:]).reshape(-1)
maxi=int(internal.max());maxc=m*n; width=maxc+1;nb=(maxi+1)*width
ranks=np.array([(xx[:,None]*yy[None,:]).reshape(-1) for xx in x for yy in y],dtype=np.int16) #64,16
features=[2,1,4,8];counts={}; t0=time.time()
for lo in range(0,1<<(m*n),1024):
 masks=np.arange(lo,min(lo+1024,1<<(m*n)),dtype=np.uint64)
 bits=((masks[:,None]>>np.arange(m*n,dtype=np.uint64))&1).astype(np.int16);B=1-2*bits
 cross=np.abs(B@ranks.T).astype(np.int16) #batch,64
 bins=internal[None,:]*width+cross
 for k,mask in enumerate(masks):
  sig=np.bincount(bins[k],minlength=nb).astype(np.uint8).tobytes()
  row=[(int(mask)>>(i*n))&15 for i in range(m)]
  code=sum(((bin(row[i]&features[i]).count('1')&1)<<i) for i in range(m))
  if sig not in counts:counts[sig]=np.zeros(16,dtype=np.int64)
  counts[sig][code]+=1
print('signatures',len(counts),'total',sum(v.sum() for v in counts.values()),'time',time.time()-t0,'maxi',maxi,'nb',nb)
np.savez('/home/math/quadra/tmp/rowcert/signatures.npz',sigs=np.frombuffer(b''.join(counts),dtype=np.uint8).reshape(len(counts),nb),counts=np.stack(list(counts.values())),internal=internal)
