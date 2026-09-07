import sys,numpy as np,math,time,collections,mpmath as mp
sys.path.insert(0,'/home/math/quadra/extremal_information/experiments');import actual_child_bridge_law_exact as ex
mp.mp.dps=80;N=9;m=4;n=5;eps=-1
As=[]
for z in [m,n]:
 c,_=ex.thermal_minimizer_classes(ex.build_signing_space(z),'2',N);As.append(np.array(c[0]['representative_matrix'],np.int8))
x=ex.projective_spins(m).astype(np.int16);y=ex.projective_spins(n).astype(np.int16);ea=ex.energies_for_matrix(As[0],x);ed=ex.energies_for_matrix(As[1],y);internal=np.abs(ea[:,None]+eps*ed[None,:]).reshape(-1);maxi=int(internal.max());maxc=m*n;width=maxc+1;nb=(maxi+1)*width
ranks=np.array([(xx[:,None]*yy[None,:]).reshape(-1) for xx in x for yy in y],dtype=np.int16)
feature_sets=[[1,1,2,2],[1,1,4,4]];counts={};t0=time.time();tot=1<<(m*n)
for lo in range(0,tot,2048):
 masks=np.arange(lo,min(lo+2048,tot),dtype=np.uint64);bits=((masks[:,None]>>np.arange(m*n,dtype=np.uint64))&1).astype(np.int16);B=1-2*bits;cross=np.abs(B@ranks.T).astype(np.int16);bins=internal[None,:]*width+cross
 for kk,mask in enumerate(masks):
  sig=np.bincount(bins[kk],minlength=nb).astype(np.uint8).tobytes();row=[(int(mask)>>(i*n))&31 for i in range(m)];codes=[]
  for fs in feature_sets:codes.append(sum(((bin(row[i]&fs[i]).count('1')&1)<<i) for i in range(m)))
  if sig not in counts:counts[sig]=np.zeros((2,16),dtype=np.int64)
  for j,c in enumerate(codes):counts[sig][j,c]+=1
 if lo%(2048*64)==0:print(lo,len(counts),time.time()-t0,flush=True)
print('done',len(counts),time.time()-t0)
np.savez('/home/math/quadra/tmp/rowcert/signatures9.npz',sigs=np.frombuffer(b''.join(counts),dtype=np.uint8).reshape(len(counts),nb),counts=np.stack(list(counts.values())),internal=internal,feature_sets=np.array(feature_sets))
