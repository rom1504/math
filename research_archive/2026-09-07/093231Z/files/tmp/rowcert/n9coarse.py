import sys,numpy as np,math,time,mpmath as mp,json
sys.path.insert(0,'/home/math/quadra/extremal_information/experiments');import actual_child_bridge_law_exact as ex
mp.mp.dps=80

def fwht(a,axis):
 a=np.swapaxes(np.array(a,float,copy=True),axis,-1);n=a.shape[-1];h=1
 while h<n:
  v=a.reshape(*a.shape[:-1],-1,2*h);x=v[...,:h].copy();y=v[...,h:].copy();v[...,:h]=x+y;v[...,h:]=x-y;h*=2
 return np.swapaxes(a,axis,-1)
for beta in [2.,4.]:
 N=9;m=4;n=5;As=[]
 for z in [m,n]:
  sp=ex.build_signing_space(z);c,_=ex.thermal_minimizer_classes(sp,format(beta,'.12g'),N);print('classes',beta,z,len(c));As.append(np.array(c[0]['representative_matrix'],np.int8))
 L,_=ex.bridge_pressures(As[0],As[1],beta,N,-1);w=np.exp(-(L-L.min()));q=w/w.sum();qT=np.transpose(q.reshape((32,)*4),tuple(reversed(range(4)))) # axes row order
 means=[]
 for i in range(4):
  axes=tuple(j for j in range(4) if j!=i);mar=qT.sum(axis=axes);means.append(fwht(mar,0))
 best={}
 for i in range(4):
  for j in range(i):
   axes=tuple(k for k in range(4) if k not in [j,i]);mar=qT.sum(axis=axes)
   # remaining axes in sorted [j,i]
   wal=fwht(fwht(mar,0),1);cov=wal-means[j][:,None]*means[i][None,:];cov[0,:]=0;cov[:,0]=0
   ix=np.unravel_index(np.argmax(abs(cov)),cov.shape);best[(j,i)]=(abs(cov[ix]),cov[ix],ix)
 print('beta',beta,'pairbest',best)
 # choose max matching among three, features from its two pairs
 matchings=[((0,1),(2,3)),((0,2),(1,3)),((0,3),(1,2))];ma=max(matchings,key=lambda mm:sum(best[tuple(sorted(e))][0] for e in mm));feat=[None]*4
 for e in ma:
  j,i=sorted(e);s,t=best[(j,i)][2];feat[j]=s;feat[i]=t
 print('matching',ma,'features',feat)
 idx=np.arange(len(L),dtype=np.uint64);rows=np.stack([((idx>>(i*n))&31).astype(np.uint8) for i in range(m)]);ys=np.stack([(np.bitwise_count(rows[i]&feat[i])%2).astype(np.int8) for i in range(m)]);code=sum(ys[i]<<i for i in range(m));Q=np.bincount(code,weights=q,minlength=16)
 logq=np.log(Q);co=[]
 for mask in range(16):co.append(sum(((-1)**(bin(s&mask).count('1')))*logq[s] for s in range(16))/16)
 print('Q',Q);print('co',[(i,x) for i,x in enumerate(co) if abs(x)>1e-10])
 open(f'/home/math/quadra/tmp/rowcert/n9b{int(beta)}.json','w').write(json.dumps({'beta':beta,'features':[int(x) for x in feat],'Q':Q.tolist(),'coef':co},indent=2)+'\n')
