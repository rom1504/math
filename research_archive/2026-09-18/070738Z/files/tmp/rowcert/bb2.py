import json,math,heapq,time,sys
D=json.load(open('/home/math/quadra/tmp/rowcert/ivq.json')); L=D['logQ_walsh_bounds']
# f coefficients negative log: c0,J pair,k cross,K
c0=(-L[0][1],-L[0][0]);J=(L[3][0],L[3][1]);k=(L[5][0],L[5][1]);K=(-L[15][1],-L[15][0])
def negH(x):
 p=(1+x)/2
 if p<=0 or p>=1:return 0
 return p*math.log(p)+(1-p)*math.log(1-p)
def mul(a,b):return (min(a[0]*b[0],a[0]*b[1],a[1]*b[0],a[1]*b[1]),max(a[0]*b[0],a[0]*b[1],a[1]*b[0],a[1]*b[1]))
def sq(a):
 l,u=a
 return (0 if l<=0<=u else min(l*l,u*u),max(l*l,u*u))
def lb(box):
 u,v=box;u2=sq(u);v2=sq(v);uv=mul(u,v);u2v2=mul(u2,v2)
 z=c0[0]-J[1]*(u2[1]+v2[1])-4*k[1]*uv[1]+K[0]*u2v2[0]
 for a in box:
  x=0 if a[0]<=0<=a[1] else min(a,key=abs)
  z+=2*negH(x)
 return math.nextafter(z,-math.inf)
def run(target):
 Q=[(lb(((0,1),(0,1))),0,((0,1),(0,1)))];n=1;proc=0;t=time.time()
 while Q:
  lo,_,b=heapq.heappop(Q);proc+=1
  if lo>=target:continue
  widths=[x[1]-x[0] for x in b];i=0 if widths[0]>=widths[1] else 1;l,u=b[i];mid=(l+u)/2
  for seg in [(l,mid),(mid,u)]:
   ch=b[:i]+(seg,)+b[i+1:];x=lb(ch)
   if x<target:heapq.heappush(Q,(x,n,ch));n+=1
  if proc%100000==0: print(proc,len(Q),Q[0][0],max(widths),flush=True)
 print('proved',target,'nodes',proc,'time',time.time()-t)
run(float(sys.argv[1]) if len(sys.argv)>1 else 1.0)
