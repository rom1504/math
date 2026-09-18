import json,math,heapq,time,itertools
D=json.load(open('/home/math/quadra/tmp/rowcert/ivq.json')); logb=D['logQ_walsh_bounds'];
# f=-logQ coefficient intervals
coef=[(-b[1],-b[0]) for b in logb]

def mulint(a,b):return (min(a[0]*b[0],a[0]*b[1],a[1]*b[0],a[1]*b[1]),max(a[0]*b[0],a[0]*b[1],a[1]*b[0],a[1]*b[1]))
def negH(m):
 p=(1+m)/2
 if p<=0 or p>=1:return 0.0
 return p*math.log(p)+(1-p)*math.log(1-p)
def lb(box):
 v=coef[0][0]
 for mask in range(1,16):
  prod=(1.,1.)
  for i in range(4):
   if mask>>i&1:prod=mulint(prod,box[i])
  c=coef[mask];v+=min(c[0]*prod[0],c[0]*prod[1],c[1]*prod[0],c[1]*prod[1])
 for l,u in box:
  z=0 if l<=0<=u else (l if abs(l)<abs(u) else u)
  v+=negH(z)
 return math.nextafter(v,-math.inf)
def exact(m):
 v=(coef[0][0]+coef[0][1])/2
 for mask in range(1,16):v+=(coef[mask][0]+coef[mask][1])/2*math.prod(m[i] for i in range(4) if mask>>i&1)
 return v+sum(negH(x) for x in m)
def run(target=1.0,maxnodes=20000000):
 start=((0.,1.),(-1.,1.),(-1.,1.),(-1.,1.));counter=0;Q=[(lb(start),counter,start)];counter+=1;processed=0;t=time.time();minwidth=2
 while Q and processed<maxnodes:
  low,_,box=heapq.heappop(Q);processed+=1
  if low>=target:continue
  # choose split maximizing width and interaction degree roughly
  widths=[u-l for l,u in box];i=max(range(4),key=lambda i:widths[i])
  l,u=box[i];mid=(l+u)/2
  for inter in [(l,mid),(mid,u)]:
   child=box[:i]+(inter,)+box[i+1:];cl=lb(child)
   if cl<target:heapq.heappush(Q,(cl,counter,child));counter+=1
  if processed%100000==0:print(processed,'queue',len(Q),'min',Q[0][0] if Q else None,'width',max(widths),'time',time.time()-t,flush=True)
 print('DONE',not Q,'processed',processed,'queue',len(Q),'min',Q[0][0] if Q else None,'time',time.time()-t)
run(float(__import__('sys').argv[1]) if len(__import__('sys').argv)>1 else 1.0)
