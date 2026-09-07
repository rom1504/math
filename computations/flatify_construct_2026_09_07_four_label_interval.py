"""Rigorous mpmath interval enclosure, capped diagnostic run if not completed."""
import itertools
import json
import time
from pathlib import Path
from mpmath import iv

iv.dps=25
signs=list(itertools.product([-1,1],repeat=3))
def det(u):
    return 256-120*(u[0]+u[1])*(u[2]+u[3])+225*u[0]*u[1]*u[2]*u[3]
def upper(box):
    a=[iv.mpf([x,y]) for x,y in box]; u=[x*x for x in a]
    # The determinant decreases in each squared amplitude on [0,1]^4.
    dl=det([iv.mpf(y)**2 for x,y in box]); dh=det([iv.mpf(x)**2 for x,y in box])
    D=iv.mpf([dl.a,dh.b])
    B={}
    B[0,0]=iv.mpf(15)/2*(8*u[2]+8*u[3]-15*u[1]*u[2]*u[3])
    B[1,1]=iv.mpf(15)/2*(8*u[2]+8*u[3]-15*u[0]*u[2]*u[3])
    B[2,2]=iv.mpf(15)/2*(8*u[0]+8*u[1]-15*u[0]*u[1]*u[3])
    B[3,3]=iv.mpf(15)/2*(8*u[0]+8*u[1]-15*u[0]*u[1]*u[2])
    B[0,1]=60*(u[2]-u[3]); B[2,3]=60*(u[0]-u[1])
    z=iv.sqrt(30)
    B[0,2]=z*(15*u[1]*u[3]-16); B[0,3]=z*(15*u[1]*u[2]-16)
    B[1,2]=z*(15*u[0]*u[3]-16); B[1,3]=-z*(15*u[0]*u[2]-16)
    diag=sum(B[i,i]*u[i] for i in range(4))
    total=iv.mpf(0)
    for tail in signs:
        s=(1,)+tail
        E=diag+sum(2*B[i,j]*a[i]*a[j]*s[i]*s[j] for i in range(4) for j in range(i+1,4))
        total+=iv.exp(-E/D)
    return total/(8*iv.sqrt(D))

stack=[tuple((0.,1.) for _ in range(4))]
done=0; examined=0; largest=0.; start=time.monotonic(); cap=180
while stack and time.monotonic()-start<cap:
    box=stack.pop()
    if box[0][1]<box[1][0] or box[2][1]<box[3][0] or box[0][1]<box[2][0]:
        continue
    val=upper(box); examined+=1
    if bool(val.b<iv.mpf(1)/12):
        done+=1; largest=max(largest,float(val.b)); continue
    dim=max(range(4),key=lambda i:box[i][1]-box[i][0])
    lo,hi=box[dim]; mid=(lo+hi)/2
    for pair in [(lo,mid),(mid,hi)]:
        child=list(box); child[dim]=pair; stack.append(tuple(child))
    if examined%1000==0:
        print(json.dumps(dict(examined=examined,accepted=done,pending=len(stack),seconds=time.monotonic()-start)),flush=True)
result=dict(status='PASS all boxes' if not stack else 'INCOMPLETE bounded run',
            examined=examined,accepted=done,pending=len(stack),largest_accepted_upper=largest,
            seconds=time.monotonic()-start)
print(json.dumps(result),flush=True)
Path('computations/results/flatify_construct_2026_09_07_four_label_interval.json').write_text(json.dumps(result,indent=2)+'\n')
