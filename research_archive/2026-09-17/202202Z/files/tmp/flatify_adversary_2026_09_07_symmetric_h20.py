"""Exact F9 Paley conference and symmetric H20 for mesoscopic-gap audit."""
import json
from pathlib import Path
import numpy as np

def mul9(x,y):
    a,b=x%3,x//3;c,d=y%3,y//3
    return (a*c-b*d)%3+3*((a*d+b*c)%3)

def sub9(x,y):return (x%3-y%3)%3+3*((x//3-y//3)%3)
squares={mul9(x,x) for x in range(1,9)}
assert len(squares)==4 and 2 in squares
chi=lambda x:0 if x==0 else (1 if x in squares else -1)
C=np.ones((10,10),dtype=np.int64);C[0,0]=0
for x in range(9):
    for y in range(9):C[x+1,y+1]=chi(sub9(x,y))
assert np.array_equal(C,C.T)
assert np.array_equal(C@C,9*np.eye(10,dtype=np.int64))
I=np.eye(10,dtype=np.int64)
H=np.block([[C+I,C-I],[C-I,-C-I]])
assert np.array_equal(H,H.T) and np.all(np.abs(H)==1)
assert np.array_equal(H@H,20*np.eye(20,dtype=np.int64))
K=np.array([[1]],dtype=np.int64)
for _ in range(6):K=np.block([[K,K],[K,-K]])
assert np.array_equal(K,K.T) and np.array_equal(K@K,64*np.eye(64,dtype=np.int64))
assert K.sum()==64 and np.trace(K)==0
assert np.sum(K-np.diag(np.diag(K)))//2==32
out=dict(status='exact PASS',field='F3[a]/(a^2+1)',nonzero_squares=sorted(squares),
         C10=C.tolist(),H20=H.tolist(),checks=['C10 symmetric','C10^2=9I','H20 full sign symmetric','H20^2=20I',
         'H64 symmetric','sum H64=64','trace H64=0','hollow H64 all-one energy=32'],
         family='b=2^r, s=64b, 5s=320b; H_(5s)=H20 tensor H16 tensor H_b')
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ['C10','H20']},indent=2))
