"""Bounded exact phase-cube constraint generation for one Paley12 matrix."""
import ctypes,json,time
from pathlib import Path
import numpy as np
root=Path(__file__).resolve().parents[1]
A=np.array(json.loads((root/'computations/results/extension_nested_m11_to_12.json').read_text())['parent_matrix'],dtype=np.int16)
n=12;k=2048;ix=np.arange(k);X=np.column_stack((np.ones(k,dtype=np.int16),1-2*((ix[:,None]>>np.arange(11))&1))).astype(np.int16)
e=np.ascontiguousarray(np.sum((X@A)*X,axis=1)//2,dtype=np.int8)
res={j*j%11 for j in range(1,11)};H=np.ones((12,12),dtype=np.int16)
H[1:,1:]=-np.array([[0 if i==j else(1 if(j-i)%11 in res else -1) for j in range(11)]for i in range(11)])-np.eye(11,dtype=np.int16)
assert np.array_equal(H@H.T,12*np.eye(12,dtype=np.int16))
K=np.ascontiguousarray(abs(X@H@X.T),dtype=np.int8)
lib=ctypes.CDLL(str(root/'tmp/flatify_construct_2026_09_07_bridge12_eval.so'));ev=lib.worst_pair
ev.argtypes=[ctypes.c_void_p,ctypes.c_void_p]+[ctypes.c_int]*4;ev.restype=ctypes.c_uint64
start=time.monotonic();records=[];budget=180
for target in [52,54,56,58,60,62,64]:
  target_feasible=False
  for pol in [1,-1]:
    alive=np.ones((k,k),dtype=bool);constraints=[];witness=None
    while alive.any() and time.monotonic()-start<budget:
      candidate=int(alive.argmax());p,q=divmod(candidate,k)
      value=int(ev(K.ctypes.data,e.ctypes.data,k,p,q,pol));cap=value>>32;i,j=divmod(value&0xffffffff,k)
      if cap<=target:
        bridge=X[p][:,None]*H*X[q][None,:]
        witness=dict(cap=cap,row_phase=p,column_phase=q,polarity=pol,bridge_matrix=bridge.tolist(),parent_matrix=np.block([[A,bridge],[bridge.T,pol*A]]).tolist());break
      alive &= abs(e[ix^i][:,None]+pol*e[ix^j][None,:])+K[i,j]<=target
      constraints.append([i,j,int(K[i,j])])
      if len(constraints)%100==0:print(json.dumps(dict(target=target,polarity=pol,constraints=len(constraints),remaining=int(alive.sum()),elapsed=time.monotonic()-start)),flush=True)
    status='FEASIBLE' if witness else('INFEASIBLE ALL PHASES' if not alive.any() else 'INCOMPLETE')
    rec=dict(target=target,polarity=pol,status=status,constraints=constraints,remaining=int(alive.sum()),witness=witness)
    records.append(rec);print(json.dumps({key:rec[key] for key in ['target','polarity','status','remaining']}),flush=True)
    if witness:target_feasible=True
    if time.monotonic()-start>=budget:break
  if target_feasible or time.monotonic()-start>=budget:break
result=dict(status='bounded exact constraint-generation; only declared fixed Paley row/column order',elapsed=time.monotonic()-start,records=records,hadamard=H.tolist(),child_matrix=A.tolist())
(root/'computations/results/flatify_construct_2026_09_07_bridge12_phase_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
