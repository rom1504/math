"""Bounded actual-optimal12+12 Hadamard bridge probe; exact accepted caps."""
import ctypes,itertools,json,time
from pathlib import Path
import numpy as np

root=Path(__file__).resolve().parents[1]
A=np.array(json.loads((root/'computations/results/extension_nested_m11_to_12.json').read_text())['parent_matrix'],dtype=np.int16)
n=12;size=1<<(n-1);ix=np.arange(size,dtype=np.int32)
X=np.column_stack((np.ones(size,dtype=np.int16),1-2*((ix[:,None]>>np.arange(n-1))&1))).astype(np.int16)
energies=np.sum((X@A)*X,axis=1)//2
assert max(abs(energies))==18
e=np.ascontiguousarray(energies,dtype=np.int8)
residues={i*i%11 for i in range(1,11)}
Q=np.array([[0 if i==j else(1 if (j-i)%11 in residues else -1) for j in range(11)] for i in range(11)],dtype=np.int16)
H=np.ones((12,12),dtype=np.int16);H[1:,1:]=-Q-np.eye(11,dtype=np.int16)
assert np.array_equal(H@H.T,12*np.eye(12,dtype=np.int16))
lib=ctypes.CDLL(str(root/'tmp/flatify_construct_2026_09_07_bridge12_eval.so'))
ev=lib.evaluate;ev.argtypes=[ctypes.c_void_p,ctypes.c_void_p]+[ctypes.c_int]*5;ev.restype=ctypes.c_uint64
rng=np.random.default_rng(9071029);start=time.monotonic();records=[];best=None;tested=0;restarts=0
output=root/'computations/results/flatify_construct_2026_09_07_bridge12_probe.json'
def save():
    result=dict(status='bounded phase/permutation search; accepted caps exact; rejected candidates can be cutoff-censored',
                child_order=n,child_cap=18,child_source='extension_nested_m11_to_12.json',
                rowregular_target=float(36*np.sqrt(23/11)),asymptotic_target=float(36*np.sqrt(2)),
                tested=tested,restarts=restarts,elapsed=time.monotonic()-start,best=best,records=records)
    output.write_text(json.dumps(result,indent=2)+'\n')

while time.monotonic()-start<120 and restarts<100:
    rp=rng.permutation(n) if restarts else np.arange(n);cp=rng.permutation(n) if restarts else np.arange(n)
    hp=H[rp][:,cp];K=np.ascontiguousarray(abs(X@hp@X.T),dtype=np.int8)
    p=int(rng.integers(size));q=int(rng.integers(size));pol=1 if restarts%2==0 else -1
    def score(p,q,pol,cut=200):
        global tested
        tested+=1;v=int(ev(K.ctypes.data,e.ctypes.data,size,p,q,pol,cut));return(v>>32,v&0xffffffff)
    current=score(p,q,pol)
    for step in range(40):
        candidate=(current,p,q,pol)
        neighbors=[(p^(1<<bit),q,pol) for bit in range(n-1)]+[(p,q^(1<<bit),pol) for bit in range(n-1)]+[(p,q,-pol)]
        rng.shuffle(neighbors)
        for pp,qq,ss in neighbors:
            sc=score(pp,qq,ss,current[0])
            if sc<candidate[0]:candidate=(sc,pp,qq,ss)
        if candidate[0]>=current:break
        current,p,q,pol=candidate
    restarts+=1
    row=dict(restart=restarts,cap=current[0],maximizer_count=current[1],row_phase=p,column_phase=q,polarity=pol)
    records.append(row)
    if best is None or current<(best['cap'],best['maximizer_count']):
        bridge=X[p][:,None]*hp*X[q][None,:]
        parent=np.block([[A,bridge],[bridge.T,pol*A]])
        direct=int(np.max(abs(energies[:,None]+pol*energies[None,:])+abs(X@bridge@X.T)))
        assert direct==current[0]
        best=dict(**row,row_permutation=rp.tolist(),column_permutation=cp.tolist(),bridge_matrix=bridge.tolist(),parent_matrix=parent.tolist())
        print(json.dumps(dict(improvement=row,elapsed=time.monotonic()-start)),flush=True);save()
    if best['cap']<=50:break
save();print(json.dumps(dict(final_cap=best['cap'],tested=tested,restarts=restarts,elapsed=time.monotonic()-start)),flush=True)
