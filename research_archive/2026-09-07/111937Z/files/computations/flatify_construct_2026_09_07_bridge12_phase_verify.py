"""Replay finite phase exclusion certificates independently of search code."""
import json,subprocess
from pathlib import Path
import numpy as np
root=Path(__file__).resolve().parents[1]
p=json.loads((root/'computations/results/flatify_construct_2026_09_07_bridge12_phase_certificate.json').read_text())
A=np.array(p['child_matrix'],dtype=np.int16);H=np.array(p['hadamard'],dtype=np.int16)
k=2048;ix=np.arange(k);X=np.column_stack((np.ones(k,dtype=np.int16),1-2*((ix[:,None]>>np.arange(11))&1))).astype(np.int16)
e=np.sum((X@A)*X,1)//2;assert max(abs(e))==18;out=[]
for r in p['records']:
  if r['target']==60:
    alive=np.ones((k,k),dtype=bool)
    for i,j,b in r['constraints']:
      assert abs(int(X[i]@H@X[j]))==b
      alive &= abs(e[ix^i,None]+r['polarity']*e[ix^j][None,:])+b<=60
    assert not alive.any()
    out.append(dict(target=60,polarity=r['polarity'],status='PASS full phase cube excluded',constraint_count=len(r['constraints'])))
  if r['target']==62 and r['witness']:
    B=r['witness']['parent_matrix'];data='24\n'+'\n'.join(' '.join(map(str,row)) for row in B)+'\n'
    replay=json.loads(subprocess.run([str(root/'tmp/flatify_construct_2026_09_07_gray')],input=data,text=True,capture_output=True,check=True).stdout)
    assert replay['cap']==62
    out.append(dict(target=62,polarity=r['polarity'],status='PASS full-parent exhaustive Gray replay',replay=replay))
(root/'computations/results/flatify_construct_2026_09_07_bridge12_phase_verify.json').write_text(json.dumps(dict(status='PASS independent replay',records=out),indent=2)+'\n')
print(json.dumps(out))
