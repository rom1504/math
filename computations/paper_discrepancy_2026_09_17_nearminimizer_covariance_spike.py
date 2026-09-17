"""Exact full-cube replay of the pinned covariance-spike construction."""
from pathlib import Path
import importlib.util
import itertools
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location("stored",ROOT/"computations/transfer_adversary_minimizer_isotropy_2026_09_06.py")
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
rng=np.random.default_rng(2026091721)
reports=[]
for name,matrix,provenance in module.cases():
    n=len(matrix)
    if n not in (4,6,8):
        continue
    m=12 if n<=6 else 14
    a=np.asarray(matrix,dtype=np.int64)
    old=np.asarray([(1,)+s for s in itertools.product((-1,1),repeat=n-1)],dtype=np.int64)
    new=np.asarray(list(itertools.product((-1,1),repeat=m)),dtype=np.int64)
    old_energy=np.einsum("bi,ij,bj->b",old,a,old)//2
    cap=int(np.abs(old_energy).max())
    assert cap==provenance["minimum_cap_imported"]
    if old_energy.max()!=cap:
        a=-a
        old_energy=-old_energy
    bridge=rng.choice((-1,1),size=(n,m))
    field=old@bridge
    bridge_cap=int(np.abs(field).sum(axis=1).max())
    f=m*(m-1)//2
    gap=f-m//2-bridge_cap
    assert gap>0
    new_sums=new.sum(axis=1)
    new_energy=(new_sums*new_sums-m)//2
    # Projective old words and all new words represent the entire parent cube.
    parent=old_energy[:,None]+new_energy[None,:]+field@new.T
    parent_cap=int(np.abs(parent).max())
    assert cap+f<=parent_cap<=cap+f+bridge_cap
    assert int((-parent).max())<=cap+m//2+bridge_cap
    windows=sorted({0,1,gap-1})
    local=[]
    for window in windows:
        selected=parent_cap-np.abs(parent)<=window
        ii,jj=np.nonzero(selected)
        assert np.all(parent[ii,jj]>0)
        assert np.all(new_sums[jj]**2>=m*m-2*(bridge_cap+window))
        assert np.all(old_energy[ii]>=cap-bridge_cap-window)
        local.append({"window":window,"projective_nearcode_size":len(ii),
                      "minimum_new_sum_square":int(np.min(new_sums[jj]**2)),
                      "covariance_rayleigh_lower_numerator":int(np.min(new_sums[jj]**2)),
                      "covariance_rayleigh_lower_denominator":m})
    reports.append({"core":name,"core_order":n,"appended_order":m,
                    "parent_order":n+m,"core_cap":cap,"bridge_cap":bridge_cap,
                    "parent_cap":parent_cap,"positive_sector_gap_lower":gap,"windows":local})
print(json.dumps({"status":"PASS exact complete parent cubes","global_core_optima_imported":True,"cases":reports},indent=2))
