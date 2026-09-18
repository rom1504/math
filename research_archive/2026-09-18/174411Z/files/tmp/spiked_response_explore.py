#!/usr/bin/env python3
import importlib.util
import itertools
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "extremal_information/experiments/actual_child_bridge_law_exact.py"
spec = importlib.util.spec_from_file_location("acl", SRC)
acl = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = acl
spec.loader.exec_module(acl)


def row_q(n, y):
    masks = np.arange(1 << n, dtype=np.uint64)
    bits = ((masks[:, None] >> np.arange(n, dtype=np.uint64)) & 1).astype(np.int8)
    b = 1 - 2 * bits
    z = b @ y / np.sqrt(n)
    e = (1 + z*z) / 2
    return e + z, e - z


def response(pressure, m, n, y):
    qp, qm = row_q(n, y)
    tensor = pressure.reshape((1 << n,) * m)
    out = []
    for v in itertools.product((-1, 1), repeat=m):
        arr = tensor
        for vi in v:
            arr = np.tensordot(qp if vi == 1 else qm, arr, axes=(0, 0)) / (1 << n)
        out.append(float(arr))
    return np.asarray(out)


def spectrum_by_degree(values, m):
    # itertools (-1,+1) order; numpy FWHT gives characters of +1-bit masks.
    f = values.copy()
    width = 1
    while width < len(f):
        a = f.reshape(-1, 2 * width)
        left = a[:, :width].copy(); right = a[:, width:].copy()
        a[:, :width] = left + right; a[:, width:] = left - right
        width *= 2
    f /= len(f)
    return [float(np.sum(f[[i for i in range(len(f)) if bin(i).count('1')==k]]**2)) for k in range(m+1)]


def main():
    import mpmath as mp
    mp.mp.dps = 80
    records=[]
    cache={}
    for N in range(4,10):
      m=N//2; n=N-m
      for beta in (0.5,1,2,4,8):
        cls=[]
        for k in (m,n):
          key=(k,beta,N)
          if key not in cache:
            cache[key]=acl.thermal_minimizer_classes(acl.build_signing_space(k), str(beta), N)
          cls.append(cache[key][0])
        for li,lrow in enumerate(cls[0]):
          A=np.asarray(lrow['representative_matrix'],dtype=np.int8)
          for ri,rrow in enumerate(cls[1]):
            D=np.asarray(rrow['representative_matrix'],dtype=np.int8)
            for eps in (-1,1):
              L,_=acl.bridge_pressures(A,D,beta,N,eps)
              best=None
              for ytail in itertools.product((-1,1),repeat=n-1):
                y=np.asarray((1,)+ytail,dtype=np.int8)
                R=response(L,m,n,y)
                rec={"y":y.tolist(),"min":float(R.min()),"max":float(R.max()),"range":float(np.ptp(R)),"sd":float(R.std()),"spectrum":spectrum_by_degree(R,m)}
                if best is None or rec['range']>best['range']: best=rec
              row={"N":N,"split":[m,n],"beta":beta,"left":li,"right":ri,"eps":eps,"best":best,"L_range":float(np.ptp(L))}
              records.append(row); print(row,flush=True)
    print(json.dumps(records))

if __name__ == '__main__': main()
