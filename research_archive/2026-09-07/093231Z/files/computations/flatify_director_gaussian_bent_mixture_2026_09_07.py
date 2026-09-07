"""Exploratory optimization of a proved-form Gaussian/bent sector MGF.

Floating-point exploration only; no asymptotic certificate is claimed.
"""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize


def logcosh(x):
    return np.logaddexp(x, -x)-np.log(2.)


def objective(v, b, entropy_coefficient):
    t=np.exp(v[0]); a=2*t+np.exp(v[1]); lam=(a-1)/2
    if lam<0:
        return 10.+100*(-lam)
    g=1-b
    value=lam*g-g*g/4*np.log(a*a-4*t*t)
    value+=g*b/2*(-np.log(a)+4*t*t/a)
    value+=b*b/4*logcosh(2*np.sqrt(2)*t)
    return (entropy_coefficient*np.log(2)+value)/(2*t)


def main():
    rows=[]
    for b in np.linspace(0,1,21):
        row={"bent_fraction":float(b),"status":"numerical only"}
        for name,ent in [("full_cube",1.),("counted",1-b+(11/32)*b)]:
            fits=[minimize(objective,[np.log(t),np.log(1.)],args=(b,ent),
                           method="Nelder-Mead",options={"maxiter":3000})
                  for t in [1.,4.,8.,16.,32.]]
            best=min(fits,key=lambda z:z.fun)
            row[name]={"q":float(best.fun),"parameters":best.x.tolist(),
                       "success":bool(best.success)}
        rows.append(row)
    out={"status":"uncertified floating-point exploration", "rows":rows}
    dest=Path("computations/results/flatify_director_gaussian_bent_mixture_2026_09_07.json")
    dest.write_text(json.dumps(out,indent=2)+"\n")
    for r in rows:
        print(r["bent_fraction"],r["full_cube"]["q"],r["counted"]["q"])


if __name__=="__main__":
    main()
