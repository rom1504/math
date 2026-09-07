"""Small exact enumeration checks for positive homogeneous permutation pressure."""
import itertools
import json
import math
from scipy.optimize import minimize_scalar

K = ((1.0, .25), (.25, .75))
out = []
for n, count in [(3,1),(4,1),(5,2)]:
    d = n-1
    p = count/d
    neighbors = [[j for j in range(n) if j != i] for i in range(n)]
    positions = [{j:a for a,j in enumerate(row)} for row in neighbors]
    choices = [tuple(int(a in ones) for a in range(d))
               for ones in itertools.combinations(range(d),count)]
    total = 0.0
    for rows in itertools.product(choices, repeat=n):
        weight = 1.0
        for i in range(n):
            for j in range(i+1,n):
                weight *= K[rows[i][positions[i][j]]][rows[j][positions[j][i]]]
        total += weight
    logz = math.log(total)-n*math.log(len(choices))
    def objective(z):
        gamma = (1-2*p+z,p-z,p-z,z)
        kernel = (K[0][0],K[0][1],K[1][0],K[1][1])
        return sum(v*(math.log(k)-math.log(v)) for v,k in zip(gamma,kernel) if v>0)
    lo,hi=max(0,2*p-1),p
    opt=minimize_scalar(lambda z:-objective(z),bounds=(lo,hi),method='bounded',
                        options={'xatol':1e-14})
    best=max(objective(lo),objective(hi),objective(opt.x))
    upper=n*d/2*best-n*math.log(len(choices))
    entropy=-p*math.log(p)-(1-p)*math.log(1-p)
    pressure=.5*best-entropy
    assert logz <= upper+1e-10
    out.append(dict(n=n,count=count,configurations=len(choices)**n,
                    exact_logz=logz,finite_entropy_upper=upper,
                    entropy_bound_slack=upper-logz,
                    logz_over_n2=logz/n**2,limiting_profile_pressure=pressure))
print(json.dumps(out,indent=2))
