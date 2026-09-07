"""Gaussian boundary diagnostics, not an actual finite-spin cap certificate."""
import itertools
import json
import math
import numpy as np

k=4
t=4.0
p=31/32
seeds={"all_plus":np.ones((k,k)),"one_negative_edge":np.ones((k,k))}
seeds["one_negative_edge"][0,1]=seeds["one_negative_edge"][1,0]=-1

def ht(s):
    s=np.asarray(s)
    z=2*t*s
    rho=2*z/(1+np.sqrt(1+4*z*z))
    return t*s*rho+.25*np.log1p(-rho*rho)

def reflection(A):
    R=np.zeros((k*k,k*k))
    for i in range(k):
        for j in range(k):
            R[i+k*j,j+k*i]=A[i,j]
    return R

def entropy(alpha):
    probs=np.full(2**k,(1-alpha)/2**k)
    probs[0]+=alpha/2
    probs[-1]+=alpha/2
    return -sum(x*math.log(x) for x in probs if x>0)

rows=[]
for alpha in [0,.1,.25,.5,.75,1.0]:
    C=(1-alpha)*np.eye(k)+alpha*np.ones((k,k))
    eig,U=np.linalg.eigh(C)
    sqrtC=(U*np.sqrt(np.maximum(eig,0)))@U.T
    sqrtV=np.kron(np.eye(k),sqrtC)
    row={"alpha":alpha,"physical_row_entropy":entropy(alpha),
         "source_entropy_per_d2":p*k*entropy(alpha)}
    for name,A in seeds.items():
        singular=np.abs(np.linalg.eigvalsh(sqrtV@reflection(A)@sqrtV))
        pressure=-t*k*k+float(sum(ht(singular)))
        row[name]={"gaussian_pressure_per_d2":pressure,
                   "entropy_plus_gaussian_pressure":p*k*entropy(alpha)+pressure}
        if alpha==1:
            seedformula=-t*k*k+float(sum(ht(np.abs(np.linalg.eigvalsh(A)))))
            assert abs(seedformula-pressure)<1e-10
    rows.append(row)
assert abs(rows[0]["all_plus"]["gaussian_pressure_per_d2"]-
           rows[0]["one_negative_edge"]["gaussian_pressure_per_d2"])<1e-12
print(json.dumps({"warning":"diagnostic only; Gaussian replacement not actual source entropy theorem",
                  "k":k,"t":t,"p":p,"rows":rows},indent=2))
