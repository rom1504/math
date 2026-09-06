"""Finite falsifiers for a proposed recursive-Hadamard entropy inequality.

No cap is computed. Numerical Sinkhorn evaluations are diagnostics, not
interval-certified inequalities. The saved input determines each test.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.special import logsumexp
from scipy.optimize import minimize_scalar, minimize


def self_cost(points, weights, t):
    points = np.asarray(points, dtype=float)
    if points.ndim == 1:
        points = points[:, None]
    weights = np.asarray(weights, dtype=float)
    keep = weights > 1e-15
    points, weights = points[keep], weights[keep]
    weights = weights / weights.sum()
    logp = np.log(weights)
    distance = ((points[:, None, :] - points[None, :, :])**2).sum(axis=2)
    logk = -t * distance
    u = np.zeros(len(weights))
    for iteration in range(1000):
        new = (u - logsumexp(logk + logp[None, :] + u[None, :], axis=1)) / 2
        if np.max(np.abs(new-u)) < 2e-13:
            u = new
            break
        u = new
    coupling = np.exp(logp[:,None]+logp[None,:]+u[:,None]+u[None,:]+logk)
    error = np.max(np.abs(coupling.sum(axis=1)-weights))
    return float(2*np.dot(weights,u)), float(error)


def entropy(weights):
    weights=np.asarray(weights)
    weights=weights[weights>0]
    return float(-np.dot(weights,np.log(weights)))


def folded_phi(points, weights, t):
    signed_points=[]; signed_weights=[]
    for point,weight in zip(points,weights):
        if weight <= 1e-15: continue
        if point == 0:
            signed_points.append(0.0);signed_weights.append(weight)
        else:
            signed_points.extend([-point,point]);signed_weights.extend([weight/2,weight/2])
    cost,error=self_cost(signed_points,signed_weights,t)
    return -cost/2,error


def ternary_first(p,t):
    delta=1-p; a=1/np.sqrt(p)
    def objective(u):
        theta=np.array([[delta-u,u],[u,p-u]])
        product=np.outer([delta,p],[delta,p])
        keep=theta>0
        divergence=float(np.sum(theta[keep]*np.log(theta[keep]/product[keep])))
        weights=[1-p/2-1.5*u,2*u,(p-u)/2]
        phi,error=folded_phi([0,a/np.sqrt(2),a*np.sqrt(2)],weights,t)
        return phi-divergence/2,error
    optimum=minimize_scalar(lambda u:-objective(u)[0],bounds=(0,min(p,delta)),
      method='bounded',options={'xatol':1e-14})
    choices=[(u,*objective(u)) for u in [0,min(p,delta),optimum.x]]
    u,value,error=max(choices,key=lambda row:row[1])
    h=entropy([delta,p/2,p/2]); gamma=1-np.sqrt(p)
    condensation=[]
    for r in [0,1,2,3,4,8]:
        phi,err=folded_phi([0,1/np.sqrt(p)],[delta,p],(2**r)*t)
        lower=-(1-2**(-r))*h+(2**(-r))*phi
        condensation.append({'depth':r,'Bellman_lower':lower,
          'certificate_exponent_lower':p*np.log(2)+lower+t*gamma,'marginal_error':err})
    return {'p':p,'t':t,'pairing_cross_mass':u,'first_Bellman_value':value,
      'first_certificate_exponent':p*np.log(2)+value+t*gamma,
      'maximum_marginal_error':error,'condensation_paths':condensation,
      'qualification':'one-level symmetry/concavity reduction requires written proof; condensation paths are lower bounds only'}


def one_level(values,weights,t):
    values=np.asarray(values); weights=np.asarray(weights)
    keep=weights>1e-12;values=values[keep];weights=weights[keep];weights/=weights.sum()
    n=len(values); pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    def details(off):
        theta=np.diag(weights.copy())
        for v,(i,j) in zip(off,pairs):
            theta[i,j]=theta[j,i]=v;theta[i,i]-=v;theta[j,j]-=v
        if np.min(theta)<-1e-9:return None
        theta=np.maximum(theta,0)
        out={}
        for i in range(n):
            for j in range(n):
                for val in [abs(values[i]+values[j])/np.sqrt(2),abs(values[i]-values[j])/np.sqrt(2)]:
                    key=round(float(val),12)
                    out[key]=out.get(key,0)+theta[i,j]/2
        ov=np.array(sorted(out));ow=np.array([out[v] for v in ov]);ow/=ow.sum()
        product=np.outer(weights,weights);positive=theta>1e-15
        divergence=float(np.sum(theta[positive]*np.log(theta[positive]/product[positive])))
        phi,error=folded_phi(ov,ow,t)
        return phi-divergence/2,divergence,theta,ov,ow,error
    def constraints(off):
        diagonal=weights.copy()
        for v,(i,j) in zip(off,pairs):diagonal[i]-=v;diagonal[j]-=v
        return diagonal
    def objective(off):
        result=details(off)
        return 1e4 if result is None else -result[0]
    if not pairs:return details([])
    starts=[np.array([weights[i]*weights[j] for i,j in pairs]),np.zeros(len(pairs))]
    candidates=[]
    for start in starts:
        answer=minimize(objective,start,method='SLSQP',bounds=[(0,min(weights[i],weights[j])) for i,j in pairs],
          constraints=[{'type':'ineq','fun':constraints}],options={'ftol':2e-11,'maxiter':300})
        result=details(answer.x)
        if result is not None:candidates.append(result)
    return max(candidates,key=lambda row:row[0])


def greedy_policy(p,t,depth):
    values=np.array([0,1/np.sqrt(p)]);weights=np.array([1-p,p])
    records=[];cost=0
    for level in range(depth):
        result=one_level(values,weights,t)
        value,divergence,theta,ov,ow,error=result
        cost+=divergence/2
        phi,_=folded_phi(ov,ow,t)
        lower=phi-cost
        records.append({'level':level+1,'input_values':values.tolist(),'input_weights':weights.tolist(),
          'theta':theta.tolist(),'output_values':ov.tolist(),'output_weights':ow.tolist(),
          'divergence':divergence,'total_cost':cost,'Bellman_lower':lower,
          'certificate_exponent_lower':p*np.log(2)+lower+t*(1-np.sqrt(p)),
          'marginal_error':error})
        values,weights=ov,ow
    return {'p':p,'t':t,'records':records,'status':'selected symmetric greedy policies; lower bounds only, not cap certificates'}


def condensation_policy(p,t,depth,grid):
    qs=np.linspace(0,p,grid+1)
    amplitude=np.sqrt((2**depth)/p)
    current=np.array([folded_phi([0,amplitude],[1-q,q],t)[0] for q in qs])
    choices=[]
    for level in range(depth-1,-1,-1):
        next_values=np.empty(grid+1);splits=np.empty(grid+1,dtype=int)
        for i,q in enumerate(qs):
            if i==0:next_values[i]=0;splits[i]=0;continue
            fractions=np.arange(i+1)/i
            hs=np.array([entropy([z,1-z]) for z in fractions])
            entries=.5*(current[:i+1]+current[i::-1])-.5*entropy([q,1-q])-.5*q*np.log(2)+.5*q*hs
            splits[i]=int(np.argmax(entries));next_values[i]=entries[splits[i]]
        choices.append(splits.tolist());current=next_values
    value=float(current[-1])
    return {'p':p,'t':t,'depth':depth,'grid':grid,'Bellman_lower':value,
      'certificate_exponent_lower':p*np.log(2)+value+t*(1-np.sqrt(p)),
      'splits_bottom_up':choices,'status':'exactly admissible rational-grid condensation policies; Sinkhorn numerical leaf values; lower bound only'}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--trials', type=int, default=3000)
    ap.add_argument('--seed', type=int, default=9061100)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--ternary',action='store_true')
    ap.add_argument('--p',type=float,default=15/16)
    ap.add_argument('--t',type=float,default=4)
    ap.add_argument('--greedy-depth',type=int,default=0)
    ap.add_argument('--condensation-depth',type=int,default=0)
    ap.add_argument('--grid',type=int,default=512)
    args = ap.parse_args()
    if args.condensation_depth:
        result=condensation_policy(args.p,args.t,args.condensation_depth,args.grid)
        args.output.write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps({key:value for key,value in result.items() if key!='splits_bottom_up'}),flush=True)
        return
    if args.greedy_depth:
        result=greedy_policy(args.p,args.t,args.greedy_depth)
        args.output.write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(result),flush=True)
        return
    if args.ternary:
        result=ternary_first(args.p,args.t)
        args.output.write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(result),flush=True)
        return
    rng = np.random.default_rng(args.seed)
    best = None
    for trial in range(args.trials):
        n = int(rng.integers(2,5))
        a = np.sort(rng.normal(size=n))
        b = np.sort(rng.normal(size=n))
        joint = rng.dirichlet(np.full(n*n, float(rng.choice([.05,.2,1,3])))).reshape(n,n)
        t = float(np.exp(rng.uniform(-3,3)))
        points = np.array([(x,y) for x in a for y in b])
        fj, ej = self_cost(points,joint.ravel(),t)
        fa, ea = self_cost(a,joint.sum(axis=1),t)
        fb, eb = self_cost(b,joint.sum(axis=0),t)
        record = {'trial':trial,'a':a.tolist(),'b':b.tolist(),'joint':joint.tolist(),
          't':t,'F_joint':fj,'F_a':fa,'F_b':fb,'violation':fj-fa-fb,
          'maximum_marginal_error':max(ej,ea,eb)}
        if best is None or record['violation'] > best['violation']:
            best = record
    result = {'seed':args.seed,'trials':args.trials,'best':best,
      'status':'diagnostic for F(joint) <= F(first marginal)+F(second marginal)'}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__ == '__main__':
    main()
