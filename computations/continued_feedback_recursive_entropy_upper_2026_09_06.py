"""Concave parity-chord relaxation of the asymmetric Bellman tree.

Solver values are diagnostics until their dual is independently enclosed.
Unlike policy maximization, this optimizes a mathematical UPPER relaxation.
The relaxation can be too weak; its value is not the true Bellman value.
"""
import argparse
import json
import math
import heapq
from fractions import Fraction
from pathlib import Path
import numpy as np
import cvxpy as cp
from scipy.special import logsumexp
from continued_convergence_asymmetric_bellman_2026_09_06 import Problem


def binary_entropy(q):
    return -sum(x*math.log(x) for x in (q,1-q) if x>0)


def rational_log_interval(value,scale=10**12):
    """Exact outward rational enclosure, with fixed denominator after rounding."""
    value=Fraction(value)
    if value<=0:raise ValueError('positive log argument required')
    power=0
    while value<1:value*=2;power-=1
    while value>=2:value/=2;power+=1
    def unit_interval(x):
        z=(x-1)/(x+1)
        partial=2*sum((z**(2*j+1))/Fraction(2*j+1) for j in range(25))
        remainder=2*z**51/(51*(1-z*z))
        return partial,partial+remainder
    lo,hi=unit_interval(value);l2,u2=unit_interval(Fraction(2))
    if power>=0:lo+=power*l2;hi+=power*u2
    else:lo+=power*u2;hi+=power*l2
    floor=(lo.numerator*scale)//lo.denominator
    ceiling=-((-hi.numerator*scale)//hi.denominator)
    return Fraction(floor,scale),Fraction(ceiling,scale)


def rational_entropy_interval(probabilities):
    lower=Fraction(0);upper=Fraction(0)
    for probability in probabilities:
        q=Fraction(probability)
        if not q:continue
        lo,hi=rational_log_interval(q)
        lower-=q*hi;upper-=q*lo
    return lower,upper


def parity_falsifier(max_depth=4):
    p=Fraction(15,16);counts=[15,2,15];denominator=32
    hroot=rational_entropy_interval(Fraction(c,denominator) for c in counts)
    log2=rational_log_interval(2);penalty=[Fraction(0),Fraction(0)]
    records=[]
    for depth in range(1,max_depth+1):
        new=[0]*(2*len(counts)-1)
        for i,a in enumerate(counts):
            for j,b in enumerate(counts):new[i+j]+=a*b
        counts=new;denominator*=denominator
        odd=Fraction(sum(counts[1::2]),denominator)
        hbit=rational_entropy_interval([odd,1-odd])
        penalty[0]+=hbit[0]/2;penalty[1]+=hbit[1]/2
        hleaf=rational_entropy_interval(Fraction(c,denominator) for c in counts)
        lo=p*log2[0]-hroot[1]+hleaf[0]/2-penalty[1]
        hi=p*log2[1]-hroot[0]+hleaf[1]/2-penalty[0]
        records.append({'depth':depth,'constant_lower_rational':str(lo),
            'constant_upper_rational':str(hi),'constant_lower':float(lo),'constant_upper':float(hi),
            'shared_odd_mass':str(odd)})
    return {'p':'15/16','records':records,
      'proved_statement':'parity-relaxed exponent >= listed lower constant + t(1-sqrt(p)) for every t>=0',
      'method':'iid pair policies, identity signed leaf transport, exact rational logarithm enclosures'}


def run(p,t,depth,boxes=None,solver="CLARABEL"):
    problem=Problem(p,t,depth)
    x=cp.Variable(problem.size, nonneg=True)
    constraints=[problem.constraints@x==problem.rhs]
    objective=-binary_entropy(p)-p*math.log(2)
    odd_expressions=[]
    for index,node in enumerate(problem.nodes):
        q=x[node['slice']];d=node['depth'];classes,M,P,Q,mult=node['maps']
        odd=np.array([(i+j)%2 for i,j,s in classes],dtype=float)@q
        odd_expressions.append(odd)
        if boxes is not None:
            lo,hi=boxes[index]
            constraints += [odd>=lo,odd<=hi]
            if hi>lo:
                slope=(binary_entropy(hi)-binary_entropy(lo))/(hi-lo)
                lower=binary_entropy(lo)+slope*(odd-lo)
            else:
                lower=binary_entropy(lo)
            objective -= (2.**(-d-1))*lower
        if d==depth-1:
            for output in (P,Q):
                weights=output@q
                n=output.shape[0]
                a=np.arange(n)/math.sqrt(p*2**depth)
                logk=logsumexp(np.stack([-t*(a[:,None]-a[None,:])**2,
                    -t*(a[:,None]+a[None,:])**2]),axis=0)-math.log(2)
                eta=cp.Variable((n,n), nonneg=True)
                constraints += [cp.sum(eta,axis=1)==weights,cp.sum(eta,axis=0)==weights]
                # C_signed/2 = H(eta)/2 + <log K,eta>/2 + (1-p0) log 2.
                objective += 2.**(-depth)*(cp.sum(cp.entr(eta))/2
                    +cp.sum(cp.multiply(logk,eta))/2+(1-weights[0])*math.log(2))
    task=cp.Problem(cp.Maximize(objective),constraints)
    options={'max_iter':300,'tol_gap_abs':1e-8,'tol_feas':1e-9,'tol_gap_rel':1e-8} if solver=='CLARABEL' else {}
    task.solve(solver=solver,verbose=False,**options)
    return {'p':p,'t':t,'depth':depth,'boxes':boxes,'status':task.status,
       'upper_relaxation_solver_value':task.value,
       'certificate_exponent_relaxation':task.value+p*math.log(2)+t*(1-math.sqrt(p)),
       'odd_masses':[float(v.value) for v in odd_expressions],
       'policy':None if x.value is None else x.value.tolist(),
       'qualification':'concave upper relaxation; numerical conic value not yet an interval-certified bound'}


class ExactEntropyBoxes:
    """Keep all pair dependences; chord only intermediate -H marginals."""
    def __init__(self,p,t,depth):
        self.problem=Problem(p,t,depth)
        self.p=p;self.t=t;self.depth=depth
        pr=self.problem
        self.x=cp.Variable(pr.size,nonneg=True)
        constraints=[pr.constraints@self.x==pr.rhs]
        self.profiles=[];self.profile_nodes=[];self.profile_slices=[];self.profile_weights=[]
        size=0
        for node in pr.nodes:
            if node['depth']==0:continue
            profile=node['maps'][1]@self.x[node['slice']]
            self.profiles.append(profile);self.profile_nodes.append(node)
            self.profile_slices.append(slice(size,size+profile.size));size+=profile.size
            self.profile_weights.extend([2.**(-node['depth'])]*profile.size)
        self.profile_vector=cp.hstack(self.profiles)
        self.lo=cp.Parameter(size,nonneg=True);self.hi=cp.Parameter(size,nonneg=True)
        self.slope=cp.Parameter(size);self.intercept=cp.Parameter(size)
        constraints += [self.profile_vector>=self.lo,self.profile_vector<=self.hi]
        objective=-binary_entropy(p)-p*math.log(2)
        objective += cp.sum(cp.multiply(np.array(self.profile_weights),
            cp.multiply(self.slope,self.profile_vector)+self.intercept))
        # Each non-root signed marginal contributes -(1-p0) log 2.
        for node,profile in zip(self.profile_nodes,self.profiles):
            objective -= 2.**(-node['depth'])*(1-profile[0])*math.log(2)
        for node in pr.nodes:
            q=self.x[node['slice']];d=node['depth'];classes,M,P,Q,mult=node['maps']
            factor=2.**(-d)
            if d<depth-1:
                objective += factor*.5*(cp.sum(cp.entr(q))+np.log(mult)@q)
            else:
                for output in (P,Q):
                    weights=output@q;n=output.shape[0]
                    out_index=np.argmax(output,axis=0)
                    # H(signed pair)-H(absolute child): a conditional entropy.
                    conditional_pair=-cp.sum(cp.rel_entr(q,weights[out_index]))+np.log(mult)@q
                    a=np.arange(n)/math.sqrt(p*2**depth)
                    logk=logsumexp(np.stack([-t*(a[:,None]-a[None,:])**2,
                        -t*(a[:,None]+a[None,:])**2]),axis=0)-math.log(2)
                    eta=cp.Variable((n,n),nonneg=True)
                    constraints += [cp.sum(eta,axis=1)==weights,cp.sum(eta,axis=0)==weights]
                    denominator=cp.reshape(weights,(n,1),order='C')@np.ones((1,n))
                    conditional_transport=-cp.sum(cp.rel_entr(eta,denominator))
                    objective += factor*.25*(conditional_pair+conditional_transport
                        +cp.sum(cp.multiply(logk,eta)))
        self.task=cp.Problem(cp.Maximize(objective),constraints)
        # True feasible profile ranges come from an LP, not the nonlinear search.
        from scipy.optimize import linprog
        map_rows=[]
        for node in self.profile_nodes:
            for row in node['maps'][1]:
                v=np.zeros(pr.size);v[node['slice']]=row;map_rows.append(v)
        self.profile_map=np.array(map_rows)
        lower=[];upper=[]
        for row in self.profile_map:
            a=linprog(row,A_eq=pr.constraints,b_eq=pr.rhs,bounds=(0,None),method='highs')
            b=linprog(-row,A_eq=pr.constraints,b_eq=pr.rhs,bounds=(0,None),method='highs')
            lower.append(max(0,a.fun-1e-10));upper.append(min(1,-b.fun+1e-10))
        self.initial_lo=np.array(lower,dtype=float);self.initial_hi=np.array(upper,dtype=float)

    def solve(self,lo=None,hi=None):
        lo=self.initial_lo if lo is None else np.array(lo)
        hi=self.initial_hi if hi is None else np.array(hi)
        def xlogx(x):
            return np.where(x>0,x*np.log(np.maximum(x,1e-300)),0.)
        slope=np.divide(xlogx(hi)-xlogx(lo),hi-lo,out=np.zeros_like(lo),where=hi>lo)
        intercept=xlogx(lo)-slope*lo
        self.lo.value=lo;self.hi.value=hi;self.slope.value=slope;self.intercept.value=intercept
        self.task.solve(solver='CLARABEL',max_iter=250,tol_gap_abs=1e-8,tol_feas=1e-9,tol_gap_rel=1e-8)
        profile=None if self.x.value is None else self.profile_map@self.x.value
        gap=None if profile is None else (slope*profile+intercept-xlogx(profile))*self.profile_weights
        return {'status':self.task.status,'value':self.task.value,
          'exponent':self.task.value+self.p*math.log(2)+self.t*(1-math.sqrt(self.p)),
          'lo':lo.tolist(),'hi':hi.tolist(),'profile':None if profile is None else profile.tolist(),
          'chord_gap':None if gap is None else gap.tolist()}

    def branch_probe(self,budget):
        """A bounded diagnostic run, NOT a verified branch certificate."""
        first=self.solve();heap=[(-first['value'],0,first)];counter=1;records=[]
        for iteration in range(budget):
            if not heap:break
            negative,_,record=heapq.heappop(heap)
            gaps=np.array(record['chord_gap'])
            index=int(np.argmax(gaps))
            lo=np.array(record['lo']);hi=np.array(record['hi'])
            split=(lo[index]+hi[index])/2
            for lower_half in (True,False):
                child_lo=lo.copy();child_hi=hi.copy()
                if lower_half:child_hi[index]=split
                else:child_lo[index]=split
                # Exact simplex implications strengthen all other chords.
                for repeat in range(3):
                    for sl in self.profile_slices:
                        l=child_lo[sl].copy();h=child_hi[sl].copy()
                        child_lo[sl]=np.maximum(l,1-h.sum()+h)
                        child_hi[sl]=np.minimum(h,1-l.sum()+l)
                if np.any(child_lo>child_hi+1e-9):continue
                try:
                    child=self.solve(np.maximum(child_lo,0),np.maximum(child_hi,child_lo))
                except cp.error.SolverError:
                    print(json.dumps({'solver_failure_at':iteration,'index':index}),flush=True)
                    continue
                if child['profile'] is None:continue
                heapq.heappush(heap,(-child['value'],counter,child));counter+=1
            if iteration%10==0 or iteration==budget-1:
                item={'iteration':iteration+1,'boxes':len(heap),'largest_relaxed_exponent':heap[0][2]['exponent']}
                records.append(item);print(json.dumps(item),flush=True)
        return {'p':self.p,'t':self.t,'depth':self.depth,'iterations':budget,
            'records':records,'top_box':heap[0][2] if heap else None,
            'qualification':'bounded floating-point branch probe only; failures/rounding not enclosed, so no certified global upper claim'}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--p',type=float,default=15/16)
    parser.add_argument('--t',type=float,default=4)
    parser.add_argument('--depth',type=int,default=3)
    parser.add_argument('--box-width',type=float,default=0.)
    parser.add_argument('--exact-boxes',action='store_true')
    parser.add_argument('--branch-budget',type=int,default=0)
    parser.add_argument('--parity-falsifier',action='store_true')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.parity_falsifier:
        result=parity_falsifier(4)
        print(json.dumps(result),flush=True)
        args.output.write_text(json.dumps(result,indent=2)+'\n')
        return
    if args.exact_boxes:
        program=ExactEntropyBoxes(args.p,args.t,args.depth)
        result=program.branch_probe(args.branch_budget) if args.branch_budget else program.solve()
        print(json.dumps(result),flush=True)
        args.output.write_text(json.dumps(result,indent=2)+'\n')
        return
    result=run(args.p,args.t,args.depth)
    print(json.dumps({k:v for k,v in result.items() if k!='policy'}),flush=True)
    if args.box_width:
        boxes=[(max(0,q-args.box_width),min(1,q+args.box_width)) for q in result['odd_masses']]
        second=run(args.p,args.t,args.depth,boxes)
        print(json.dumps({k:v for k,v in second.items() if k!='policy'}),flush=True)
        result={'unrestricted':result,'local_box_only':second}
    args.output.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    main()
