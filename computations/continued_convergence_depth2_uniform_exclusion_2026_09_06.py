"""Interval proof excluding the depth <=2 criterion for every p,t.

Each terminal cost uses an explicit self-coupling, so this is a LOWER bound
on a selected Bellman policy. Adaptive boxes cover the compact parameter
region; separate analytic bounds cover its complement.
"""
import argparse
from fractions import Fraction
import json
from pathlib import Path
from mpmath import iv


def law(p):
    d=1-p
    return [d**4+3*p*p*d*d+3*p**4/8,
            4*p*d**3+3*p**3*d,
            3*p*p*d*d+p**4/2,
            p**3*d,p**4/8]


def entropy(p):
    return -p*iv.log(p)-(1-p)*iv.log(1-p)


def condensation(p,t):
    return p*iv.log(2)/8-7*entropy(p)/8+t*(1-iv.sqrt(p))


def coupling_coefficients(pm,tm):
    import numpy as np
    from scipy.special import logsumexp
    p=np.array(law(pm));a=np.arange(5)/(2*np.sqrt(pm))
    logk=logsumexp(np.stack([-tm*(a[:,None]-a[None,:])**2,
                            -tm*(a[:,None]+a[None,:])**2]),axis=0)-np.log(2)
    lp=np.log(p);u=np.zeros(5)
    for iteration in range(500):
        nxt=.5*(u-logsumexp(logk+lp[None,:]+u[None,:],axis=1))
        if max(abs(nxt-u))<2e-13:u=nxt;break
        u=nxt
    q=np.exp(lp[:,None]+lp[None,:]+u[:,None]+u[None,:]+logk)
    coeff=q/(p[:,None]*p[None,:]);np.fill_diagonal(coeff,0)
    return coeff


def evaluate(box,save=False,fixed_coefficients=None):
    pa,pb,ta,tb=box
    p=iv.mpf([str(pa),str(pb)]);t=iv.mpf([str(ta),str(tb)])
    lower_c=condensation(p,t)
    if lower_c.a>0:return float(lower_c.a),None,'condensation'
    pm=(pa+pb)/2;tm=(ta+tb)/2
    prob=law(p);D=2**40
    if fixed_coefficients is None:
        coeff=coupling_coefficients(pm,tm)
        # Scale by an interval-valid bound before rational rounding.
        ratio=1.
        for i in range(5):
            total=iv.mpf(0)
            for j in range(5):
                if i!=j:total+=iv.mpf(str(coeff[i,j]))*prob[j]
            if total.b>1:ratio=min(ratio,1/float(total.b))
        ratio*=1-1e-12
        integers=[[0]*5 for _ in range(5)]
        for i in range(5):
            for j in range(i+1,5):
                integers[i][j]=integers[j][i]=max(0,int(coeff[i,j]*ratio*D))
    else:
        integers=fixed_coefficients
        assert len(integers)==5 and all(len(row)==5 for row in integers)
        for i in range(5):
            assert integers[i][i]==0
            for j in range(5):assert isinstance(integers[i][j],int) and integers[i][j]>=0 and integers[i][j]==integers[j][i]
    remainders=[]
    for i in range(5):
        rem=1-sum(iv.mpf(integers[i][j])/D*prob[j] for j in range(5) if i!=j)
        if rem.a<0:return -1e6,None,'wide'
        remainders.append(rem)
    cost=iv.mpf(0)
    for i in range(5):
        for j in range(5):
            logk=-t*(i-j)**2/(4*p)+iv.log(1+iv.exp(-t*i*j/p))-iv.log(2)
            if i==j:
                r=remainders[i]
                if r.a<=0:return -1e6,None,'wide'
                q=prob[i]*r
                cost+=q*(iv.log(r/prob[i])-logk)
            elif integers[i][j]:
                c=iv.mpf(integers[i][j])/D;q=c*prob[i]*prob[j]
                cost+=q*(iv.log(c)-logk)
    lower_i=p*iv.log(2)-cost/2+t*(1-iv.sqrt(p))
    if lower_i.a>0:
        return float(lower_i.a),integers if save else None,'iid'
    return max(float(lower_c.a),float(lower_i.a)),None,'split'


def complement_checks():
    # g(p)=p log2+log(1-p)/4 is concave and g(0)=0.
    # Therefore positivity at .92 excludes every 0<p<=.92 by the
    # already-proved Gaussian profile obstruction.
    a=iv.mpf('0.92');b=iv.mpf('0.981');c=iv.mpf('0.945')
    gaussian=a*iv.log(2)+iv.log(1-a)/4
    high_p=condensation(b,iv.mpf(0))
    # C_5 is strictly convex; its tangent at c lower-bounds it throughout
    # [.92,.981]. This range is all that remains after the p exclusions.
    derivative=iv.log(2)/8-7*iv.log((1-c)/c)/8-5/(2*iv.sqrt(c))
    large_t=condensation(c,iv.mpf(5))+derivative*(iv.mpf(['0.92','0.981'])-c)
    assert gaussian.a>0 and high_p.a>0 and large_t.a>0
    return {'gaussian_at_p_point92':str(gaussian),'condensation_at_p_point981_t0':str(high_p),
      'convex_tangent_lower_for_t_ge5':str(large_t)}


def run(precision,max_boxes):
    iv.dps=precision;outside=complement_checks()
    stack=[(.92,.981,0.,5.)];leaves=[];minimum=1.;visited=0
    while stack:
        box=stack.pop();visited+=1
        bound,coeff,kind=evaluate(box,save=True)
        if bound>0:
            minimum=min(minimum,bound)
            record={'box':box,'lower_bound':bound,'kind':kind}
            if coeff is not None:record['coupling_coefficients_numerator']=coeff
            leaves.append(record)
        else:
            pa,pb,ta,tb=box
            if (pb-pa)*80>(tb-ta):
                mid=(pa+pb)/2;stack.extend([(pa,mid,ta,tb),(mid,pb,ta,tb)])
            else:
                mid=(ta+tb)/2;stack.extend([(pa,pb,ta,mid),(pa,pb,mid,tb)])
        if visited%1000==0:print(json.dumps({'visited':visited,'pending':len(stack),'leaves':len(leaves),'minimum':minimum}),flush=True)
        if visited>=max_boxes:break
    return {'status':'PASS' if not stack else 'INCOMPLETE','precision':precision,'visited':visited,
      'leaves':len(leaves),'minimum_interval_lower_bound':minimum,'complement':outside,
      'coefficient_denominator':2**40,'cover':leaves,'pending':stack,
      'scope':'criterion excluded for every fixed p in (0,1), t>0 at depth<=2; not an actual ensemble cap obstruction'}


def verify(data,precision):
    iv.dps=precision;outside=complement_checks()
    assert data['status']=='PASS' and not data['pending'] and data['coefficient_denominator']==2**40
    boxes=[];minimum=1.
    for record in data['cover']:
        assert record['kind'] in {'iid','condensation'}
        box=tuple(Fraction(str(x)) for x in record['box']);pa,pb,ta,tb=box
        assert Fraction('0.92')<=pa<pb<=Fraction('0.981') and 0<=ta<tb<=5
        boxes.append(box)
        coeff=record.get('coupling_coefficients_numerator')
        if record['kind']=='iid':assert coeff is not None
        # A condensation box does not need coefficients and must verify
        # before any coefficient generation would be reached.
        if record['kind']=='condensation':
            p=iv.mpf([str(record['box'][0]),str(record['box'][1])])
            t=iv.mpf([str(record['box'][2]),str(record['box'][3])])
            bound=condensation(p,t)
            assert bound.a>0
            minimum=min(minimum,float(bound.a))
        else:
            bound,_,_=evaluate(record['box'],fixed_coefficients=coeff)
            assert bound>0;minimum=min(minimum,bound)
    assert sum((pb-pa)*(tb-ta) for pa,pb,ta,tb in boxes)==Fraction('0.061')*5
    for i,a in enumerate(boxes):
        for b in boxes[i+1:]:
            assert not (max(a[0],b[0])<min(a[1],b[1]) and max(a[2],b[2])<min(a[3],b[3]))
    return {'status':'PASS','boxes':len(boxes),'cover':'exact rational area and pairwise interior-disjointness PASS',
      'minimum_interval_lower_bound':minimum,'precision':precision,'complement':outside}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--precision',type=int,default=35)
    ap.add_argument('--max-boxes',type=int,default=30000);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--verify',action='store_true')
    args=ap.parse_args()
    if args.verify:
        print(json.dumps(verify(json.loads(args.output.read_text()),args.precision),indent=2));return
    result=run(args.precision,args.max_boxes)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['cover','pending']},indent=2))


if __name__=='__main__':main()
