#!/usr/bin/env python3
"""Finite independent direct-spin audit of cycle repair and near-order bounds.

This is a test of the proof's pointwise norm/deletion/completion steps, not
an exhaustive proof of its universal statement and not a minimizer search.
"""
import json
from functools import lru_cache
from fractions import Fraction
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
SEED=2026091901


@lru_cache(None)
def spins(n):
    masks=np.arange(1<<n,dtype=np.int64)
    return 1-2*((masks[:,None]>>np.arange(n))&1)


def beta(a):
    return int(np.max(np.sum(np.abs(spins(len(a))@a.T),axis=1)))


def q(a):
    x=spins(len(a))
    energies=np.einsum('bi,ij,bj->b',x,a,x)//2
    return int(np.max(np.abs(energies)))


def row_norm(a,indices):
    return int(np.max(np.sum(np.abs(spins(len(a))@a.T)[:,list(indices)],axis=1)))


def twist(a,p,s):
    s=np.asarray(s,dtype=np.int64)
    return a[np.ix_(p,p)]*s[:,None]*s[None,:]


def doubled(a,b,d):
    c=b+np.diag(d)
    return np.block([[a,c],[c,-a]])


def repair(p,removed):
    """Independently follow each surviving vertex until the next survivor."""
    p1=np.arange(len(p),dtype=np.int64)
    for i in range(len(p)):
        if i in removed:continue
        j=int(p[i])
        while j in removed:j=int(p[j])
        p1[i]=j
    assert sorted(p1.tolist())==list(range(len(p)))
    return p1


def bernoulli_audit(a,p,size):
    """Exact expectation over every reservoir; compare radicals by squaring."""
    n=len(a);subsets=(1-spins(n))//2
    subset_norms=np.max(subsets@np.abs(spins(n)@a.T).T,axis=1)
    cardinalities=np.sum(subsets,axis=1)
    denominator=n**n;numerator=0;large_numerator=0;good_numerator=0
    beta_a=beta(a)
    inverse=np.argsort(p)
    # Indicator of p(S) at vertex j is indicator of S at p^{-1}(j).
    image_masks=np.sum(subsets[:,inverse]*(1<<np.arange(n)),axis=1)
    paired_numerator=0
    for mask in range(1<<n):
        k=int(cardinalities[mask]);weight=(2*size)**k*(n-2*size)**(n-k)
        numerator+=int(subset_norms[mask])*weight
        paired=int(subset_norms[mask]+subset_norms[image_masks[mask]])
        paired_numerator+=paired*weight
        if k>=size:
            large_numerator+=weight
            excess=Fraction(paired)-48*Fraction(size*beta_a,n)
            if excess<=0 or excess*excess<=48**2*n**2*size:good_numerator+=weight
    expectation=Fraction(numerator,denominator)
    density=Fraction(2*size,n)
    excess=expectation-density*beta_a
    radical_squared=n**2*2*density*(1-density)*(n-1)
    assert excess<=0 or excess*excess<=radical_squared
    assert paired_numerator==2*numerator
    assert Fraction(large_numerator,denominator)>=Fraction(1,6)
    assert good_numerator>0
    return dict(density=str(density),exact_expected_L=str(expectation),
        exact_linear_bound_term=str(density*beta_a),radical_bound_term_squared=str(radical_squared),
        probability_size_at_least_r=str(Fraction(large_numerator,denominator)),
        probability_joint_good_reservoir=str(Fraction(good_numerator,denominator)),
        all_reservoir_subsets=1<<n)


def audit_case(n,case,rng):
    upper=np.triu(rng.choice([-1,1],size=(n,n)),1)
    a=upper+upper.T
    p=rng.permutation(n);s=rng.choice([-1,1],size=n)
    size=int(rng.integers(1,n//4+1))
    removed=set(map(int,rng.choice(n,size=size,replace=False)))
    kept=sorted(set(range(n))-removed);rr=sorted(removed)
    p1=repair(p,removed);s1=s.copy()
    for i in range(n):
        if p1[i]!=p[i]:s1[i]=rng.choice([-1,1])
    changed=set(map(int,np.flatnonzero((p1!=p)|(s1!=s))))
    pinv_removed={i for i in range(n) if int(p[i]) in removed}
    union=removed|{int(p[i]) for i in removed}
    assert changed<=removed|pinv_removed
    assert {int(p[i]) for i in changed}<=union
    assert {int(p1[i]) for i in changed}<=union
    assert {int(p1[i]) for i in kept}==set(kept)
    assert {int(p1[i]) for i in rr}==removed
    b0=twist(a,p,s);b1=twist(a,p1,s1);difference=b0-b1
    unchanged=sorted(set(range(n))-changed)
    assert not np.any(difference[np.ix_(unchanged,unchanged)])
    lu=row_norm(a,union);lr=row_norm(a,removed)
    lt0=row_norm(b0,changed);lt1=row_norm(b1,changed)
    assert lt0==row_norm(a,{int(p[i]) for i in changed})
    assert lt1==row_norm(a,{int(p1[i]) for i in changed})
    row_difference=row_norm(difference,changed)
    error_beta=beta(difference)
    assert row_difference<=lt0+lt1<=2*lu
    assert error_beta<=3*row_difference<=6*lu
    incident=a.copy();incident[np.ix_(kept,kept)]=0
    assert beta(incident)<=3*lr
    assert 2*q(incident)<=3*lr
    ac=a[np.ix_(kept,kept)]
    assert 0<=q(a)-q(ac)<=q(incident)
    assert 0<=beta(a)-beta(ac)<=beta(incident)
    record=dict(case=case,n=n,removed=rr,A=a.tolist(),p=p.tolist(),s=s.tolist(),
        repaired_p=p1.tolist(),repaired_s=s1.tolist(),changed_columns=sorted(changed),
        union=sorted(union),L_R=lr,L_union=lu,beta_bridge_difference=error_beta,
        L_T_difference=row_difference,L_T_B0=lt0,L_T_B1=lt1,
        Q_incident=q(incident),beta_incident=beta(incident),
        Q_A=q(a),Q_restricted=q(ac),beta_A=beta(a),beta_restricted=beta(ac))
    record['bernoulli_checks']=bernoulli_audit(a,p,size)
    if n<=6:
        d=rng.choice([-1,1],size=n)
        original=doubled(a,b0,d);repaired=doubled(a,b1,d)
        kept_parent=kept+[n+i for i in kept]
        restricted=repaired[np.ix_(kept_parent,kept_parent)]
        local={v:i for i,v in enumerate(kept)}
        pc=np.asarray([local[int(p1[i])] for i in kept])
        assert np.array_equal(restricted,doubled(ac,twist(ac,pc,s1[kept]),d[kept]))
        qo,qp,qc=q(original),q(repaired),q(restricted)
        assert qc<=qp<=qo+error_beta<=qo+6*lu
        # Reverse direction: an independently chosen child twist is extended
        # by the identity on R. This checks the pointwise completion bound,
        # stronger than checking it only at an unknown optimal child twist.
        small_p=rng.permutation(len(kept));small_s=rng.choice([-1,1],size=len(kept))
        small_d=rng.choice([-1,1],size=len(kept))
        extended_p=np.arange(n);extended_s=rng.choice([-1,1],size=n)
        extended_d=rng.choice([-1,1],size=n)
        for i,v in enumerate(kept):
            extended_p[v]=kept[int(small_p[i])]
            extended_s[v]=small_s[i];extended_d[v]=small_d[i]
        small=doubled(ac,twist(ac,small_p,small_s),small_d)
        extended_b=twist(a,extended_p,extended_s)
        extended=doubled(a,extended_b,extended_d)
        bridge_incident=extended_b.copy();bridge_incident[np.ix_(kept,kept)]=0
        assert beta(bridge_incident)<=3*lr
        qs,qe=q(small),q(extended)
        explicit_cost=2*q(incident)+beta(bridge_incident)+size
        assert qe<=qs+explicit_cost<=qs+6*lr+size
        record['parent_checks']=dict(d=d.tolist(),Q_original=qo,Q_repaired=qp,
            Q_repaired_restriction=qc,small_p=small_p.tolist(),small_s=small_s.tolist(),
            small_d=small_d.tolist(),extended_p=extended_p.tolist(),
            extended_s=extended_s.tolist(),extended_d=extended_d.tolist(),
            Q_small=qs,Q_extended=qe,beta_incident_extended_bridge=beta(bridge_incident),
            explicit_completion_cost=explicit_cost,
            full_spin_states_per_order_n_parent=1<<(2*n))
    return record


def main():
    rng=np.random.default_rng(SEED);records=[]
    for n in range(4,11):
        for case in range(15):records.append(audit_case(n,case,rng))
        print(f'PASS order{n}:15 cases',flush=True)
    result=dict(schema='twisted-chiral-near-order-pointwise-audit-v1',random_seed=SEED,
        proof_source='artifacts/twisted_chiral_near_order_2026_09_19.md',
        exact_integer_arithmetic=True,case_count=len(records),parent_case_count=sum('parent_checks' in r for r in records),
        all_assertions_passed=True,scope='Finite independent test of pointwise cycle containment, exact bilinear norms, restriction, and completion bounds; not a proof by exhaustive finite testing and not an F minimization',
        records=records)
    out=ROOT/'computations/results/twisted_chiral_near_order_audit_2026_09_19.json'
    out.write_text(json.dumps(result,indent=2)+'\n');print(out)


if __name__=='__main__':main()
