"""Exact finite audit of the corrected Walsh eigenlaws and their responses."""

import itertools
import json
import math
from pathlib import Path

import numpy as np

from paper_discrepancy_2026_09_17_involution_ground import matchings, parity


SEED = 2026091719
rng = np.random.default_rng(SEED)
reports = []
matching_cf_tests = 0
matching_mgf_tests = 0
identity_block_tests = 0
portable_kernel_tests = 0
general_kernel_witness_tests = 0
sharpness_cases = []
mixed_tail_tests = 0
centered_increment_tests = 0
rare_jump_cases = []
for p in (4, 8):
    n = p*p
    chi = np.array([[(-1)**parity(u&v) for v in range(p)]
                    for u in range(p)], dtype=np.int64)
    coords = list(itertools.product(range(p), repeat=2))
    f = np.array([[(-1)**(parity(u&t)+parity(v&s))
                   for s,t in coords] for u,v in coords], dtype=np.int64)
    a = f-np.eye(n,dtype=np.int64)
    pm = {-1: [], 1: []}
    all_matchings = list(matchings(tuple(range(p))))
    for pairs in all_matchings:
        pi = np.zeros(p,dtype=np.int64)
        for v,w in pairs:
            pi[v],pi[w] = w,v
        for sigma in (-1,1):
            for eps in itertools.product((-1,1),repeat=p//2):
                z = np.ones(p,dtype=np.int64)
                for (v,w),sign in zip(pairs,eps):
                    z[v],z[w] = sign,sigma*sign
                pm[sigma].append((chi[:,pi]*z[None,:]).reshape(n))
    pm = {s:np.array(words,dtype=np.int64) for s,words in pm.items()}
    ident = np.array([(chi*np.array(eps)[None,:]).reshape(n)
                      for eps in itertools.product((-1,1),repeat=p)],dtype=np.int64)
    m = len(pm[1])
    mi = len(ident)
    covplus_numer = pm[1].T@pm[1]
    covminus_numer = pm[-1].T@pm[-1]
    covid_numer = ident.T@ident
    eye = np.eye(n,dtype=np.int64)
    assert np.array_equal((p-1)*covminus_numer,m*((p-1)*eye-a))
    # Clear both atomic denominators in the corrected positive pole.
    assert np.array_equal((p-1)*mi*covplus_numer+2*m*covid_numer,
                          m*mi*((p+1)*eye+a))
    mixed_numer = mi*(p-1)*(covplus_numer+covminus_numer)+2*m*covid_numer
    assert np.array_equal(mixed_numer,2*p*m*mi*eye)
    # Deterministic identity-label scheduling has exactly isotropic aggregate.
    cov0 = (covplus_numer+covminus_numer)/(2*m)
    covid = covid_numer/mi
    for q in (p,2*p,5*p):
        r = q//p
        assert np.allclose((q-r)*cov0+r*covid,q*eye)
    for sigma in (-1,1):
        assert np.array_equal(pm[sigma]@f,sigma*p*pm[sigma])
        assert np.all(pm[sigma].sum(axis=0)==0)
    assert np.array_equal(ident@f,p*ident)
    all_law=np.concatenate([pm[-1],pm[1]])
    for _ in range(30):
        xx,yy=rng.choice((-1,1),size=(2,n))
        ww=xx-yy
        bb=ww.reshape(p,p).T@chi
        coeff=np.array([[bb[v,w]+sigma*bb[w,v]
                        for v,w in itertools.combinations(range(p),2)] for sigma in (-1,1)])
        big=float(np.abs(coeff).max())
        mass=float(np.max(np.sum(coeff*coeff,axis=1)))
        difference=np.abs(all_law@xx)-np.abs(all_law@yy)
        centered=difference-difference.mean()
        linear=all_law@ww
        for lam in (0.1/p,0.3/p,0.8/p):
            lhs=float(np.mean(np.exp(lam*centered)))
            middle=float(np.mean(np.cosh(2*lam*linear)))
            assert lhs<=middle+1e-9
            if big>0:
                upper=math.exp(math.e*mass/(p*big*big)*(math.cosh(2*lam*big)-1))
                assert middle<=upper+1e-8
            else:
                assert np.all(centered==0)
            centered_increment_tests+=1
    for _ in range(40):
        ww=rng.normal(size=n)
        bb=ww.reshape(p,p).T@chi
        coeff=np.array([[bb[v,w]+sigma*bb[w,v]
                        for v,w in itertools.combinations(range(p),2)] for sigma in (-1,1)])
        big=float(np.abs(coeff).max())
        mass=float(np.max(np.sum(coeff*coeff,axis=1)))
        assert mass<=2*p*float(ww@ww)+1e-9
        assert big<=2*np.max(np.sum(np.abs(ww.reshape(p,p)),axis=0))+1e-9
        vals=all_law@ww
        for lam in (0.1/p,0.5/p,1.0/p):
            logbound=math.e*mass/(p*big*big)*(math.cosh(lam*big)-1)
            assert np.mean(np.exp(lam*vals))<=math.exp(logbound)+1e-9
            center=np.abs(vals)-np.mean(np.abs(vals))
            assert np.mean(np.exp(lam*center))<=math.exp(2*logbound)+1e-9
            mixed_tail_tests+=1
    xx0=pm[1][0].copy()
    xx1=xx0.reshape(p,p).copy()
    v,w=all_matchings[0][0]
    xx1[:,[v,w]]*=-1
    xx1=xx1.reshape(n)
    diff=xx1-xx0
    assert np.count_nonzero(diff)==2*p
    assert np.array_equal(f@xx1,p*xx1)
    vals=all_law@diff
    assert set(vals)=={-4*p,0,4*p}
    for val in (-4*p,4*p):
        assert int(np.sum(vals==val))*4*(p-1)==len(all_law)
    assert np.isclose(np.mean(vals*vals),8*p*p/(p-1))
    rare_jump_cases.append(dict(p=p,hamming_distance=2*p,jump=4*p,
                                each_jump_probability=1/(4*(p-1))))
    if p==4:
        probes = np.array([(1,)+word for word in itertools.product((-1,1),repeat=n-1)],
                          dtype=np.int64)
    else:
        probes = np.concatenate([rng.choice((-1,1),size=(1000,n)),
                                 pm[-1][:250],pm[1][:250],ident[:50]]).astype(np.int64)
    # Chunk to avoid an unnecessarily large dense temporary matrix.
    max_jensen_slack = 0.0
    max_eigen_response = {-1:0.0,1:0.0}
    for start in range(0,len(probes),256):
        x = probes[start:start+256]
        e = np.einsum('bi,ij,bj->b',x,a,x)/2
        plus = np.abs(x@pm[1].T).mean(axis=1)
        minus = np.abs(x@pm[-1].T).mean(axis=1)
        identity = np.abs(x@ident.T).mean(axis=1)
        resp = (p-1)/(2*p)*(plus+minus)+identity/p
        rhs = ((p+1)/(2*p)*np.sqrt(np.maximum(0,n+2*e/(p+1)))
               +(p-1)/(2*p)*np.sqrt(np.maximum(0,n-2*e/(p-1))))
        assert np.all(resp <= rhs+1e-10)
        max_jensen_slack = max(max_jensen_slack,float(np.max(rhs-resp)))
        for sigma in (-1,1):
            mask = np.all(x@f==sigma*p*x,axis=1)
            if np.any(mask):
                max_eigen_response[sigma] = max(max_eigen_response[sigma],float(resp[mask].max()/p))
        b = np.einsum('buv,ua->bva',x.reshape(-1,p,p),chi)
        for sigma in (-1,1):
            eigen = np.all(x@f==sigma*p*x,axis=1)
            assert np.all(b[eigen]==sigma*b[eigen].transpose(0,2,1))
        assert np.all(np.sum(b*b,axis=2)==p*p)
    reports.append(dict(p=p,n=n,matching_law_atoms=m,identity_atoms=mi,
                        tested_queries=len(probes),all_projective_queries=(p==4),
                        maximum_eigen_response=max_eigen_response,
                        largest_jensen_slack=max_jensen_slack))
    edge_list = list(itertools.combinations(range(p),2))
    edge_index = {e:j for j,e in enumerate(edge_list)}
    matching_indices = np.array([[edge_index[tuple(sorted(e))] for e in pairs]
                                 for pairs in all_matchings])
    test_queries = probes[rng.choice(len(probes),size=80,replace=False)]
    fourier_coeff = np.einsum('buv,ua->bva',test_queries.reshape(-1,p,p),chi)/p
    for b in fourier_coeff:
        for sigma in (-1,1):
            c = np.array([b[v,w]+sigma*b[w,v] for v,w in edge_list])
            assert np.max(np.abs(c))<=2+1e-12
            assert c@c<=2*p+1e-12
            degree = np.zeros(p)
            for (v,w),cc in zip(edge_list,c):
                degree[v]+=cc*cc
                degree[w]+=cc*cc
            maxdegree = degree.max()
            for t in (0.1,0.25,0.5,1.0,2.0):
                phi = np.prod(np.cos(t*c[matching_indices]),axis=1).mean()
                ww = 1-np.cos(t*c)
                psi = math.exp(-float(ww.sum())/(p-1))
                d = maxdegree*t*t/2
                finite_bound = 30/p*(d+d*d)*math.exp(2*d)
                assert abs(phi-psi)<=finite_bound+1e-12
                assert psi>=math.exp(-t*t*float(c@c)/(2*(p-1)))-1e-12
                mgf = np.prod(np.cosh(t*c[matching_indices]),axis=1).mean()
                mgf_bound = math.exp(math.e/2*(math.cosh(2*t)-1))
                assert mgf<=mgf_bound+1e-10
                matching_cf_tests+=1
                matching_mgf_tests+=1
    # Also audit the generic weighted-degree CF lemma independently of Fourier arrays.
    for _ in range(80):
        c = rng.normal(size=len(edge_list))/math.sqrt(p)
        degree = np.zeros(p)
        for (v,w),cc in zip(edge_list,c):
            degree[v]+=cc*cc
            degree[w]+=cc*cc
        for t in (0.1,0.5,1.0):
            phi = np.prod(np.cos(t*c[matching_indices]),axis=1).mean()
            psi = math.exp(-float(np.sum(1-np.cos(t*c)))/(p-1))
            d = degree.max()*t*t/2
            assert abs(phi-psi)<=30/p*(d+d*d)*math.exp(2*d)+1e-12
            matching_cf_tests+=1
    if p==4:
        for _ in range(20):
            r=3
            ee=rng.choice((-1,1),size=(p,r))
            physical=np.stack([(chi*ee[:,j][None,:]).reshape(n) for j in range(r)],axis=1)
            fullcap=np.max(np.sum(np.abs(probes@physical),axis=1))
            smallwords=np.array(list(itertools.product((-1,1),repeat=p)))
            reducedcap=p*np.max(np.sum(np.abs(smallwords@ee),axis=1))
            assert fullcap==reducedcap
            identity_block_tests+=1
    # Portability to nonsymmetric switched/permuted Hadamard kernels.
    for _ in range(12):
        uu=chi[rng.permutation(p)][:,rng.permutation(p)].copy()
        uu*=rng.choice((-1,1),size=p)[:,None]
        uu*=rng.choice((-1,1),size=p)[None,:]
        ff=np.array([[uu[u,t]*uu[s,v] for s,t in coords] for u,v in coords])
        assert np.array_equal(ff,ff.T)
        assert np.array_equal(ff@ff,n*eye)
        assert np.all(np.diag(ff)==1)
        pairs=all_matchings[rng.integers(len(all_matchings))]
        pi=np.zeros(p,dtype=np.int64)
        for v,w in pairs:
            pi[v],pi[w]=w,v
        for sigma in (-1,1):
            zz=np.ones(p,dtype=np.int64)
            for v,w in pairs:
                zz[v]=rng.choice((-1,1))
                zz[w]=sigma*zz[v]
            hh=(uu[:,pi]*zz[None,:]).reshape(n)
            assert np.array_equal(ff@hh,sigma*p*hh)
        portable_kernel_tests+=1
    for _ in range(40):
        uu=rng.choice((-1,1),size=(p,p)).astype(np.int64)
        ff=np.array([[uu[u,t]*uu[s,v] for s,t in coords] for u,v in coords])
        xx=uu.reshape(n)
        gram=uu.T@uu
        assert int(xx@ff@xx)==int(np.trace(gram@gram))
        assert int(xx@ff@xx)>=p**3
        general_kernel_witness_tests+=1

kappa = math.sqrt(2/math.pi)
window=[]
for c in (0.433322111664,0.46,0.48,0.493608094,0.5):
    target = 1.5*c/kappa
    tau = math.sqrt(1-(2*target*target-1)**2)/(2*c)
    window.append(dict(c=c,required_tau=tau,energy_feasibility_ceiling=2*c,
                       spectral_norm_ceiling=1/tau))
for r in (2,4,8):
    p=r*r
    uu=np.array([[(-1)**parity(u&v) for v in range(p)] for u in range(p)],dtype=np.int64)
    dd=np.array([(-1)**parity((u//r)&(u%r)) for u in range(p)],dtype=np.int64)
    xx=dd[:,None]*uu
    assert np.array_equal(uu@xx.T@uu,p*xx)
    b=xx.T@uu
    assert np.all(np.abs(b)==r)
    assert np.all(np.diag(b)==r)
    def mu(j):
        return sum(math.comb(j,k)*abs(2*k-j) for k in range(j+1))/2**j
    exact=((p-1)*mu(p//2)+mu(p))/(p*np.sqrt(p))
    sharpness_cases.append(dict(p=p,n=p*p,exact_response=float(exact),
                                limiting_response=1/math.sqrt(math.pi)))
report=dict(status='PASS',seed=SEED,cases=reports,
    matching_characteristic_tests=matching_cf_tests,
    matching_exponential_moment_tests=matching_mgf_tests,
    exact_identity_block_tests=identity_block_tests,
    portable_hadamard_kernel_tests=portable_kernel_tests,
    arbitrary_sign_kernel_witness_tests=general_kernel_witness_tests,
    sharpness_cases=sharpness_cases,
    mixed_tail_tests=mixed_tail_tests,
    centered_increment_tests=centered_increment_tests,
    rare_jump_cases=rare_jump_cases,
    jensen_threshold=6/math.sqrt(145),
    gaussian_response_threshold=math.sqrt(72*math.pi/(256+81*math.pi**2)),
    window=window,
    assertions=['exact both-pole covariances','exact isotropic corrected mixture',
                'exact Boolean eigensectors','uniform finite Jensen response',
                'Fourier symmetry of all tested eigenqueries',
                'finite signed-matching characteristic comparison',
                'all-query uniform exponential moment',
                'exact deterministic aggregate isotropy',
                'exact identity-column full-bridge reduction'])
output=Path('tmp/paper_portfolio_2026_09_17/localization/nonlocal_sign_response_audit.json')
output.parent.mkdir(parents=True,exist_ok=True)
output.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
