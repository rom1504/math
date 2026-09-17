"""Exact finite replay of the all-isotropic-law Hadamard response dual."""

from fractions import Fraction
import itertools
import json

import numpy as np

from paper_discrepancy_2026_09_17_involution_ground import matchings, parity


SEED = 2026091710
rng = np.random.default_rng(SEED)
reports = []
for p in (4, 8):
    n = p*p
    coords = list(itertools.product(range(p), repeat=2))
    walsh = np.array([[(-1)**parity(u&a) for a in range(p)]
                     for u in range(p)], dtype=np.int64)
    fourier = np.array([[(-1)**(parity(s&v)+parity(t&u))
                        for u,v in coords] for s,t in coords], dtype=np.int64)
    assert np.array_equal(fourier@fourier, p*p*np.eye(n, dtype=np.int64))
    laws, sectors, permutations = [], [], []
    for pairs in matchings(tuple(range(p))):
        perm = np.zeros(p, dtype=np.int64)
        for left,right in pairs:
            perm[left],perm[right] = right,left
        permutations.append(perm)
        for sector in (-1,1):
            for eps in itertools.product((-1,1), repeat=p//2):
                z = np.ones(p, dtype=np.int64)
                for (left,right),sign in zip(pairs,eps):
                    z[left],z[right] = sign,sector*sign
                laws.append([walsh[u,perm[v]]*z[v] for u,v in coords])
                sectors.append(sector)
    laws = np.array(laws, dtype=np.int64)
    sectors = np.array(sectors, dtype=np.int64)
    permutations = np.array(permutations)
    assert np.array_equal(laws@fourier, p*sectors[:,None]*laws)
    cap = n*(p+1)//2
    energies = sectors*p*n//2-n//2
    assert set(cap-np.abs(energies)) == {0,n}
    if p == 4:
        probes = np.array([(1,)+s for s in itertools.product((-1,1), repeat=n-1)], dtype=np.int64)
    else:
        probes = rng.choice((-1,1), size=(600,n)).astype(np.int64)
        # Exceptional identity chirps and genuine eigenwords are included.
        probes = np.concatenate((probes,laws[:120],walsh.reshape(1,n),-walsh.reshape(1,n)))
    coefficients = np.einsum('huv,ua->hva', probes.reshape(-1,p,p),walsh)
    diagonal = np.diagonal(coefficients, axis1=1, axis2=2)
    d = np.sum(diagonal**2,axis=1)
    w = p**3-d
    assert np.array_equal(np.sum(coefficients**2,axis=(1,2)), np.full(len(probes),p**3))
    response_numerators = np.sum(np.abs(probes@laws.T),axis=1)
    # Exact squared version of f >= W/[p(p-1)sqrt(6C_p)].
    lhs = 18*(p-2)*p*p*(p-1)**2*response_numerators**2
    rhs = (p-3)*w*w*len(laws)**2
    assert np.all(lhs >= rhs)
    identity_law = np.array([
        (walsh*np.array(eps)[None,:]).reshape(n)
        for eps in itertools.product((-1,1),repeat=p)],dtype=np.int64)
    assert np.array_equal(identity_law@fourier,p*identity_law)
    identity_response = np.abs(probes@identity_law.T).mean(axis=1)
    assert np.all(3*identity_response**2 >= d-1e-10)
    c_value = 2+p/(p-3)
    aa = p*p/((p-1)*np.sqrt(6*c_value))
    dd = p**1.5/np.sqrt(3)
    weight = aa/(aa+dd)
    constant = aa*dd/(aa+dd)
    mixed_response = (1-weight)*response_numerators/len(laws)+weight*identity_response
    assert np.all(mixed_response >= constant-1e-10)
    count = len(permutations)
    for index in rng.choice(len(probes), size=min(250,len(probes)),replace=False):
        b = coefficients[index]
        m = np.array([sum(int(b[v,perm[v]])**2 for v in range(p))
                      for perm in permutations],dtype=np.int64)
        assert int(m.sum())*(p-1) == count*int(w[index])
        # E M² <= p² C_p E M, clearing denominators.
        assert int(m@m)*(p-3) <= p*p*(3*p-6)*int(m.sum())
    # Every row of F is a physical sign word, and uniform rows are isotropic.
    nu_words = fourier
    assert np.array_equal(nu_words.T@nu_words,n*np.eye(n,dtype=np.int64))
    nu_b = np.einsum('huv,ua->hva',nu_words.reshape(-1,p,p),walsh)
    nu_d = np.sum(np.diagonal(nu_b,axis1=1,axis2=2)**2,axis=1)
    assert int(nu_d.sum()) == n*p*p
    chirp = walsh.reshape(n)
    assert np.all(laws@chirp == 0)
    cp = Fraction(3*p-6,p-3)
    reports.append(dict(p=p,n=n,law_atoms=len(laws),matching_count=count,
                        physical_sign_queries=len(probes),
                        all_projective_queries_exhausted=(p==4),
                        Cp=str(cp), squared_cover_constant=str(1/(24*cp)),
                        exact_dual_squared_slack_min=int((lhs-rhs).min()),
                        exact_exceptional_chirp_overlap=0,
                        identity_law_atoms=len(identity_law),
                        all_sign_laws_pointwise_constant=float(constant),
                        minimum_tested_mixed_response=float(mixed_response.min())))

print(json.dumps(dict(status='PASS',seed=SEED,cases=reports,
    assertions=['both exact Boolean eigensectors',
                'exact cap and nearlevel deficits',
                'column Parseval and quadratic dual for every tested query',
                'matching first and second moments',
                'isotropic cancellation of the quadratic correction',
                'exceptional chirp invalidates uncorrected pointwise constant'],
    scope='Uniform proof is analytic; finite tests do not certify original minimizer applicability'),indent=2))
