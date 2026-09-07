"""Exact integer covariance checks on the first nonlinear bipartite tree."""
import itertools
import json
import numpy as np

rng=np.random.default_rng(2047)
trials=0
energy_tests=0
for m in range(1,5):
    cross_blocks=(np.array(values,dtype=np.int64).reshape(m,m)
                  for values in itertools.product([-1,1],repeat=m*m)) if m<=3 else (
                  rng.choice([-1,1],size=(m,m)) for _ in range(30))
    for C in cross_blocks:
        A=np.block([[np.zeros((m,m),dtype=int),C],
                    [C.T,np.zeros((m,m),dtype=int)]])
        N=2*m
        triples=list(itertools.combinations(range(N),3))
        idx={t:k for k,t in enumerate(triples)}
        coeff=np.zeros((N,len(triples)),dtype=np.int64)
        for i in range(N):
            for j in range(N):
                if A[i,j]==0:
                    continue
                for k,l in itertools.combinations(range(N),2):
                    if i in (k,l) or j in (k,l):
                        continue
                    coeff[i,idx[tuple(sorted((j,k,l)))]]+=A[i,j]*A[j,k]*A[j,l]
        gram=coeff@coeff.T
        aa=A@A
        for i in range(N):
            for j in range(N):
                if i==j:
                    assert 2*gram[i,j]==m*(m-1)*(m-2)
                elif (i<m)==(j<m):
                    assert 2*gram[i,j]==(m-2)*(m-3)*aa[i,j]
                else:
                    assert gram[i,j]==0
        spins=np.array(list(itertools.product([-1,1],repeat=N)),dtype=np.int64)
        x1=spins@A.T
        monomials=np.array([[int(np.prod(spin[list(t)])) for t in triples]
                            for spin in spins],dtype=np.int64)
        x3=monomials@coeff.T
        # X1=X1_integer/sqrt(m), X3=sqrt(2)*X3_integer/m^(3/2).
        assert np.array_equal(x1.T@x1,(1<<N)*aa)
        assert np.array_equal(x3.T@x3,(1<<N)*gram)
        assert not np.any(x1.T@x3)
        shore_spins=np.array(list(itertools.product([-1,1],repeat=m)),dtype=np.int64)
        bilinear=int(np.max(np.abs(shore_spins@C.T).sum(axis=1)))
        beta=int(np.max(np.abs(spins@A.T).sum(axis=1)))
        cap_twice=int(np.max(np.abs(np.sum(spins*(spins@A.T),axis=1))))
        assert beta==2*bilinear
        assert cap_twice==2*bilinear
        for _ in range(4):
            f=rng.integers(-2,3,size=N,dtype=np.int64)
            h=2-np.abs(f)
            af=A@f
            orientation=np.sign(af).astype(np.int64)
            plus=f+h*orientation
            minus=-f+h*orientation
            assert max(np.max(np.abs(plus)),np.max(np.abs(minus)))<=2
            difference=int(plus@A@plus-minus@A@minus)
            assert difference==4*int(h@np.abs(af))
            assert abs(int(plus@A@plus))<=8*bilinear
            assert abs(int(minus@A@minus))<=8*bilinear
            energy_tests+=1
        trials+=1
print(json.dumps(dict(status='PASS',full_sign_cross_blocks=trials,exhaustive_m=[1,2,3],sampled_m4=30,feasible_rational_energy_tests=energy_tests,arithmetic='exact integer',tests=['X1 Gram','X3 variance','same-shore Q covariance','cross-shore zero covariance','X1-X3 orthogonality','beta(A)=2L(C)','cap(A)=L(C)','feasible mean energy difference'])))
