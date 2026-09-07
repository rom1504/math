# Failed isotropic-ground shortcut; exact finite certificate

2026-09-07. **Falsified finite hypothesis; no asymptotic conclusion.**

Root tested whether cap below the spectral half scale could exclude an
isotropic probability law on the absolute ground states. Such an exclusion
would have supplied a potential route to balancing those states. The generic
claim is FALSE. The stored order-six signing has cap5 and twelve projective
absolute grounds whose uniform covariance is exactly I_6. The accompanying
program enumerates every projective spin, reconstructs the cap independently,
and checks G^T G=12I with integers. No solver optimum claim is needed for
this falsification. The signing is also the recorded exact order-six optimum.

Its cap5 is below both6sqrt(5)/2 and6^(3/2)/2. If G has any isotropic
ground law mu, its factorization norm is maximal: gamma2(G)=sqrt(n).
Indeed for G=LV with column norms of V at most one and row norms of L
at most gamma, the matrix diag(sqrt(mu))G has all n singular values one.
The nuclear/Frobenius inequality gives

    n<=||diag(sqrt(mu))L||F ||V||F<=gamma sqrt(n).

The identity factorization gives the reverse bound. Thus even a small
finite exact minimizer can have maximally difficult balanced-ground data.
This does not prove that a liminf sequence has the same property.

A separate screen of the four normalized atlas representatives attaining
cap10 at order8 found no isotropic absolute-ground law by a floating LP.
That is reconnaissance, not a rational infeasibility certificate and not
an asymptotic theorem. It must not rescue the false general hypothesis.

There is also a useful orientation warning. Isotropy plus local stability
gives E||Ax||1=2Q, while Cauchy--Schwarz gives2Q<=n sqrt(n-1), an UPPER
bound, not a lower bound. Reversing that inequality would falsely turn the
spectral half scale into a necessary cap. All finite data and parameters
are saved by `computations/principle_director_ground_isotropy_probe_2026_09_07.py`.
