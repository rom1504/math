# Uniform covariance retained by actual optimized Gibbs laws

Status: proved below; independently reconstructed by `transfer_seeds`.
This is a theorem about actual finite-temperature optimizers. It is NOT a
cross-order recurrence or a claim that variance controls exponential pressure.

## Statement

Let A be a hollow symmetric signing, J=beta A/sqrt(n), and let

    mu_g(s,x) proportional to exp[s sum_(i<j) Jij xi xj + h g s],

where g is standard Gaussian (h=0 is allowed). Choose A globally minimizing
E_g log Z_A(g), before drawing g. All expectations below include g and its
Gibbs law. Let D=diag(d_i) satisfy D>=J and D>=-J in PSD order, with d_i>=0.
Then

    E[xx^T] >= diag( exp(-2 beta^2)/(1+d_i+2 beta^2) ).       (1)

No spectral flatness, bounded operator norm, pure-state assumption, or
pointwise-in-g edge optimality is used. The only optimality used in (1) is
optimal replacement of each entire incident row in the QUENCHED objective.

If Q(A)<=C n^(3/2), the simultaneous Grothendieck diagonal majorant gives
sum d_i<=4 K_G beta C n. Therefore the diagonal weights c_i in (1) satisfy

    sum_i c_i >= n exp(-2 beta^2)/(1+2 beta^2+4 K_G beta C).  (2)

Here any valid real Grothendieck constant may be used, for example
K_G=pi/(2 asinh(1)). Its simultaneous-majorant proof, including the factor
four for Q=max|x^TAx|/2, is reconstructed in
`transfer_director_sparse_restriction_audit_2026_09_06.md`, Section 1.

For fixed beta and h=o(n), an actual optimizer has some uniform finite C:
Jensen gives log Z_A(0)<=E log Z_A(g), while comparison with a fixed bounded-cap
signing gives E log Z_A(g)<=beta C_0 n+(n+1)log2+h sqrt(2/pi).
Thus Q(A)/n^(3/2)<=C_0+(1+1/n)log2/beta+o(1).

## 1. Row optimality supplies average conditional variance

Delete vertex i and let nu_g be the Gibbs law of (s,x_-i) in the remaining
system. Put h_i=sum_(j!=i) Jij xj and C_i(g)=E_(nu_g) cosh h_i.
Integrating xi gives Z_A(g)=2 Z_cavity(g) C_i(g).

Replace the signs incident to i by independent uniform signs. Global
optimality, averaging this legitimate row replacement, and concavity of log
give

    E_g log C_i(g)
      <= E_g E_randomrow log C_i^random(g)
      <= (n-1)log cosh(beta/sqrt(n)) <= beta^2/2.             (3)

The cavity law is independent of these replaced signs. At fixed g,

    v_i(g):=E_(mu_g) sech^2 h_i
       = E_(nu_g) sech h_i / C_i(g) >= C_i(g)^(-2),

using E cosh h_i E sech h_i>=1. Jensen in g and (3) imply

    v_i:=E_g v_i(g) >= exp(-beta^2).                         (4)

This is averaged in g; no unjustified pointwise quenched optimum is asserted.

## 2. Score covariance and an explicit second-order error

Set eta_i=xi-s tanh h_i. Conditional expectation in xi gives

    E[eta_i xj]=1{i=j} v_i,  E eta_i^2=v_i<=1.              (5)

For i!=j write b_j=sum_(k!=i,j) Jjk xk. Conditional-spin integration gives

    E eta_i eta_j
       =-E[s sech^2 h_i {tanh(b_j+Jij)-tanh(b_j-Jij)}/2].

The derivative of sech^2 has absolute value at most 2. Hence, for xi=+-1,

    |{tanh(b+J)-tanh(b-J)}/2-J sech^2(b+Jxi)|<=2J^2.

Let V=diag(sech^2 h_i). The off-diagonal score covariance is consequently
-E[s VJV]+R, with |Rij|<=2Jij^2 and Rii=0. Symmetry holds because this is
the difference of symmetric covariance and Gram matrices. The diagonal of
E[s VJV] is zero. Since D>=+-J and 0<=V<=I,

    -E[s VJV] <= E[VDV] <=D,
    ||R||op <=2 max_i sum_j Jij^2 <=2 beta^2.

Thus

    E[eta eta^T] <= I+D+2 beta^2 I.                         (6)

The joint second-moment matrix of (x,eta) is PSD. Replacing its lower-right
block by the larger positive diagonal matrix in (6) preserves PSD. Taking
the Schur complement and using (4),(5) proves (1). The x law is globally
spin-reversal symmetric, so its second moment is also its covariance.
Jensen for 1/(a+d_i) and the trace majorant prove (2).

## 3. A consequence and its exact limitation

For two independently sampled optimized child laws satisfying (2), and ANY
flat rectangular sign bridge B of size m by n,

    E(x^T B y)^2 = tr(B Sigma_y B^T Sigma_x)
                 >=(sum_i c_i)(sum_j c'_j)>=kappa_beta,C mn. (7)

The first comparison follows successively from PSD order, and the final
equality uses Bij^2=1. More generally any bridge gives sum Bij^2 c_i c'_j.
This excludes a collapse of the actual product-child covariance onto a
bridge-annihilating subspace, without storing its full covariance.

It does NOT prove an extensive logarithmic moment-generating function:
symmetry only immediately gives

    log E exp(t x^T B y) >= log[1+t^2 kappa mn/2].

At t=beta/sqrt(m+n) this is only order log(m+n). A two-point variable of
size sqrt(mn) illustrates why variance alone cannot supply order-(m+n)
pressure. Uniform control under the nonzero bridge tilt would be a further
theorem, not an implication silently drawn from (7).

## Scope and originality assessment

The Schur/score covariance argument is classical matrix Cauchy--Schwarz.
The useful input here is exact whole-row optimality, combined with a
simultaneous diagonal majorant rather than a full operator-norm hypothesis.
This removes a spectral assumption from a concrete statement about actual
optimizer Gibbs laws. It is structural progress, not a solution of the
original convergence problem. No claim of a new general information theory.
