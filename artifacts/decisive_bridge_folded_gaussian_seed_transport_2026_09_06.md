# Folded/unfolded transport equality and a seed-sensitive Gaussian boundary

Date: 2026-09-06. Status: exact variational identities. The director and
bridge agent independently noticed the folding equality in discussion.
This does not identify actual finite-spin source entropy with Gaussian
source entropy, nor prove an actual cap bound from the Gaussian formula.

## 1. Exact folding identity

Let nu be centrally symmetric on R^q, R symmetric orthogonal, and t>0.
For a positive kernel K define

`Psi_nu(K)=sup_(pi coupling(nu,nu)) [integral logK dpi-D(pi||nu tensor nu)]`.

Assume finite second moment, so the Gaussian-distance kernels below have
finite reference expectation. Define b(x,y)=2t x^T R y and

`K_single=exp[-t(||x||²+||y||²)+b]`,
`K_fold=exp[-t(||x||²+||y||²)]cosh b`.

Then `Psi_nu(K_single)=Psi_nu(K_fold)` EXACTLY. There is no log2 loss.
This is a variational equality, not a pointwise kernel equality.

Proof in the folded-to-single direction: independently reverse the signs
of x and y in any folded coupling and average. Marginals are unchanged,
the folded energy is unchanged, and relative entropy decreases. Therefore
restrict to couplings gamma invariant under both reversals. Define

`dpi = (exp(b)/cosh(b)) dgamma`.

This is a coupling of nu with itself: conditional on x, average the two
signs of y to see that its density has mean one, and vice versa. The
density is at most 2. Every independently sign-invariant observable has
the same expectation under pi and gamma. Expanding the relative entropy
of pi using this density shows

`E_pi logK_single-D(pi||nu²)
 = E_gamma logK_fold-D(gamma||nu²)`.

In the reverse direction first simultaneously sign-symmetrize any single
coupling pi, which preserves its energy and decreases entropy. Put
`gamma=(pi+reflection_y pi)/2`; it is independently sign-symmetric.
Let pi_star=(exp(b)/cosh(b))gamma. The same expansion, now with the
additional nonnegative divergence D(pi||pi_star), gives

`E_pi logK_single-D(pi||nu²)
 = E_gamma logK_fold-D(gamma||nu²)-D(pi||pi_star)`.

Thus single cannot exceed folded. Infinite-entropy couplings do not improve
the supremum; finite entropy and second moment justify each expression by
truncation if necessary. Points x=0 or y=0 cause no exception: b=0 there
and the density multiplier is 1. No choice of projective representatives
or a noncompact Sinkhorn existence theorem is required.

## 2. Explicit Gaussian formula

For nu=N(0,V), V positive semidefinite, define
`W=V^(1/2) R V^(1/2)` on the support of V, with singular values s_j.
For s>=0 let rho(s) be the unique number in [0,1) solving
`2t s=rho/(1-rho²)` and put

`h_t(s)=t s rho(s)+(1/4)log(1-rho(s)²)`.

Then half the transport pressure is

`(1/2)Psi_nu(K_fold)=(1/2)Psi_nu(K_single)
   =-t Tr(V)+sum_j h_t(s_j)`.

To verify, normalize the nonzero Gaussian coordinates so a coupling has
cross-correlation C. Gaussian maximum entropy gives
`I(X;Y)>=-(1/2)log det(I-C C^T)`; equality is attained by the jointly Gaussian
coupling with that cross-correlation. Singular-value trace inequality reduces
the maximization to independent scalars rho_j. Their stationary equations
are precisely the equation above. One can alternatively diagonalize the
symmetric W and choose C with eigenvalues sign(lambda_j(W))*rho(s_j),
which attains the bound. The case s=0 contributes zero; singular V is
handled by restricting to its support. The folded optimizer is generally
a mixture of the Gaussian optimizer with its y-reflection, not Gaussian.

The scalar notation in the campaign has `h_t(s)=t s+g_t(s)`.

## 3. A literal seed dependence in the Gaussian joint-column family

For a symmetric full sign seed A of order k, let
`(R_A X)_ij=A_ij X_ji`. Let X have independent Gaussian columns, each with
row covariance C; thus V=I_(columns) tensor C in column-major convention.

If C=I, then V=I, all k² singular values of W equal 1, and the Gaussian
pressure is seed-independent.

If C=J=11^T, write X=1 z^T with z standard Gaussian in R^k. The orthonormal
embedding E:z->(1/sqrt(k))z^T has V=k E E^T. Direct computation gives

`k E^T R_A E=A`.

Therefore the nonzero spectrum of W is EXACTLY the spectrum of A, and

`(1/2)Psi_(N(0,I tensor J))(K_(R_A))
 = -t k² + sum_(j=1..k) h_t(|lambda_j(A)|)`.

This retains a spectral seed invariant and is not a scalar |A|/entrywise
norm expression. It still depends only on the seed spectrum in this
rank-one-column family, not directly on its Boolean quadratic cap.
The entropy cost of forcing k physical spins to agree is substantial:
instead of k ell independent signs, one has only ell signs, losing
(k-1)ell log2 choices per group before any additional source constraint.
Whether pressure gain offsets that loss in the full union bound requires
the actual source enumerator and normalization; this identity alone does
not establish it.
