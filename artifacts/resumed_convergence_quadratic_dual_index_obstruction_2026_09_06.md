# A precise restriction on quadratic dual-bent index permutations

Date: 2026-09-06. Status: an elementary exact restriction, with a
quantitative Hamming consequence, for the QUADRATIC vectorial
dual-bent realization family. This is not an obstruction to the
profile-level Fourier matching theorem, and is not a bound on R.

## 1. Setup and inverse polar matrices

Let Phi,Psi:F_2^n -> F_2^m, n even, satisfy the exact component
duality identity

    H_n (-1)^(alpha dot Phi)
                  =(-1)^(sigma(alpha) dot Psi),           (1)

for all nonzero alpha, where sigma permutes F_2^m and fixes zero.
Assume every component of Phi is quadratic. Every component of
Psi is then quadratic as well: the dual of a nondegenerate binary
quadratic function is quadratic, and sigma is onto.

Write the alternating polar matrices of the components as

    M_alpha=sum_i alpha_i M_i,
    N_beta =sum_j beta_j N_j.

All M_alpha,N_beta are nonsingular for nonzero indices, because
their corresponding quadratic functions are bent. The exact
duality relation gives

    N_(sigma(alpha))=M_alpha^(-1),       alpha!=0.          (2)

For completeness, if q has polar matrix M, then translation in
its Walsh sum gives

    q*(u+Mh)+q*(u)=q(h)+q(0)+u dot h.

Subtract the same formula at u=0. The polar pairing of q* at
(u,Mh) is u dot h, so its polar matrix is M^(-1). This proves
(2), without importing a classification of quadratic bent maps.

## 2. An invertible bilinear compatibility form

Fix any coordinate r in F_2^n and define the m-by-m binary matrix

    Q_ij=(M_i N_j)_(rr).

Taking the (r,r) entry of M_alpha N_(sigma(alpha))=I gives

    b(alpha,sigma(alpha))=1   for every alpha!=0,
    b(alpha,beta)=alpha^T Q beta.                          (3)

The matrix Q is necessarily nonsingular. A nonzero vector in
its left kernel would contradict (3) directly. Equivalently a
nonzero vector in its right kernel, pulled back by sigma, would
give the same contradiction. Thus every quadratic dual-bent
index permutation satisfies a nondegenerate bilinear edge rule.

This is only a NECESSARY rule. It is much weaker than the full
inverse-linear-space identity (2); an arbitrary permutation
satisfying (3) need not come from any quadratic bent family.

## 3. Exact and quantitative disagreement with linear maps

Let H be a d-dimensional linear subspace of F_2^m, d>=2, and
let L:F_2^m -> F_2^m be any linear map, not necessarily invertible.
Then

    #{alpha in H\{0}: sigma(alpha)!=L alpha}
                                      >=2^(d-2)-1.        (4)

Indeed, in coordinates on H the Boolean polynomial

    p(alpha)=1+b(alpha,L alpha)

has algebraic degree at most two and p(0)=1. Every nonzero
point where p=1 is a disagreement by (3). A nonzero binary
polynomial of degree at most two on d variables has support
size at least 2^(d-2). One elementary proof is the standard
induction: write p=p_0+x_d p_1; if p_1=0 its weight doubles,
and otherwise the sum of the two slice weights is at least
the weight of p_1, to which the degree-one-smaller induction
applies. Subtract the point zero to obtain (4).

In particular sigma cannot agree with any linear map on all
nonzero points of a three-dimensional linear subspace. This
also follows directly because the nonzero indicator on F_2^3
has algebraic degree three, whereas b(alpha,L alpha) has
degree at most two.

For m>=4 no quadratic dual-bent sigma can fix a whole linear
hyperplane pointwise. More quantitatively, it must differ from
the identity on at least 2^(m-3)-1 points of that hyperplane.

## 4. Why this does NOT obstruct the matching norm theorem

The enlarged-model Fourier-matching theorem used permutations
supported on an affine frequency hyperplane a=1 and fixed the
complement a=0 pointwise. Section 3 says that this PARTICULAR
extension is unavailable from quadratic dual-bent indices in
dimension m>=4. It does not say the active matching is
unavailable: both profile Fourier arrays vanish on a=0, so the
permutation on that complement can be changed freely, provided
the active hyperplane is preserved.

Furthermore, even a positive Hamming distance from a desired
index permutation does not imply positive matching cost for
Gaussian point clouds. Many different frequency permutations
can match nearby profile values. The bilinear rule (3) itself
is dense and may be compatible with arbitrarily accurate
profile-level transport. The actual quadratic-family gap is
the much stronger simultaneous inverse-linear-space
realization (2), not the fixed-complement convention.

No assertion is made here about nonquadratic dual-bent maps,
arbitrary Maiorana--McFarland constructions, all available R
operators, or the original convergence problem.
