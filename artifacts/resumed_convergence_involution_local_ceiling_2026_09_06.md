# An original-class ceiling for fixed equivariant first-order rules

Date: 2026-09-06. Status: theorem-level reduction, with its algorithmic scope
spelled out. This is not a bound on the Boolean optimum of a Hadamard
matrix, and it does not settle the original convergence problem.

## 1. Statement

Let `B_n` be deterministic real symmetric matrices such that

    B_n^2=I,     max_ij |(B_n)_ij|=n^(-1/2+o(1)).                 (1)

Consider a fixed finite computation using multiplication by `B_n` and
coordinatewise globally Lipschitz functions of previously computed
vectors and finitely many independent Gaussian/Rademacher seed vectors.
All functions, dimensions, and iteration counts are independent of `n`.
Assume signed-permutation equivariance: conjugating `B_n` by a signed
permutation and applying the same signed permutation to the seed vectors
applies that signed permutation to every vector passed to a matrix
multiplication and to the final vector. In the usual all-spin-register
description this says that each coordinate update is jointly odd in its
spin-type arguments. Even local gates are allowed inside such an odd
update. A coordinatewise map need NOT be equivariant merely because it is
the same at every coordinate.

If the output `m_n` belongs to `[-1,1]^n`, then

    limsup_n E |m_n^T B_n m_n|/(2n) <= sqrt(15)/8
                                      =0.4841229182759271... .    (2)

The same conclusion holds for each fixed finite Lipschitz approximation
to such a rule, and therefore for a rule which is approximable by them
in normalized mean-square output error, uniformly along the matrix
sequence and its Haar comparator.

In particular, (2) applies to normalized symmetric Hadamard matrices.
Their diagonal-deleted hollow signings differ by operator norm `o(1)`
after normalization. Fixed globally Lipschitz computations and their
normalized energies are stable under this perturbation, so the same
algorithmic ceiling holds for these actual hollow signings.

There is NO assertion that every injective graph-polynomial response is
automatically a first-order rule. A finite injective-tree formula requires
an actual normalized-L2 approximation by such a rule, or a separate
traffic argument. Nor is there any assertion about depths, dimensions,
functions, or seed counts growing with `n`, matrix-dependent initial
vectors, arbitrary non-equivariant rules, or optimizing over the entire
Boolean cube.

## 2. Primary theorem and verification of its matrix hypotheses

The primary source is Wang--Zhong--Fan, *Universality of approximate
message passing algorithms and tensor networks*,
[arXiv:2206.13037v5](https://arxiv.org/pdf/2206.13037), specifically
Definition 2.6, Proposition 2.7(b)(2), Theorem 2.8, and Lemma 2.14 with
Remark 2.15. These statements were read directly from the primary PDF.

Proposition 2.7(b)(2) says that for `W_n=Pi_n B_n Pi_n^T`, with `Pi_n`
a uniform independent signed permutation, it suffices that every fixed
power satisfy, for every fixed epsilon>0 eventually,

    max_i |(B_n^r)_ii - Tr(B_n^r)/n| < n^(-1/2+epsilon),
    max_(i!=j) |(B_n^r)_ij| < n^(-1/2+epsilon).                   (3)

The resulting generalized-invariant ensemble has the same limiting
diagonal distribution as the orthogonally invariant ensemble with the
same limiting spectrum. Under (1), even powers are exactly `I` and odd
powers are exactly `B_n`, so (3) follows directly (absorb the factor 2
by using epsilon/2 in (1)). Also

    |Tr B_n|/n <= max_i |B_ii|=o(1),

so the spectral law tends to `(delta_1+delta_-1)/2`; the operator norm is
exactly one. No averaged-to-uniform regularization, puncturing, or
approximate-involution theorem is used here.

Theorem 2.8 supplies state evolution for the orthogonally invariant AMP
prescription, with continuous nonlinearities of polynomial growth which
are Lipschitz in the state arguments, and nonsingular state covariance
matrices. The mixed Gaussian/Rademacher seed law satisfies its
Assumption 2.1: all moments exist and polynomials are dense in L2.
Lemma 2.14 and Remark 2.15 separately give universality of polynomial
first-order empirical moments with arbitrary fixed correction constants.
Moment universality alone is NOT being used to assert convergence of
an arbitrary clipped polynomial; high-degree polynomial laws need not be
moment determinate.

## 3. Passing from the AMP prescription to Lipschitz first-order rules

A finite first-order computation can be serialized as

    y_t=B u_t,
    u_(t+1)=f_(t+1)(y_1,...,y_t; seed data),

with finite history and globally Lipschitz `f_t`. Introduce the AMP
coordinates

    z_t=y_t-sum_(s<=t) b_ts u_s,

where `b_ts` are the Haar/Bernoulli AMP prescription. Inductively recover
each `y_t` from the available `z` history and already recovered `u`
history. This defines coordinatewise Lipschitz AMP nonlinearities and is
an exact finite-dimensional reparameterization, not an approximation.
The coefficients can be chosen recursively from the Haar state law.

For completeness, covariance degeneracy can be removed before applying
Theorem 2.8. Add `delta eta_t` to the input of the t-th matrix
multiplication, using a fresh independent standard Gaussian seed
`eta_t` which earlier updates do not use. The Bernoulli limiting spectral
law has mean zero, so the AMP diagonal coefficient `b_tt` is zero.
Consequently the t-th AMP residual contains a fresh term
`delta B eta_t`, independent of the preceding residuals conditional on
the matrix and earlier seeds. Since `B` is orthogonal, this term is an
independent standard Gaussian vector times `delta`, even after removing
the conditioning on the matrix. Its cross empirical inner products with
the preceding residuals vanish and its empirical squared norm tends to
`delta^2`. Induction therefore makes every finite state covariance
positive definite. A dummy independent initial Gaussian multiplication
may be inserted if needed to place an arbitrary first seed-based update
after the initialization required by Theorem 2.8.

Apply Theorem 2.8 at fixed `delta>0`. Remove the perturbations afterward:
the finite Lipschitz recursion and `||B||op=1` give

    (E ||m_delta-m||_2^2/n)^(1/2) <= C_algorithm delta.           (4)

Use a final coordinatewise clipping if necessary to keep the perturbed
output feasible; clipping is 1-Lipschitz and fixes the original feasible
output. Equation (4) holds for both the deterministic generalized-
invariant sequence and the Haar comparator, with the same constant.
Thus the universality conclusion extends to the original possibly
degenerate finite Lipschitz rule.

One extra matrix multiplication exposes the energy as the empirical
coordinate test `m_i (B m)_i`. Wasserstein-2 state convergence controls
this quadratic test. If this extra step has degeneracy, the same fresh
noise argument applies, and its energy contribution vanishes with
`delta`. Feasibility and bounded operator norm make the normalized energy
uniformly bounded, so convergence in probability also passes to its
expected absolute value.

Finally, equivariance removes the auxiliary signed permutation exactly:

    m(Pi B Pi^T, seeds)=Pi m(B,Pi^T seeds).

The transformed iid centrally symmetric seed rows have the same joint
law as the original seeds. Hence the energy law for `Pi B Pi^T` is
exactly the energy law for `B`, not merely asymptotically the same.

## 4. Exact Haar union bound and its scope

Let `U=O diag(I_(n/2),-I_(n/2)) O^T` for Haar orthogonal `O`, with `n`
even. For each fixed Boolean vector `x`,

    x^T U x/n = 2V-1,       V~Beta(n/4,n/4).

For each fixed `0<t<1`, integration of the beta density, or its elementary
Laplace bound, gives

    P(|x^T U x|/n>=t) <= exp{(n/4)log(1-t^2)+o(n)}.

Union over at most `2^n` Boolean vectors. Whenever `t>sqrt(15)/4`,

    log 2 + (1/4)log(1-t^2)<0,

and therefore

    P(max_x |x^T U x|/(2n)>sqrt(15)/8+eta) <= exp(-c_eta n)      (5)

for every fixed eta>0 with a nonvacuous threshold. Thus the limsup
expected cap is at most `sqrt(15)/8`. This is an UPPER bound only, not
an exact Haar ground-state value. Balanced dimension can be replaced by
eigenspace ranks with relative dimensions tending to one half: the beta
large-deviation rate converges to the same displayed rate.

For nonzero diagonals, maximizing a quadratic form on `[-1,1]^n` is not
literally the Boolean maximization. Remove the diagonal first. The
diagonal-deleted quadratic form is multi-affine and its absolute maximum
over the cube is at a Boolean vertex. Adding the diagonal back costs
at most `||diag U||op/2=o(1)` in normalized energy. Haar concentration
gives `max_i |U_ii|=o(1)` in probability, while the deterministic sequence
has this property directly from (1). This supplies the feasible-vector
version of (5) needed in (2).

## 5. Consequence and archive reconciliation

The Haar cap calculation itself was already present in the historical
campaign. The WZF theorem and its generalized-diagonal hypotheses were
also audited historically for a DIFFERENT purpose: high-temperature
free energy and the failure of Frobenius-near-conference structure to
imply uniform delocalization. See
`finite_temperature_universality_relaxation_audit.md`.

The present consequence uses exact involutions, where the required
uniform hypotheses are automatic, and exact gauge equivariance to give
a fixed-rule performance ceiling inside the original flat sign class.
It does not need high-temperature assumptions or an Ising free-energy
formula. It rules out proving a universal lower bound of 1/2 solely by
optimizing over this particular collection of fixed local algorithms.
It leaves all nonlocal constructions, dimension-growing algorithms,
matrix-adapted initialization, and the actual Hadamard Boolean cap open.
