# Global coherent contraction with a deterministic energy-tail budget

Date: 2026-09-06. Status: exact finite-dimensional inequality and uniform
truncation reduction. This is NOT yet an improved cap theorem: the remaining
coherent variational expression must still be bounded using the seed.

## Setup and deterministic normalization

There are n vertices and independent random outgoing arrays
`(X_ij: j != i)` in R^q. Each array is exchangeable in its n-1 slots, and
almost surely `sum_j ||X_ij||² <= C n`. No typical-spectrum assumption is
made. Let R be a symmetric orthogonal map on R^q and t>0. Define

`K_R(x,y)=exp[-t(||x||²+||y||²)] cosh(2t x^T R y)`

and `Z_R=E product_(i<j) K_R(X_ij,X_ji)`. Each kernel lies in [0,1].
In the grouped sign weave of the companion vector-seed artifact, q=k²,
R=R_A, n=d, X=U/sqrt(ell), and C=k²: Parseval gives at each group
`sum_beta ||U^(alpha,beta)||²/ell <= k m = k² d`.
Independent uniform permutations of the d spectral COLUMN BLOCKS at each
group provide the required exchangeability (the self slot is marginalized).
This is compatible with the literal Hadamard weave and with every retained
spin configuration. Sharing the block permutation across the k rows retains
the joint seed source. This statement does not itself count all spins.

## Exact bounded-energy reduction

Choose a cutoff B>0, delta=C/B²<1/2. At vertex i mark precisely those slots
with ||X_ij||>B; call their set S_i. Its size is at most floor(Cn/B²).
For a pattern S=(S_i), retain only edges ij for which j notin S_i and
i notin S_j; denote the resulting graph G_S. Discarding the other factors
increases the integrand, because 0<=K_R<=1. The patterns partition the
sample space, so

`Z_R <= sum_S E[ product_i 1_{pattern S_i} product_(ij in G_S) K_R(X_ij,X_ji)]`.

There are at most `P_n=[sum_(r<=floor(delta n)) binom(n-1,r)]^n` patterns.
For fixed delta<1/2, `limsup n^-2 log P_n <= h(delta)`.
Patterns of zero probability can be ignored. This bound pays the pattern
entropy, not the energy of discarded slots; energy can concentrate there.

## Multiplicative finite Fock truncation

Put a=2t B² and choose L with

`sum_(r>L) a^(2r)/(2r)! <= epsilon`.

The feature vector is

`phi_L(x)=exp(-t||x||²) direct_sum_(r=0..L)
           sqrt((2t)^(2r)/(2r)!) x^(tensor 2r)`

in the symmetric tensor spaces. Its dimension is
`D=sum_(r=0..L) binom(q+2r-1,2r)`, independent of n.
Let R_L be the orthogonal involution induced by the even tensor powers of R.
For ||x||,||y||<=B,

`0 <= K_R(x,y) <= (1+epsilon) <phi_L(x), R_L phi_L(y)>`.

Indeed the latter inner product is the positive even partial sum of cosh
times its positive exponential prefactor; its partial sum is at least 1.
The tail is at most epsilon. This is a POINTWISE multiplicative comparison,
not an additive approximation of a signed contraction.

## Coherent-state global bound

For each pattern define the real tensor, of degree d_i=deg_G(i),

`T_i^S=E[1_{pattern S_i} tensor_(j:ij in G_S) phi_L(X_ij)]`.

It is symmetric in its surviving legs: permutations of those slots preserve
the exact pattern event and the exchangeable source law. The other slots,
including low slots deleted because of their opposite endpoint, are simply
integrated out. Independence across vertices makes the expected truncated
product exactly the contraction of these tensors along G_S with edge R_L.

For any symmetric degree-d tensor on C^D, normalized complex-sphere measure
gives the exact resolution

`T = binom(D+d-1,d) integral u^(tensor d) <u^(tensor d),T> dmu(u)`.

This follows by the unitary-invariant resolution of identity on Sym^d(C^D)
(or directly by integrating complex monomials). Explicitly, for multiindices
alpha,beta of total degree d, the sphere integral of u^alpha times the
conjugate of u^beta is zero unless alpha=beta, and then equals
`(D-1)! product_j alpha_j! / (D+d-1)!`. In the normalized symmetric monomial
basis, multiply by `d!/product alpha_j!` to get exactly the reciprocal of
`binom(D+d-1,d)` on every diagonal entry. Degree zero is included.
Taking the absolute value only AFTER inserting all resolutions gives

`|contract_G(T_i;R_L)| <= product_i binom(D+d_i-1,d_i)
   sup_(||u_i||=1) product_i |<u_i^(tensor d_i),T_i>|
                   product_(ij in G) |u_i^T R_L u_j|`.

Complex conjugations in the coefficient depend on inner-product convention;
the edge contraction is bilinear and the displayed absolute values are
unambiguous. Real-sphere resolution cannot replace complex-sphere resolution
without additional trace terms.

Consequently Z_R is bounded by `(1+epsilon)^(n(n-1)/2)` times the SUM of
the displayed global coherent bounds over patterns S. Alternatively replace
the sum by P_n times its largest term. For fixed q,t,B,L, the frame factor
has logarithm at most `n log binom(D+n-2,n-1)=O_D(n log n)=o(n²)`.
The total leading overhead in the max-pattern formulation is at most

`h(C/B²) + (1/2)log(1+epsilon) + o(1)`

after division by n². First take n to infinity for fixed B,L, then B to
infinity and epsilon to zero, increasing L as needed. Thus a target-order
independent feature dimension SUFFICES with arbitrarily small leading loss.

## What remains, and what this avoids

This bound retains R on every surviving edge and the entire joint vertex
polynomial. It does not replace edges by local Hilbert norms. It handles
ALL arrays satisfying the deterministic energy bound, including concentrated
spectra. However it supplies no seed-dependent upper bound for the final
supremum and no reason that the bound is sharp. Its maximizing coherent
vectors may erase the seed through another mechanism. Also the max over
patterns is a real adversarial relaxation: a useful seed estimate must be
robust to deletion of at most delta n² edge incidences. No convergence of
actual signing optima or equality of ensemble pressure is claimed.
