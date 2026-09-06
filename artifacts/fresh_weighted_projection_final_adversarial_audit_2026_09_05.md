# Final adversarial audit of the weighted original lower-bound mechanism

Date: 2026-09-05. This is a fresh reconstruction of the delicate steps,
not an assertion that earlier agreement substitutes for proof. The
cubic theorem and every fixed odd-degree extension pass this audit.

## 1. Exact signing hypotheses remain in use

Throughout the energy proof, `B=A/sqrt(n-1)` is symmetric and hollow,
has every off-diagonal entry of modulus `1/sqrt(n-1)`, and has a fixed
operator bound `L`. In particular its rows are unit vectors. The
regularized principal submatrix has these same properties after its
own order is used in the normalization.

No step extends the own-coordinate argument to nonhollow matrices.
No step replaces a deterministic signing by an arbitrary random-matrix
model, and no step assumes a Gaussian law for an arbitrary sign-dependent
field. The cosquare full-sign construction is a separate structural
example, not an input to this lower-bound proof.

## 2. The global derivative bound imports exactly the needed theorem

I reread the actual primary statement of Theorems 2.1 and 2.4 in
[Bandeira--Lucca--Nizic-Nikolac--van Handel, *Matrix Chaos Inequalities and
Chaos of Combinatorial Type*](https://web.math.princeton.edu/~rvan/chaosconf241224.pdf).
They give square-free decoupling followed by an expected operator bound
in terms of complementary coefficient flattenings that keep the two
matrix indices on opposite sides. This is precisely the flattening
family controlled in the tree derivative module.

The original global tensor has all marked vertices and the external root
as indices. Crossing-edge label-copy maps are isometries. The internal
forest count `I<=E_internal` cancels every possible `sqrt(n)` factor from
a crossing-isolated vertex. Collision deletion has global Hilbert norm
`O(1)`, sufficient for global flattenings. When the external root is
fixed, the corresponding collision norm is `O(n^-1/2)`.

Differentiation contracts fixed slots against the row tensors `b_j^r`.
Their Gram matrix is `Q^(circ r)`, whose norm is at most `L^2`. Putting
all contracted slots on the side containing the new derivative index
`j` is a legitimate composition of an original global flattening. The
remaining stochastic indices are still distinct, so the primary
square-free hypothesis is satisfied literally. Fixed higher moments
follow from Banach-valued hypercontractivity; they are not inferred from
an expectation estimate alone.

## 3. The one-level replacement preserves exact own-coordinate exclusion

The field used in the difficult integration by parts is

`Y_T,i=sum_l B_il N_l h_T(X_l)`,

with the child fields `X_l` still the actual injective fields. This is
not a fully recursive raw-product substitution. Hermite products cancel
full child pairings; surviving proper contractions have norm
`O(n^-1/2)`. The Gaussian own-root creation inequality transports the
remaining error at the same order because every child error excludes
`N_l` exactly. Cross-child and external-root collisions also have
fixed-root Hilbert norm `O(n^-1/2)`.

Thus `Y-X=O_Lp(n^-1/2)` uniformly at each root. Fixed unit-direction and
ordinary-coordinate derivatives retain that order by the finite-degree
Gaussian derivative inequality. Bounded-operator energy continuity
allows the one-level replacement before the chain calculation.

For `K_r(l,j)=D_(b_j)^r h_T(X_l)`, the entire row `l` excludes `N_l`.
The exact matrix identity is

`D^r Y_T=B diag(N) K_r+rB(B circ K_(r-1))`.

The two root-not-hit traces, after rearranging all indices, are

`sum_l E N_l [B diag(F_T) B diag(M) B K_r^T]_ll`,

and

`sum_l E N_l [B diag(M_T(BF)) B K_r^T]_ll`.

In each diagonal entry only row `l` of `K_r` occurs. Therefore its
own-coordinate derivative is exactly zero. This is the essential
cancellation; replacing it by a purported small operator norm of the
whole third-derivative matrix would be incorrect.

The all-on-M case uses the Frobenius coordinate Jacobian of
`M_T(Y_i)(BF)_i`. Its bound sums `||(BF)_i||_2^2` over roots and uses
fixed-degree hypercontractivity. It does not require a false uniform
pointwise bound on every transported response `(BF)_i`.

## 4. Every derivative split is accounted for

When both response factors receive derivatives, the entire term is

`n^-1 <B circ D^a M, B D^b F>_F`, `a,b>=1`.

Flatness makes the first Frobenius norm polylogarithmic; the second is
at most `sqrt(n)` times a polylogarithm. These terms vanish.

Within one response, a chain partition with at least two nonedge factors
has pointwise size `O(n^-1)`. With exactly one nonedge factor and another
factor, the latter must be an edge first derivative and supplies a power
of `Q_ij`; its square sum is `O(n)`. The sole nonedge high derivative is
the own-root trace handled above. The all-edge term on the local weight
vanishes using `sum_j |Q_ij|^r<=sum_j Q_ij^2` for odd `r>=3`.

The only surviving term is the all-edge derivative on `F`. The two-root
moment theorem factors its local weight with error `O(|Q_ik|)+o(1)`.
Its signed weights are `d_i B_ik(B Q^(circ r))_ik`; their absolute sum
is `O(n)` and the extra `|Q_ik|` sum is `O(sqrt(n))`. Their exact signed
row sums are `d_i(B Q^(circ r)B)_ii`. No cancellation between separately
maximized sign energies enters this calculation.

## 5. Sign input is compared at the complete functional level

For a hybrid Gaussian/Rademacher input and any fixed polynomial `f`,
product Poincare/Efron--Stein gives

`Var(sum_j w_j f(G_j))<=C_(f,L)||w||_2^2`.

The continuous gradient is `B[w f'(G)]`. Each Rademacher central
difference differs from that derivative by `O_Lp(||w||_1/(n-1))`;
summing its squared errors is bounded by a constant times `||w||_2^2`.
This estimate is valid in the actual hybrid law. Even Hermite means are
`O(1/n)` by scalar fourth-order replacement, and odd means are zero.

These facts prove all ordinary input-derivative scales for
`Z_r=B h_r(G)`. The only exceptional scale through order four is

`partial_k^3 Z_(3,i)=sqrt(6)Q_ik/(n-1)`.

It has the required `O(n^-1)` vector norm but not a uniform pointwise
`O(n^-3/2)` bound. The complete energy calculation uses the vector norm.
Freezing a replacement coordinate is justified by first making it
Gaussian and projecting onto its finite one-coordinate Hermite basis;
one cannot recover those coefficients from only two sign endpoints.

The complete normalized weighted energy has fourth input derivative
`O_L1(n^-2)`. Moment matching through degree three then costs `O(n^-1)`
after all replacements. This is the energy bridge; rootwise convergence
alone is not substituted for it.

## 6. The local Gaussian limit used for cutoff errors is separately proved

With Gaussian input the higher odd transport is
`I_r(K_(r,i))/sqrt(r!)`, where
`K_(r,i)=sum_j B_ij b_j^(tensor r)`.
Every proper flattening is bounded by `L^2/sqrt(n-1)` and its Hilbert
norm is bounded. Hence its nontrivial contractions vanish. Equal-degree
old-tree covariances are controlled by the explicit child contraction
recursion; old-tree collision errors are paired with the already summed
kernel, not summed column by column.

Repeated unmarked slots have Hilbert norm `O(n^-1/2)`, using their exact
`1/(n-1)` contraction identity. The sign-row Hermite/multilinear error
is `O_Lp(1/n)` for fixed degree, so bounded-operator transport gives a
uniform root error `O(n^-1/2)`. The remaining fixed-degree multilinear
kernels have vanishing influences. Their invariance and Gaussian chaos
limit give the asserted local joint independence from the old family.
Roots with vanishing limiting variance are covered by their vanishing
second moment, not by dividing by that variance.

## 7. Final verdict and the remaining original convergence obligation

The cubic weighted theorem, its arbitrary fixed odd-degree extension,
the normalized local gain, and the uniform scalar-family escape survive
this adversarial reconstruction. Degrees and smooth approximations remain
fixed until the matrix-size limit. The later finite-cover argument fixes
a finite degree set before choosing a near-optimal scalar target.

The strongest precise remaining sufficient obligation toward convergence
is the already stated **all-order one-profile upper recovery**: for an
action limit of a bounded-operator, spectrally regularized minimizing
subsequence, construct hollow symmetric flat sign matrices at every
sufficiently large order whose one-profiles lie in vanishing neighborhoods
of the target profile. The quadratic energy is continuous in that signed-
action topology under the fixed operator cap, so this would convert the
liminf subsequence into an all-order limsup upper bound.

No result in this audit supplies that recovery construction. The exact
cosquare scalable gap now proves that retaining only the squared
operator cannot substitute for the signed-action data in that obligation.
