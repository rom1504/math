# Coherent polynomial Walsh inputs: summable covariance and transported cuts

Date: 2026-09-06. Status: structural proofs reconstructed independently
with the audit agent; submitted as a bounded-degree module. This does
not assert closure under arbitrary feedback iteration.

## 1. Setup

Let `S` have `N` independent sign coordinates, where `N/n` is bounded
(a fixed number of seed colors is allowed). Let `M_1,...,M_k` be
`n`-by-`N` matrices with operator norms at most a fixed `L_0`. Fix a
polynomial `P` of degree at most `D`, independently of `n`, and put

```math
d_i=P((M_1S)_i,\ldots,(M_kS)_i).
```

Let `d_q` denote its EXACT Boolean Walsh component of degree `q`.
The constants below may depend on `P,k,D,L_0`, but not on matrix order.
The same assertions hold for cross-covariances of two such inputs.

Let `B` be a symmetric hollow normalized signing on the `n` output
coordinates with a fixed operator cap. Its exact flat-entry property
is used below. The degree-`q` transported field is `Y_q=B d_q`.

## 2. Exact row-tensor form

Partition the original input positions in each monomial by their equal
seed labels. Odd blocks remain marked after `S_a²=1`; even blocks are
summed out. Inclusion-exclusion separates the disappeared blocks from
the remaining marked labels. At each fixed Walsh degree, this yields
a finite linear combination of kernels of the form

```math
D_c\,\Pi_{\rm distinct}
       [\text{row tensor of }N_1,\ldots,N_q],                 (1)
```

where `D_c` is a bounded row-scalar diagonal, and each factor matrix
`N_l` is a nonempty entrywise product of some `M_j`. The marked-slot
distinctness projection is retained EXACTLY.

Every `N_l` has bounded operator norm. A one-factor product inherits
its bound. A product with at least two factors has bounded absolute
row AND column sums, by retaining two factors in Euclidean
Cauchy--Schwarz and bounding the remaining entries. The disappeared
even blocks have at least two factors, so their row scalars are bounded
by the same calculation.

The global root map of (1) has bounded operator norm: before projection,
its Gram is a Schur product of bounded-op Gram matrices with bounded
diagonal. Orthogonal projection on the common marked tensor space does
not increase this norm. Thus every positive Walsh degree has a bounded
global root map. Constant Walsh degree zero is excluded from this
statement; its uncentered covariance can have operator norm of order `n`.

## 3. Higher Walsh covariance has bounded absolute row sums

For every fixed `q>=2`,

```math
\max_i\sum_j|\mathbb E[d_{q,i}d_{q,j}]|=O(1),\qquad
\max_j\sum_i|\mathbb E[d_{q,i}d_{q,j}]|=O(1).                 (2)
```

To prove this, expand the two kernels in the finite form (1), pair
their marked slots, and apply inclusion-exclusion to the remaining
distinct-label constraint. Each resulting term is a Schur product of
Gram factors, one for each remaining free label block.

If there are at least two free blocks, retain two Gram factors. Their
absolute row sum is bounded by Cauchy--Schwarz and their bounded
operator norms; the same holds for column sums. Further Gram factors
have bounded entries.

If there is exactly one block, the term has the form `U V^T`, where
both `U` and `V` are entrywise products of the `q>=2` row-factor
matrices on their respective sides. Both `U,V` have bounded absolute
row and column sums. Therefore `U V^T` has the same property. Bounded
row scalars and the finite sum preserve these estimates, proving (2).

This proof also handles cross-covariance of two different degree-`q`
coherent inputs. Different Walsh degrees are orthogonal exactly.
In particular, for `q>=2`,

```math
\|\operatorname{Cov}(d_q,\widetilde d_q)B\|_{\max}
=O(n^{-1/2}).                                                (3)
```

## 4. Only the first Boolean Walsh degree contributes coherent energy

Let `K_ia=E[S_a d_i]`. For a globally odd coherent polynomial, the
constant component is zero. From (2), flatness of `B`, and exact
orthogonality of Walsh degrees,

```math
\frac{\mathbb E[d^{\mathsf T}Bd]}{2n}
=\frac{\operatorname{Tr}(B K K^{\mathsf T})}{2n}+O(n^{-1/2}).  (4)
```

Indeed every fixed degree `q>=2` contributes at most its absolute
row-sum bound times `1/sqrt(n-1)` after normalization.
This is a BOOLEAN first-chaos statement, not a Gaussian replacement of
the coherent return. In particular `K` can retain strong input atoms.

A bounded coherent response can inherit (4) whenever it is approximated
in averaged `L²` by fixed coherent polynomials in the stated ordered
limit. First Walsh projection is an `L²` contraction, so the right side
is continuous under that approximation as well. No uniform operator
bound on the first-Walsh kernel of an arbitrary bounded response is
asserted solely by this continuity argument.

## 5. Proper transported flattenings are small

For every fixed `q>=2`, each proper fixed-root flattening of the exact
degree-`q` kernel of `Y_q=B d_q` has operator norm `O(n^{-1/2})`.

Use inclusion-exclusion to expand the distinctness projection in (1),
and fix a nontrivial left/right cut of the `q` formal slots. Consider
one equality partition term.

If no equality block straddles the cut, collect the row factors in
each side. Their global root maps have bounded operator norm. The
transported flattening is their product with middle diagonal
`diag(b_i)`, whose norm is `1/sqrt(n-1)`.

If a block straddles the cut, its common label appears on both sides,
so the flattening is block-diagonal in that label (and in any other
straddling labels). Fix these shared labels. The row weight of the
root sum contains the product of at least TWO original row factors
at one shared label. Its absolute COLUMN sum over the root index is
bounded by Cauchy--Schwarz. Multiply by the flat `B` root coefficient;
the sum of the absolute weights is `O(n^{-1/2})`. Remaining shared
factors are entrywise bounded, and the residual left/right row tensors
have bounded Hilbert norms. Each diagonal block therefore has operator
norm `O(n^{-1/2})`, uniformly in its fixed shared labels. Taking the
supremum over the blocks proves the same bound for the flattening.

There are only finitely many partition terms at fixed degree. Their
sum proves the claim. The global root map of `Y_q` is also bounded,
by multiplying that of `d_q` by the bounded operator `B`.

## 6. Relevance and exact limitation

The odd-degree channels `Y_q`, `q>=3`, satisfy two key estimates
used by the retained-coherent energy argument: proper transported cuts
are small, and a complete higher-degree input covariance followed by
`B` is entrywise small. They therefore provide candidate fixed
polynomial channels beyond scalar or independently colored functions
of `BS`. Extending the full energy theorem to them still requires the
mixed-star bookkeeping with their EXACT source Walsh projections.

In particular, one must not restore an unrestricted pre-Walsh source
tensor at small Hilbert cost. On an actual involution comparator `Q=I`,
the input `(QS)_i³=S_i` has zero degree-three Walsh part, whereas the
unrestricted Gaussian transported kernel `sum_a B_ia e_a^{tensor3}`
has squared row Hilbert norm one, entirely on repeated slots. The
small proper flattening in Section 5 does not make that source
collision error small. Its exact partition proof keeps the Walsh
projection and avoids this error.

For example, the actual own-spin polynomial
`S_i h_2((BS)_i)` uses only the bounded-op matrices `I,B`. Its exact
positive input degree is three because `B_ii=0`; after transport it is
the first nontrivial canonical old-tree creation field, up to the
already controlled output-root collision.

The lemma DOES NOT say that an arbitrary bounded-op mixture of such
polynomial channels, followed by an arbitrary nonlinear coordinate
response, remains a polynomial of bounded-op LINEAR seed fields.
In particular, `Q[S K(X)]` is not automatically in the coherent class
defined in Section 1. Nor does bounded covariance alone give the
row-tensor structure used in Sections 3 and 5. A complete new feedback
iteration must check those structures for its new coherent inputs,
all input degrees, and all mixed star contractions separately.
