# Marked lower bound for fixed finite-species variance profiles

Date: 2026-09-07. Status: **Pending independent audit**. This is a proposed
extension of the audited marked-response lower theorem, not favorable
flatification and not equality of optimized profile values.

## 1. Precise proposed theorem

Fix q species with positive limiting proportions lambda_s, and a symmetric
nonnegative matrix v_st satisfying

```math
\sum_t\lambda_t v_{st}=1\qquad(1\le s\le q).
```

For each N partition the coordinates into species V_s with
`|V_s|/N -> lambda_s`. Let W_N be hollow symmetric, with entries

```math
(W_N)_{ij}=\sqrt{v_{st}}\,a_{ij},\qquad
a_{ij}\in\{-1,1\},\quad i\in V_s,\ j\in V_t,
```

where zero v_st means a zero block. One may also use block amplitudes
converging to the displayed fixed values, provided the rows have squared
norm `(1+o(1))N` uniformly. Define the energy half-width by

```math
\mathcal W(W)=\frac{\max_x H_W(x)-\min_x H_W(x)}2.
```

The proposed conclusion is

```math
\liminf_{N\to\infty}\frac{\mathcal W(W_N)}{N^{3/2}}
\ge c_*=0.4333221116640807534\ldots.
\tag{1}
```

It implies the same lower bound for Q(W_N). The constant is the existing
certified marked-response value; no numerical certificate is changed.
This excludes sub-c_* bounded-amplitude counterexamples from every fixed
block-constant row-regular profile, including arbitrary two-block profiles.
It does NOT show that the optimized profile infimum equals the full-sign
infimum, which may exceed c_*.

## 2. Fixed-spectrum setup

First assume `B=W_N/sqrt(N)` has a fixed operator bound L. Its entries are
O(N^(-1/2)), its row squared norms tend uniformly to one, and
`Q=B²` is positive semidefinite with diagonal tending uniformly to one.
For exact row regularity, its diagonal is exactly one after the harmless
choice of normalization. The original proof tolerates the uniform o(1)
row error: each fixed tree has finitely many leaf summations, and the
correlation normalization changes only o(1) at fixed L.

Condition every fixed diagram on the species assigned to each vertex.
There are finitely many such assignments, depending only on q and the
fixed diagram. On a species assignment each doubled edge contributes its
constant v_st/N, while each single parity edge contributes an actual
signed block entry sqrt(v_st) a_ij/sqrt(N).

If a single parity edge has two free endpoints, sum it first with all
other labels fixed. The remaining factors at those endpoints are bounded
unary functions, after conditioning their species. A restricted block of
W has bilinear cube norm at most beta(W), by zero extension and cube
rounding; beta(W)<=N||W||op=O_L(N^(3/2)). For v_st>0, the sign block has
norm at most beta(W)/sqrt(v_st); a zero block makes the diagram vanish.
Thus every such normalized parity contribution gains O(N^(-1/2)). All
constants may depend on the finitely many positive v_st.

If no parity edge survives in a leading local moment, the same graph-
theoretic copy-pairing argument forces a doubled tree. Sum leaves using
`sum_j B_ij²=1+o(1)`. It gives the same Wick covariance and rooted
automorphism normalization as in the original theorem. At two outputs,
the whole branch pairing is exactly `sum_a B_ia B_ja=Q_ij`; summing the
species recovers Q, not a substitute block-averaged covariance.

Therefore the finite old Gaussian frame, the inverse isometry U, and the
whole-pairing covariance main term remain unchanged.

## 3. Generic cuts and nonlinear covariance remain unchanged

The global/proper tree-cut estimates use only bounded ||B||op, uniformly
bounded row norms, maxentry O(N^(-1/2)), and tree combinatorics. These
are precisely the hypotheses checked in Section 4 of
`decisive_audit_bipartite_marked_response_lower_2026_09_07.md` and the
underlying nonlinear covariance sources. They do not delete squared
edge weights and do not need the complete-bipartite support.

Hence genuinely partial forest matchings have the same Frobenius O(1)
and nuclear o(N) bounds. Whole matches give the exact Q-Schur main term:

```math
\left\|\operatorname{Cov}(R)-
\sum_{k\ge3,\ k\ {m odd}}w_k Q^{\circ k}\right\|_*=o(N).
```

The same covariance Cauchy--Schwarz estimates transfer raw local Wick
errors and collision deletion. The only profile-sensitive remaining
step is the exceptional no-free-free-edge old/new contraction.

## 4. Exceptional contractions use subcorrelations, not false scalarization

For each species let P_s be its diagonal coordinate projection and put

```math
Q_s=B P_s B.
```

Every Q_s is PSD, `diag Q_s<=diag Q=1+o(1)`, and `||Q_s||op<=L²`.
This replaces the unjustified assertion that a weighted doubled diagram
always reduces to an entrywise power of Q.

In the old-input/nonlinear-source pairing, if a free-free parity edge
survives, Section 2 supplies the same discrepancy gain. Otherwise the
parity graph has p length-two paths a-u-j and possibly an edge a-j,
with p>=2. After fixing all species, a length-two path through species s
contributes `(Q_s)_aj`. Any surviving doubled-edge weights are fixed
constants. Summing the finite species assignments therefore gives a
finite linear combination of kernels of the forms

```math
N^{-1/2} P_r
(Q_{s_1}\circ\cdots\circ Q_{s_p})P_t,
\qquad
P_r\big[B\circ Q_{s_1}\circ\cdots\circ Q_{s_p}\big]P_t,
\tag{2}
```

with the diagonal and injectivity collision terms treated as in the
original proof. Coefficients in this combination depend only on the
fixed diagram and profile.

Schur multiplication by a PSD matrix with diagonal at most 1+o(1) is
bounded on operator and Frobenius norms. Since `||Q_s||F<=L² sqrt(N)`,
the first kernel in (2) has Frobenius O(1); the second has Frobenius O(1)
by `max|B_ij|=O(N^(-1/2))`. Thus the total squared old/return covariance
remains O(1), as required. It is not necessary to identify these small
error kernels with Q^(circ p).

For the root-hit full contraction, retain diag(Q J B). In its doubled
diagram the same combinatorial argument leaves a single parity edge
between free labels. Condition species and use Section 2. Its squared
Frobenius contribution is O(sqrt(N)), with collision terms O(1), hence
E||J||F²=o(N). Bounded root transports finish that exceptional case.

## 5. Unchanged scalar certificate

At exact unit row norm, Q is a correlation matrix. The existing odd-Schur
variance averaging theorem applies directly, without a mixing assumption
on the species transition matrix `(lambda_t v_st)`:

```math
\sum_i\sqrt{\left(B\left[\sum_k w_kQ^{\circ k}\right]B\right)_{ii}}
\ge N\sqrt{\sum_k w_k}.
```

Uniform o(1) row-normalization errors are removed first at fixed L by
normalizing Q's diagonal and the same finite-tree continuity argument.
All subsequent marked/unmarked tests, fixed-noise soft-sign limits, and
literal cube means F+-H sign(BF) have their original normalizations.
The original 21-anchor/two-Gaussian certificate therefore yields (1) in
the fixed-spectrum class. Both extremal endpoints are retained; this is
a width lower bound rather than a one-sided surrogate.

## 6. Spectral deletion with species proportions retained

On a sequence with Q(W_N)<=C N^(3/2), the simultaneous Grothendieck
diagonal majorant deletes at most epsilon N heavy coordinates and leaves
operator norm O_C(sqrt(N)/epsilon). Let lambda_min=min_s lambda_s.
Choose a common retention fraction `theta=1-O(epsilon/lambda_min)`.
In each species, further delete arbitrary coordinates until its retained
count is theta times its original count, rounded by O(1). This is possible
for sufficiently large N and retains the limiting proportions exactly.

Every remaining block still has its original constant amplitude and sign
entries. After normalization by the square root of the retained total
order N', its row squared norms are uniformly 1+o(1) and its operator norm
is bounded at fixed epsilon. Restriction cannot increase the width:
averaging omitted independent spins puts every retained quadratic value
inside the original interval. Apply the fixed-spectrum result, pay the
factor `(N'/N)^(3/2)`, and send epsilon to zero.

If normalized Q(W_N) is unbounded but width is the object being minimized,
first note that a hollow quadratic has mean zero, so Q<=2 width. Thus any
sequence potentially violating (1) has the bounded cap needed for deletion.

## 7. Exact audit obligation and nonclaims

The crucial new claim needing an independent diagram audit is that after
species conditioning the no-free-free-edge exceptions have precisely
the finite Schur-subcorrelation form (2), including any residual doubled
edges. The generic cut estimates and scalar numerical policy are reused
only at their stated hypotheses.

No extension to an arbitrary bounded entrywise variance profile follows
merely from graphon approximation: cut-norm errors may be too large for
the required Frobenius/nuclear exceptional-contraction bounds. No favorable
flatification, optimized profile equality, or convergence result is asserted.
