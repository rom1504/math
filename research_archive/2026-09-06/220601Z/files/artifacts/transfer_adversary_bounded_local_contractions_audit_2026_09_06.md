# Independent audit: bounded local full contractions and weighted tails

2026-09-06. **PASS at the local, fixed-degree scope** of
`transfer_seed_bounded_local_full_contractions_2026_09_06.md`, including
its trigonometric weighted-tail extension. The new tail argument answers
a specific limitation raised during the audit; it is not a two-root
Frobenius estimate.

## 1. Exact bounded chain and degree isolation

For a fixed-degree homogeneous noise X of degree P and coherent W of
degree at most M<P, finite differences give

```
Delta_a A(W)=sum_l A_l'(W) Delta_a W_l+R_a,
|R_a|<=C ||Delta_a W||^2.
```

Thus the first terms in Gamma(X,A(W)) are bounded row coefficients
times the fixed-degree polynomial Gamma(X,W_l). Its proper-contraction
L2 bound is O(epsilon_n), and fixed-degree hypercontractivity preserves
that small factor in every needed fixed higher moment. For the error,

```
|sum_a Delta_a X R_a|
 <= C max_a|Delta_a X| sum_a ||Delta_a W||^2.
```

The first factor has L4 norm `O(epsilon_n polylog(n))`, by the small
noise influences and fixed-degree hypercontractivity at a logarithmic
moment order. The second has bounded L4 norm by Minkowski and the
finite total derivative energy of W. Holder proves the claimed small
L2 norm of the exact Gamma. High individual influences in W cause no
problem for this LOCAL statement.

At output degree d, the product formula has only right degrees
`e=d-P+2r`, for 1<=r<=P. Thus all but the r=P term are bounded by a
proper X cut; their right norms are at most ||A(W)||2. Orthogonal
projection of Gamma onto d isolates the full contraction with A_(d+P).
Infinite right degrees cannot enter this fixed output degree.

The possible symmetrization concern is resolved exactly: a FULL
contraction leaves no unmatched X slots. All remaining slots belong
to the symmetric A_q tensor, so the resulting tensor is already
symmetric. There is no discarded antisymmetric component. Proper
contractions may be symmetrized and projected, but their norms were
already small before those contractions.

## 2. Why the trigonometric Sobolev estimates are uniform

For every fixed derivative order s, the exact cover expansion of
exp(itW) expresses each derivative as a finite sum of products of
factors bounded by a constant times |Delta_V W|. The map from a full
U tuple into its covered subtuples is injective. For a fixed root, the
sum of squares over U is therefore at most the product of the factor
row-Hilbert sums, before Holder averages over the seeds.

Each derivative row-Hilbert norm of W has every fixed Lp moment bounded
by fixed-degree hypercontractivity and

```
sum_(|V|=r) E|Delta_V W|^2
 =sum_q binom(q,r)||W_q||2^2 <= binom(M,r) Var(W).
```

This proves a bound independent of n for
`sum_q binom(q,s)||A_q||2^2`. Constants may grow with s and frequency,
which is harmless because these are fixed before the order limit.

## 3. The multiplication bound is polynomial in the large degree

For homogeneous X_P and Y_q, the contribution with r shared labels
has normalized coefficient factor

```
c_(P,q,r)=r! binom(P,r)binom(q,r)
          sqrt((P+q-2r)!/(P! q!)).
```

At fixed P this is `O_P((q+1)^(P/2))`, uniformly over r<=min(P,q).
For large q, the binomial factors contribute q^r and the square-root
factor contributes q^(P/2-r); the finitely many small q are absorbed
in the constant. The underlying Hilbert contraction has norm at most
||X||2||Y_q||2. Boolean distinctness is a projection, not a norm loss
in the opposite direction.

Each output degree receives at most P+1 choices of q. Cauchy--Schwarz
over that fixed number yields

```
||XY||2^2 <= C_P ||X||2^2 sum_q (q+1)^P ||Y_q||2^2.
```

There is no hypercontractive factor exponential in q. Combining this
inequality with the arbitrarily high fixed Sobolev tail estimates proves
uniform weighted L2 convergence of X times the original-degree tails
of a fixed trig response. Finite tails are Cauchy in L2 and converge
in probability to the intended product, identifying the limit.

This justifies ordered LOCAL product surgery: fix a large original
degree cutoff, take n to infinity using the fixed-degree contraction
bounds, then remove the weighted tail. It does not upgrade one local
O(n^-1/2) factor into the two-root small norm required after unrestricted
root summations. That additional matrix/tensor issue remains separate.

## 4. Rate-free averaged version and the exact old-rich dependency

The source's subsequently added averaged estimate also passes. Put
`delta_i=max_a Inf_a(X_i)` and
`gamma_i^2=sum_l ||Gamma(X_i,W_li)||2^2`, with fixed degrees and bounded
row variances. Then

```
E max_a|Delta_a X_i|^4
 <=sum_a E|Delta_a X_i|^4
 <=C sum_a Inf_a(X_i)^2 <=C delta_i.
```

The old derivative energy has bounded L4, hence
`||Gamma(X_i,A(W_i))||2<=C(gamma_i+delta_i^(1/4))`.
Squaring, averaging and Jensen give
`C(average gamma_i^2+sqrt(average delta_i))`. No unproved logarithmic
rate is needed in this averaged version.

Read Sections 10--12 of
`resumed_bound_audit_full_nonlinear_covariance_trace_2026_09_06.md`
to check the proposed source of gamma control for exceptional rich
primitive returns. Section 11 proves deterministic proper/full kernel
contraction estimates for the fully injective nonlinear forest versus
the exact old marked-tree return. These imply averaged Boolean Gamma
L2 bounds directly: the same diagonal-free tensors occur for Boolean
and Gaussian inputs, and the Boolean product formula imposes only the
corresponding symmetrizations and distinct remaining-label projections.
Fixed-degree main/raw errors pass under the stated row-moment controls.

The correct bridge is those deterministic contraction norms. Section
12's sign-replacement estimate with the exceptional return entering
only LINEARLY would NOT, by itself, justify the quadratic Gamma norm.
This extraction retains the old theorem's fixed-L, exact injective
forest and marked-tree hypotheses. It does not generalize a Gaussian
local law merely by analogy.

## 5. Remaining rich returned-field assembly caution

The alternative mixed-covariance theorem retains the literal rest
R=BF-Z_P. This is safe for the dual split, which separately checks its
original-degree components. It is NOT automatically safe to treat the
enlarged list (W,R) as the primitive list in the high-degree trig
transport theorem: a finite source approximation to R can have degree
D>kP, and an output degree q>kP need not exceed D.

The returned-field assembly must either use the full-noise Hermite
expansion relative to the ORIGINAL primitive list of degree at most
M<P, justified by the old full local separation, or explicitly peel
off the one-factor high-degree R terms into the deterministic linear
returns controlled by Hall. The mean coefficient a is deterministic
and is not the issue. This is an actual hypothesis-matching condition,
not an assertion that the remaining rich theorem is false.

The seed's frozen rich-extension checkpoint keeps the full returned-
field conclusion open. The separately proved local and trig modules
must not be presented as though this final assembly were already a
theorem.
