# Adversarial verification of the known-model validation draft

**Scope.** This is an independent audit of
`phase2_known_model_validation.md`, with emphasis on KM.2, KM.3, and the
claimed dimensional obstruction.  It reconstructs the proofs rather than
using the research narrative.  It makes no surface-file edits and does not
assert literature novelty.

## 1. Verdict

The formal cores of KM.2 and KM.3 are correct for arbitrary norms.  The
Hausdorff direction, the Lipschitz consequences for both suprema and infima,
and the `A/B` obstruction all check.  KM.3 admits an immediate strengthening:
the number of exceptional summands is controlled by the **effective affine
difference rank**

```math
r=\dim\operatorname{span}\bigcup_i(E_i-E_i),
```

not necessarily the displayed ambient dimension `p`.  The draft's bound by
`p` is valid but is not the intrinsic version.  Its stopping condition should
therefore be described as growing effective response rank, not merely growing
ambient dimension.

There is one important qualification to the conceptual claim.  A convex
response body is demonstrably a strict sub-landscape state only when that
body itself has a sub-landscape description (as a generator list for a
zonotope, a fixed-dimensional polytope of controlled complexity, or an
oracle declared as the state).  KM.3 by itself does not bound the description
complexity of arbitrary `conv(E_i)` or of their Minkowski sum.  It is a real
approximation theorem, but not a general compression theorem without this
additional representation hypothesis.

The Gaussian response-law minimality and the homogeneous max-plus cyclicity
claims are correct in their stated scopes, subject to the scope clarifications
in Sections 7--8 below.

| Claim | Audit result |
|---|---|
| KM.2(a), exact support roof | correct for the full labelled linear-field query family |
| KM.2(b), target-distance state | correct; the zero set recovers the signed-sum set |
| KM.2(c), floating-variable bound | correct, including the sharper `sum of r largest norms` constant |
| KM.3, Shapley--Folkman bound | correct; strengthen `p` to effective rank `r` |
| Lipschitz sup/inf transfer | correct in both directions because `E subset K` and `K` is within the stated distance of `E` |
| Lifted-energy interpretation | correct for additive components and uniformly Lipschitz aggregate queries; it needs bounded component diameters and includes the energy coordinate in the effective rank |
| `A/B` obstruction | exact; the gap is `2d`, a positive fraction of either natural extensive scale |
| Gaussian full-field minimality | correct for labelled Gaussian ensemble laws, not quenched samples modulo self-averaging |
| Max-plus cyclicity | correct for finite-entry (hence irreducible) homogeneous `K`; transient and an admissible eventual period are part of the data |

## 2. Effective-rank Shapley--Folkman lemma

The version actually used can be stated intrinsically.

### Lemma V.1

Let `E_1,...,E_n` be nonempty compact subsets of a finite-dimensional normed
space.  Fix `e_i^0 in E_i` and put

```math
L=\operatorname{span}\bigcup_{i=1}^n(E_i-E_i),
\qquad r=\dim L.
```

For every

```math
z\in\sum_i\operatorname{conv}E_i
```

there is a set `J subset {1,...,n}` with `|J|<=min(r,n)` and a representation

```math
z=\sum_{i\notin J}e_i+\sum_{i\in J}y_i,
\qquad e_i\in E_i,\quad y_i\in\operatorname{conv}E_i.       \tag{V.1}
```

#### Reconstruction

Translate each set by `e_i^0`.  All translated variation lies in the common
`r`-dimensional space `L`.  Applying the standard Shapley--Folkman lemma in
`L` yields (V.1) with at most `r` exceptional sets.  Equivalently, expand
each `y_i` as a finite convex combination and choose an extreme feasible
coefficient vector.  There is at least one positive coefficient per group;
the `n` normalization equations and `r` vector equations allow at most `r`
additional positive coefficients, so at most `r` groups can remain
nonintegral.  This also shows why ambient dummy coordinates must not count.

The original KM.3 statement uses `p>=r`, so it remains true as written.
Replacing `p` by `r` makes its mathematical boundary precise.

## 3. Hausdorff constant for arbitrary compact nonconvex sets

Let

```math
E=\sum_iE_i,
\qquad K=\sum_i\operatorname{conv}E_i=\operatorname{conv}E,
```

and let `Delta_i=diam(E_i)` in an arbitrary fixed norm, sorted decreasingly.
Take `z in K` and a representation (V.1).  For every exceptional `i`, select
any `e_i' in E_i`.  If

```math
y_i=\sum_j\alpha_je_{ij},
```

then

```math
\|y_i-e_i'\|
\le\sum_j\alpha_j\|e_{ij}-e_i'\|
\le\Delta_i.                                      \tag{V.2}
```

Replacing every exceptional `y_i` by `e_i'` produces `e in E` with

```math
\|z-e\|
\le\sum_{i\in J}\Delta_i
\le\sum_{i=1}^{\min(r,n)}\Delta_i.                \tag{V.3}
```

There is no hidden Euclidean step: (V.2)--(V.3) use only convexity of the
norm and the triangle inequality.

The Hausdorff direction is also correct.  Since `E subset K`,

```math
\sup_{e\in E}\operatorname{dist}(e,K)=0.
```

Equation (V.3) controls the other directed distance `K -> E`, hence

```math
d_H(E,K)
\le\sum_{i=1}^{\min(r,n)}\Delta_i.                \tag{V.4}
```

Compactness guarantees that the extrema and nearest points invoked in the
draft exist.  Boundedness plus closure could replace compactness in suitable
finite-dimensional variants, but the stated hypothesis is clean and enough.

## 4. Lipschitz suprema and infima

Put `delta=d_H(E,K)` and let `Psi` be `L`-Lipschitz.  Because `E subset K`,

```math
\sup_E\Psi\le\sup_K\Psi,
\qquad
\inf_K\Psi\le\inf_E\Psi.                         \tag{V.5}
```

For each `z in K`, choose `e in E` with `||z-e||<=delta` (or use an
arbitrarily small slack).  Then

```math
\Psi(z)\le\Psi(e)+L\delta,
\qquad
\Psi(z)\ge\Psi(e)-L\delta.                        \tag{V.6}
```

Taking the appropriate extrema and combining with (V.5) gives

```math
|\sup_E\Psi-\sup_K\Psi|\le L\delta,
\qquad
|\inf_E\Psi-\inf_K\Psi|\le L\delta.             \tag{V.7}
```

Thus KM.17 is correct for both signs.  The conclusion is useful only when
the query's Lipschitz constant is uniformly controlled on the chosen scale;
if `L` itself grows like the number of factors, the nonaccumulating geometric
defect need not yield a subleading response error.

## 5. Reconstruction of KM.2 and its sharper constant

For vectors `v_i` define

```math
P_z=\{t\in[-1,1]^n:\ \sum_it_iv_i=z\}.
```

This is a nonempty compact polytope for every `z in Z(V)`.  At an extreme
point of `P_z`, at most

```math
r=\dim\operatorname{span}\{v_i\}
```

coordinates can lie strictly inside `[-1,1]`: if more did, the corresponding
vectors would have a nonzero linear dependence, and a sufficiently small
positive or negative perturbation along it would remain in the cube while
fixing `z`, contradicting extremality.

Round only those fractional coordinates to their nearest signs.  For each
one, `|t_i-epsilon_i|<=1`, so for every norm

```math
\left\|z-\sum_i\epsilon_iv_i\right\|
\le\sum_{i\text{ fractional}}\|v_i\|
\le\sum_{i=1}^{\min(r,n)}a_i.                     \tag{V.8}
```

This verifies KM.14.  It is sharper than applying KM.3 blindly to
`E_i={-v_i,+v_i}`, whose diameter is `2||v_i||`; the nearest-sign rounding
uses the radius of the segment from an interior coefficient, not its full
diameter.

Since `S(V) subset Z(V)`, for every target `t`,

```math
\operatorname{dist}(t,Z(V))\le D_V(t).
```

Choose a closest `z in Z(V)` and then a signed sum within the right side of
(V.8).  The triangle inequality yields KM.15.  No strict convexity or inner
product structure is used.

The other two parts are also exact:

* with base energy zero, every feasible barycentric fibre has expected
  energy zero, so the upper roof is the constant zero function on `Z(V)`;
  its full linear support response is `h_Z(theta)=sum_i|<theta,v_i>|`;
* the zero set of the complete target-distance function is exactly `S(V)`,
  so that query family recovers the holes that the zonotope forgets.

The word "minimal" in the support statement is valid for the full labelled
linear-field experiment, up to one-to-one recoding of the support function
or closed convex body.  It should not be extended to a restricted set of
directions without checking that those directions determine the body.

## 6. Lifted energy sets: what KM.3 does and does not give

For additive component landscapes `H_i:X_i -> R` with features
`phi_i:X_i -> R^d`, set

```math
E_i=\{(\phi_i(x),H_i(x)):x\in X_i\}\subset\mathbb R^{d+1}.  \tag{V.9}
```

Then `sum_i E_i` is exactly the reachable table of aggregate feature and
total energy over the product state space.  Its convexification is the
Minkowski sum of the component lifted bodies, so composition is exact at the
convex-body level.  For a query

```math
\Psi(u,h)=h+g(u),                                  \tag{V.10}
```

KM.3 compares the true discrete extremum with the convex relaxation whenever
`Psi` is uniformly Lipschitz in the chosen product norm.  For fixed `u`, the
coefficient of `h` in (V.10) is positive, so optimizing over the convex body
uses only its upper boundary.  This justifies the response-roof
interpretation.

Three limitations should be explicit.

1. The effective rank is the dimension of the span of all differences in
   (V.9), generally as large as `d+1`; the energy coordinate does not
   disappear merely because only the upper boundary is queried.
2. The component diameters must be bounded on the normalization of interest,
   and the Lipschitz constant of `g` must not erase the gain.
3. For arbitrary compact `E_i`, storing their exact convex hulls or total
   Minkowski sum can itself have high descriptive complexity.  A strict
   sub-landscape state follows in the zonotope/fixed-complexity cases, not
   from Shapley--Folkman alone.

Within those limits this is a genuine, non-tautological composition result:
the discrete holes create a total error charged to at most `r` components,
not one error per component.

## 7. Exact check of the dimensional obstruction

For each coordinate, the `A` pair contributes

```math
2\epsilon_1e_j+2\epsilon_2e_j\in\{-4e_j,0,4e_j\},
```

whereas the `B` pair contributes

```math
3\epsilon_1e_j+\epsilon_2e_j
\in\{-4e_j,-2e_j,2e_j,4e_j\}.
```

Their convex hulls are both the coordinate interval `[-4e_j,4e_j]`.
Minkowski addition across coordinates therefore gives exactly

```math
Z(A)=Z(B)=[-4,4]^d.
```

The `A` set contains zero by cancelling each equal pair.  Every coordinate
of every `B` signed sum has absolute value at least two, and choices are
coordinatewise independent.  Hence

```math
\min_{s\in S(A)}\|s\|_1=0,
\qquad
\min_{s\in S(B)}\|s\|_1=2d.                       \tag{V.11}
```

This is a macroscopic gap.  There are `n=2d` component vectors, the maximum
component norm is `R=3`, and the gap is `(1/3)nR`.  Alternatively, both lists
have total generator mass `4d`, and the gap is half that scale.  Thus no
dimension-free error `o(nR)` can follow from the zonotope alone.  More
intrinsically, the effective variation rank here is `r=d=Theta(n)`; adding
unused ambient coordinates would not create an obstruction.

The example proves exactly the claimed information loss: every linear
support response agrees, while the nonlinear target query differs.  It does
not prove that every high-rank structured family has a large gap; it only
rules out a uniform dimension-free theorem for arbitrary lists.

## 8. Gaussian response-law minimality

KM.1 is correct with labelled state coordinates and the response process
indexed by all deterministic fields.  For a fixed coordinate `x`, the
pinning field that is zero at `x` and `-M` elsewhere satisfies

```math
V_G(u^{x,M})\longrightarrow G_x
```

almost surely.  Taking the vector of such fields for finitely many `x`
recovers the joint law of every coordinate subvector.  A Gaussian law is
then uniquely determined by its mean and covariance.  Conversely, mean and
covariance determine the Gaussian vector and hence every measurable response.

Two scope clauses matter.

* If state labels are quotiented by permutations, the corresponding minimal
  parameter is likewise only defined up to that relabelling.  The draft uses
  labelled fields, so its literal statement is fine.
* This is minimality for the complete **ensemble law**.  Neither covariance
  nor the GREM hull reconstructs a quenched realization.  The draft states
  this distinction correctly and does not evade it through self-averaging.

The REM formula and the tree recursion are direct consequences of
independence.  The claimed `O(k)` GREM parameter count requires, as the draft
does say, that the rooted hierarchy, branching data, and state labels are
declared side information.

## 9. Max-plus cyclicity scope

For `K in R^{Q times Q}`, every transition has finite weight, so the
associated directed graph is complete and hence irreducible.  The standard
irreducible max-plus cyclicity theorem supplies a transient `T`, maximum
cycle mean `lambda`, and an admissible eventual period `gamma` (for example,
a suitable period derived from the critical components) such that

```math
K^{\odot(t+\gamma)}
=K^{\odot t}+\gamma\lambda\mathbf 1
```

for all sufficiently large `t`.  Thus the formal equation in the draft is
correct.  To avoid convention dependence, `gamma` should be called an
admissible ultimate period, or the draft should define "critical cyclicity"
as the relevant common period when the critical graph has several
components.

The finite representation of the complete length family consists of
`lambda`, the eventual residue kernels, **and the finite transient**.  The
draft's phrase "after the transient" is accurate; it should not be read as
discarding the initial powers.  If forbidden transitions `-infinity` are
allowed, irreducibility must be separately assumed (or the reducible theory
used).  If the factors vary with time, fixed `Q^2` boundary-kernel closure
still holds but cyclicity does not.  Queries that expose internal vertices
also lie outside the endpoint interface.  These latter two exclusions are
already stated.

## 10. Director-level assessment

### What is substantive

The vector-balancing package establishes a real statement that is stronger
than "the roof answers roof queries": under repeated additive composition,
convexification approximately answers **nonlinear** uniformly Lipschitz
aggregate extrema, and the error is charged to the effective feature rank
rather than the number of factors.  The `A/B` example identifies a matching
qualitative failure regime.  This supplies a nontrivial model outside the
motivating signing problem with a composable state strictly smaller than its
enumerated extremal landscape.

### What is classical repackaging

The proof mechanism is precisely Shapley--Folkman/floating variables, and the
trellis and Gaussian mechanisms are classical.  The response terminology
organizes them usefully but does not make the underlying theorems new.  The
strongest justified classification is therefore:

* **mathematically generative validation:** yes, in the limited sense that a
  known theorem yields a new declared-query consequence and a sharp stopping
  example;
* **new mathematical mechanism:** no;
* **automatic general compression theorem:** no, because response-body
  description complexity and effective rank still need separate control.

This is more than vocabulary, but it is best presented as a Level-2 synthesis
with one theorem-level application of classical machinery, rather than as a
standalone Level-3 theory result.

### Required corrections before promotion

1. Strengthen or annotate KM.3 with the effective rank
   `r=dim span union_i(E_i-E_i)`; describe the obstruction in terms of this
   rank rather than ambient dimension.
2. Qualify "exactly composable, strict sub-landscape state": this is proved
   for vector balancing through its generator representation and for other
   families whose convex bodies have controlled descriptions, not for
   arbitrary compact component response sets.
3. State explicitly that the Lipschitz constant and component diameters are
   uniform on the normalization used.
4. In the lifted setting, count the energy coordinate unless an affine
   relation removes it.
5. Treat `gamma` as a specified admissible eventual period and retain the
   transient when describing the complete max-plus family.

With these changes, no theorem in the draft needs to be withdrawn.
