# Adversarial audit: dense-transversal composition growth

**Scope.**  This is an independent audit of
[`phase3_transversal_composition_growth.md`](phase3_transversal_composition_growth.md).
It reconstructs the five theorem claims from the definitions and separates
the exact state update from the approximate response decoder.  The companion
finite audit is
[`verify_phase3_transversal_composition_audit.py`](../experiments/verify_phase3_transversal_composition_audit.py).

## Verdict

The mathematical cores of TC.1--TC.5 survive the audit.  In particular:

* the constants `8` and `10` in TC.1 follow from the exact count of
  nondegenerate quotient pairs and from a two-letter consolidation overhead;
* the product in TC.2 is well defined on affine-quotient states and its
  decoder has the advertised additive error;
* TC.3 is an exact identity, not only a lower bound;
* the two families in TC.4 have equal source count and isomorphic linear
  dependence data, while their radii are `D` and `D/2`; and
* the affine-circuit-rank decoder in TC.5 has error at most `(r+1)/2`
  against arbitrary appended futures, independently of the number of exact
  linear graph sources.

Two scope qualifications should be made in the main note.

1. The product in TC.2 is an **exact algebra of affine feature states**, not
   an exact algebra of word metrics or future responses.  Its response
   decoder remains additive-`ell` approximate.  Calling it simply an
   "exact algebra" without this qualifier risks contradicting TC.4.
2. Exact support reconstruction is automatic for the very strong query
   class consisting of all appended-fragment radius queries: append every
   nonzero group element except `x`; the resulting radius is one exactly
   when `x` belonged to the original support, and otherwise is two (in
   ambient dimension at least two).  Thus the exact-semiltattice observation
   is correct, but by itself is not a new mixed-cycle lower bound.  The new
   content of TC.4 is the **macroscopic approximate** separation already at
   the empty future, despite equal affine state and equal depth.

There is also a presentational qualification: `ell=o(D)` is what makes the
`O(ell)` response error submacroscopic.  It does not by itself imply that the
`D^2+Dk`-bit affine encoding is shorter than an `ell Dk`-bit list of already
linearized maps.  The state is nevertheless polynomial in `(D,k)` and is a
strict quotient of the full rooted response landscape; bit compression of a
particular input representation depends on the parameter regime.

A further scope restriction is essential.  The rank-only TC.5 estimate is a
theorem about unions of **exact linear graph supports**.  If arbitrarily many
nonlinear cycle-contracting transversals are first synchronized one at a
time, TC.1 still contributes `10 ell`.  No theorem in the note replaces that
term by the affine rank of the synchronized centres.  Consequently the
opening and final research judgment should not state the arbitrary-depth
rank law for raw cycle-contracting sources.

## TC.1 reconstruction

Write `N=|Q|`, `e=f+L`, and

```math
A=\sum_{q\in Q}|e(q)|_B.
```

The vector BLR estimate and the fact that the defect vanishes on the
`2N-1` degenerate ordered pairs give

```math
A\le \frac{3(N-1)(N-2)}N.
```

For fixed nonzero `q`, every `a` outside `{0,q}` gives the three-cycle
inequality

```math
|e(q)|_B\le3+|e(a)|_B+|e(q+a)|_B.
```

After summing, both error sums equal `A-|e(q)|_B`, so

```math
N|e(q)|_B\le3(N-2)+2A<9N.
```

Integrality yields `|e(q)|_B<=8`.  This checks the slight improvement from
the generic uniform BLR constant nine.

For response replacement, a word over the linear graph consolidates all of
its graph letters exactly to zero or one graph letter.  The latter can be
changed to its `f`-letter using at most eight basis letters.  Conversely, a
set of `r` `f`-letters with nonzero quotient sum can first be replaced by the
single `f`-letter over that sum and at most `r+1` basis letters.  This costs at
most two more letters, after which synchronization costs at most eight.
The zero-sum case costs no overhead.  Arbitrary future letters are untouched,
so the bounds are genuinely all-context and telescope source by source.

## TC.2 reconstruction

In the quotient by

```math
V=\sum_{i>1}\operatorname{Im}(L_i+L_1),
```

all graph generators have the same image.  Its exact quotient word metric is

```math
d_B(w+L_1(q),V)+\mathbf 1_{q\ne0}.
```

This proves the lower bound.  Any minimizing `v in V` has a representation

```math
v=\sum_{i>1}(L_i+L_1)(q_i).
```

Choosing `q_1=q+sum_{i>1}q_i` realizes the quotient correction using at most
`ell` graph letters, proving the upper bound.

If `L` is perturbed by a map into `V` and `K` by a map into `Z`, then

```math
V+Z+\operatorname{Im}(L+K)
```

is unchanged.  Hence the product is well defined on the quotient classes.
The reference-free family description proves associativity, commutativity,
and idempotence.  The finite audit independently enumerates the quotient
states at `(D,k)=(2,2)` and checks these laws, as well as the decoder after
min-plus convolution with every future at `(D,k)=(2,1)`.

## TC.3 reconstruction

Every shortest word for `t=sum_i b_i` is a subset of the support.  If `R`
is its set of nonbasis letters, its quotient sum is zero.  Once `R` is fixed,
the basis correction is unique and has size

```math
D-\left|\sum_{s\in R}s\right|_B.
```

Minimizing over `R` is exactly the same as subtracting the maximum excess,
which proves TC.3 with no hidden spanning or parity hypothesis.  When each
source contracts its own quotient cycles, positive excess in a union must
come from a subset not contained in any one source.  This is the precise
sense in which it is mixed.

## TC.4 reconstruction

For even `D`, the linear transformation `e_i mapsto t+e_i` is invertible:
if `sum_i c_i(t+e_i)=0`, every coordinate gives `c_i=sum_j c_j`, and evenness
then forces every coefficient to vanish.  Thus the two `D+1` point sets are
linearly isomorphic and both generate state `[0,W]`.

For the sparse family, `(t,1)` costs exactly `D`, and every other element
costs at most `D`.  For the mixed family, let `D=2m`.  If `r` nonzero lifts
are used, their kernel sum has weight `r` for even `r` and `D-r` for odd
`r`.  Against any target kernel vector of weight `m`, the reverse Hamming
triangle inequality makes lift cost plus correction cost at least `m`.
The two one-lift/no-lift constructions in TC.4 give the matching upper
bound.  Hence the radii are exactly `D` and `m`.

Since the data `(D,k,ell,[L,V])` agree, any single estimator based on those
data incurs error at least half the radius gap, namely `D/4`, on one of the
two instances.  This is a valid equal-depth lower bound and establishes that
the `O(ell)` order in the affine decoder cannot be improved uniformly from
affine-state information alone.

## TC.5 reconstruction

Let the affine hull of the exact linear maps in `F` be `L_0+cal U`, and
choose maps `L_j=L_0+A_j` from `F` so that `A_1,...,A_r` is a basis of
`cal U`.  The kernel subspace used in the lower profile is

```math
V({\cal U})=\sum_{A\in\cal U}\operatorname{Im}A
            =\sum_{j=1}^r\operatorname{Im}A_j.
```

The equality on the right is important.  One inclusion is immediate.  For
the other, every `A in cal U` is a binary linear combination of the `A_j`,
so every value `A(q)` lies in the sum of their images.  It follows that any
`v in V(cal U)` can indeed be written with one evaluation per basis map,

```math
v=\sum_{j=1}^r A_j(q_j).
```

With `q_0=q+sum_j q_j`, the `r+1` available graph generators over the
chosen affine basis sum to `(L_0(q)+v,q)`.  This validates the nontrivial
realization step in the upper bound.  Quotienting by `V(cal U)` gives the
matching lower profile.

If `c=r+1`, then pointwise `delta<=lambda<=delta+c`; hence
`h=delta+c/2` satisfies `||h-lambda||_infty<=c/2`.  Min-plus convolution by
the word metric of an arbitrary future and subsequent maximization are both
nonexpansive in uniform norm.  This proves the midpoint all-future claim,
including its half-integer interpretation when `r` is even.

The join formula is exactly the affine hull of a union in the vector space
`Hom(Q,W)`.  It is independent of representatives and inherits all
semilattice laws.  This is, again, an exact update of the **affine-subspace
feature**, not of the response profile.  The finite audit enumerates all
affine subspaces of `F_2^3`, checks the join laws and representative
independence, and checks the rank decoder against every future in a small
Cayley instance.

Finally, TC.4 has affine direction rank `r=D` in both families and gives an
empty-future radius gap `D/2`.  Any common affine-state estimate therefore
errs by at least `D/4=r/4` on one family.  This validates both the order of
the TC.5 error and the claim that source count does not repair the state.

## Cross-audit of the mixed-circuit continuation

The exact gluing statement MC.4 in
[`phase3_mixed_circuit_hierarchy.md`](phase3_mixed_circuit_hierarchy.md)
also survives a direct sequence check.  If `U_j` is the quotient span of
fragment `j`, then compatible changes of local representatives form
`bigoplus_j Hom(U_j,W)`.  Quotienting by restrictions of one global shear
gives

```math
\left(\bigoplus_j\operatorname{Hom}(U_j,W)\right)
 /\operatorname{Hom}(U,W)
\cong
\operatorname{Hom}\!\left(\ker(\bigoplus_jU_j\to U),W\right).
```

The kernel is `Z/Z_loc` of dimension `kappa`, so the labeled gluing count is
exactly `2^(D kappa)`.  This counts gauge information, not automatically the
minimal information for one scalar response, a limitation the source note
states correctly.

MC.6 correctly converts one mixed channel into response entropy.  With two
singleton fragments over the same nonzero quotient, `kappa=1` and the rooted
kernel profile is `F_v(u)=min(|u|,2+|u+v|)`.  A constant-rate Hamming packing
of weight-at-least-`D/2` vectors and mutual distance `D/4-O(1)` is separated
by the query `u=v`, giving `Omega(D)` bits for distortion below `D/8`.
This is a lower bound for the rooted endpoint-query family; it should not be
silently restated as a radius-only lower bound.

## What the result teaches

The surviving law is not that affine quotients form an exact response
algebra.  It is the more precise dichotomy:

* **within-source cycle control** gives constant synchronization and a
  bounded-depth affine quotient; but
* **cross-source zero-sum cycles** are newly generated under union, and their
  excess is exactly the loss at the antipodal kernel target.

That is a theorem-level mechanism for feature-algebra growth.  Any proposed
sublinear-depth continuation must retain enough weighted mixed-cycle data to
control this excess; exact all-future sufficiency is too strong because it
already separates every support.
