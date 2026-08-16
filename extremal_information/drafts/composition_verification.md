# Independent verification of tensor-span composition

Date: 2026-08-16.

Scope: Section 4 of `sufficiency_axioms_report.md`, especially Theorem 4.1,
(4.7)--(4.9), and the associativity/cocycle discussion.  This report does not
audit the broader theory or alter the main project state.

## Verdict

Theorem 4.1 and inequalities (4.7)--(4.9) are correct for the finite,
common-interface setting stated in the report.  No finite counterexample was
found, and the proof reduces to two exact facts: separately affine functions
commute with independent barycentres, and taking a maximum is 1-Lipschitz in
the uniform norm.

The associativity and cocycle identities are also the right pointwise
conditions for exact bracket-independent composition.  Three qualifications
should be added:

1. the relevant feature domain must be closed under the proposed composition
   (at roof level, it is safest to use a convex domain);
2. the associated bilinear algebra is naturally defined on
   `span{(1,u):u in D}`, not automatically on all of `R^(d+1)`;
3. normalized magnetizations and normalized energies require size-dependent
   weighted composition laws, or a size/mass coordinate in the state.

The result is non-tautological: it identifies a genuine algebra of future
operations under which upper concavification loses no optimum-value
information.  It is not, by itself, a compression theorem; a roof can still
have exponentially large description complexity.

## 1. Audit of Theorem 4.1

Let `P_i=conv(phi_i(Omega_i))`.  For fixed parent field `eta`, combine the
cross energy and queried parent feature as

```math
G_\eta(u,v)=c_\eta+p_\eta^T u+q_\eta^T v+u^T B_\eta v.
```

For arbitrary `u in P_1` and `v in P_2`, choose roof-attaining laws
`lambda,mu` with the prescribed feature means.  Finiteness guarantees such
laws.  Under the product law,

```math
E[H_1(X)+H_2(Y)+G_\eta(phi_1(X),phi_2(Y))]
=\widehat H_1(u)+\widehat H_2(v)+G_\eta(u,v),
```

because `G_eta` is affine in either variable separately and `X,Y` are
independent.  An expectation cannot exceed the largest pure queried parent
energy.  Conversely, choosing feature means of any pure pair gives a value
at least that pair's energy, since a roof dominates every generator in its
feature fibre.  This proves (4.5).

There is no hidden assumption that parent roof mixtures are product laws.
Product laws prove the response identity; concave Fenchel--Moreau inversion
of the complete response function then recovers all correlated parent
mixtures in the parent roof.  This distinction is important and the report
handles it correctly.

The argument needs either:

- finite state spaces, as currently stated; or
- for an infinite extension, compact feature/energy data and proper
  upper-semicontinuous concave roofs, with attainment replaced by an
  approximation argument.

The maps `C,F`, the child feature maps, and the energy normalization must be
common to the landscapes being compared.  Affine relations among feature
coordinates cause no problem, but the bi-affine representation is then not
unique and all statements should be understood on the realized affine hull.

A random stress test of 2--5-state scalar children, integer
features/energies, and integer bi-affine coefficients checked 2,000 cases and
found no product-relaxation value above the pure maximum.  This is only a
sanity check; the preceding argument is the proof.

## 2. Audit of the data-processing inequalities

For (4.7), order the maximum first by the second pure state:

```math
V_{P(H)}(\eta)-c_\eta
=\max_y\{H_2(y)+q_\eta^Tv
          +V_{H_1}(p_\eta+B_\eta v)\}.
```

Replacing `H_1` by `G_1` costs at most
`d_{Theta_1(eta)}(H_1,G_1)`.  Reorder the resulting maximum by the first pure
state and replace `H_2` by `G_2`; this costs at most the analogous second
term.  The triangle inequality gives (4.7).  It is enough that the local
direction sets in (4.6) use **pure** feature values.  Enlarging those sets to
the feature polytopes is valid but can give a strictly larger metric.

For (4.8), if `f,g` are closed concave roofs on the same compact polytope,
the response transform and its inverse are both contractions in vertical
sup norm:

```math
\|V_f-V_g\|_\infty\le\|f-g\|_\infty,
\qquad
\|f-g\|_\infty\le\|V_f-V_g\|_\infty.
```

The second inequality follows by applying the infimum formula to both
responses.  Hence equality holds.  The common-polytope hypothesis is
essential: singleton features at `0` and `1`, both with energy zero, have
identical vertical roof values on their respective domains but responses
`0` and `theta`, whose uniform distance over all fields is infinite.

Applying (4.7) uniformly in `eta`, then (4.8) to the common parent feature
polytope, proves (4.9).  For a weighted parent energy
`alpha_1 H_1+alpha_2 H_2+C`, the correct right side is
`|alpha_1| epsilon_1+|alpha_2| epsilon_2`, not simply the sum.  This matters
when passing to normalized asymptotic energies.

The induction claim is correct for a fixed composition tree and common
operations.  If randomized sketches are used, the corresponding expected
distortion statement additionally requires access to decoders for both
child summaries (product decoders suffice); it does not require statistically
independent approximation errors.

## 3. Associativity, cocycles, and normalization

For pure feature values `u,v,w`, the two bracketings have features

```math
F(F(u,v),w),\qquad F(u,F(v,w)),
```

and accumulated cross energies

```math
C(u,v)+C(F(u,v),w),
\qquad
C(v,w)+C(u,F(v,w)).
```

Thus (4.10) and (4.11) are sufficient, and are necessary if exact equality
of feature and energy is required for every realizable triple.  They need
only hold on the realizable domain; requiring them on a larger ambient set is
a convenient stronger condition.  For equality only after passing to a
roof, they need not be logically necessary because invisible lower states or
a feature-preserving relabelling can conceal pointwise disagreement.

Homogenize a bi-affine `F` by defining a bilinear product on augmented
coordinates whose first coordinate is `alpha beta` and whose second
coordinate is the homogenization of `F`.  Equation (4.10) makes this product
associative on

```math
E=span\{(1,u):u\in D\}.
```

It yields an associative algebra on all of `R^(d+1)` only if these augmented
points span the whole space, equivalently if `D` has full affine span.  A
small warning example is `D={0}` with an ambient formula `F(u,v)=2u`:
(4.10) holds on every realizable triple, while the homogenized product is not
associative away from the one-dimensional span of `(1,0)`.  There is no
failure on the actual feature algebra; only the ambient claim would be too
strong.

For iterative roof composition, take `D` convex (or replace it by its convex
hull) and require `F(D,D) subset D`.  Separate affinity then extends identities
verified on pure feature values to their convex hulls.

Normalization is the other important caveat.  Raw magnetization uses
`F(u,v)=u+v`, which is associative on an additive domain and for which
`C(u,v)=beta uv` obeys the cocycle identity.  The naively normalized rule
`F(u,v)=(u+v)/2` is not associative:

```math
F(F(u,v),w)=(u+v+2w)/4,
\qquad
F(u,F(v,w))=(2u+v+w)/4.
```

For blocks of sizes `m,n`, the correct rule is
`F_{m,n}(u,v)=(mu+nv)/(m+n)`.  Its associativity and the corresponding energy
cocycle are size-indexed.  One must either retain the sizes externally or
include mass as a state coordinate.  Analogous weights are required for
energies normalized by a power of system size.

## 4. Is this a genuine composition theorem?

Yes, with a limited but precise scope.  Defining the roof makes one-query
sufficiency almost formal; Theorem 4.1 is stronger.  It proves that an entire
class of new parent queries can be answered after discarding all child
points below the upper concave envelopes, and that the discard can be
iterated whenever the feature algebra closes.  The nonlinear-feature example
(4.13)--(4.14) shows that this closure is not automatic and marks a real
boundary of the theorem.

What the theorem does **not** establish is equally important:

- fixed feature dimension does not bound the number of roof faces;
- it gives no rate--distortion or all-order realization bound;
- a full-spin interface may encode the entire landscape despite having only
  `n` coordinates;
- a new feature outside the tensor span can resurrect discarded states.

It should therefore be called an exact composition/closure theorem, not yet
an extremal compression theorem.

## 5. Application outside the motivating signing problem

Consider a finite-spin mean-field model.  Site `i` has a finite state set,
an arbitrary local energy `h_i(s)`, and a bounded feature
`phi(s) in Z^d`.  For a symmetric matrix `J`, define

```math
E(s_1,\ldots,s_N)
=\sum_i h_i(s_i)
 +\sum_{i<j}\phi(s_i)^T J\phi(s_j).
```

For a block, retain the upper roof over its total feature
`u=sum_i phi(s_i)`.  Merging two blocks uses

```math
F(u,v)=u+v,
\qquad C(u,v)=u^T Jv.
```

Both are bi-affine, `F` is associative, and `C` satisfies (4.11).  Theorem
4.1 therefore gives an exact, bracket-independent dynamic program for the
ground-state value under arbitrary hierarchical merging.  If `d` and the
single-site feature alphabet are fixed, a block of size `N` has only
polynomially many attainable total features; the roof can be represented by
at most that many lifted generators even though the spin landscape has
exponentially many states.  This covers finite-rank Curie--Weiss and
mean-field Potts-type ground-state models with site-dependent local fields.

This application demonstrates real compression and an exact composition
algebra outside the quadratic-signing problem.  It also exposes why the same
theorem may fail to compress a full Boolean bridge: there the natural feature
polytope has exponentially many exposed pure states, so fixed coordinate
dimension does not imply small extremal information.
