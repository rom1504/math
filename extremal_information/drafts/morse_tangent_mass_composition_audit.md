# Audit: Morse tangent-mass composition

Scope: independent audit of
`drafts/morse_tangent_mass_composition.md` and
`experiments/verify_morse_tangent_mass.py`.  I checked the global and local
hypotheses, the lattice count, both powers of `n`, the Gaussian amplitude,
uniform-query and repeated-composition claims, the Vandermonde and
multinomial applications, the quartic falsifier, the information content of
the proposed carrier, and the literature boundary.  I made no edits to either
audited file.

## Verdict

**REPAIR.**  The two analytic conclusions are correct as **fixed-query
standard-lattice discrete Laplace estimates**, after reading (TM.5) as an
eventual statement along the admissible integers.  In the coordinates used in
the draft, the exponent `alpha+beta+d/2` and the constant

```math
{(2\pi)^{d/2}a(x_z)b(z-x_z)\over\sqrt{\det J_z}}
```

are exactly right.

The stronger interpretation is not yet established.  In particular:

1. the full binomial and multinomial arrays do not satisfy the global
   `n^{-d/2}` envelope (TM.1) at boundary types, so the stated corollary needs
   a localization lemma;
2. TM.2 is pointwise and does not prove the uniform output expansion needed
   to feed its result into another application of TM.2;
3. `(f,alpha,a)` contains two arbitrary continuum functions.  It is a finite
   *list of fields*, not a finite-dimensional, finite-cardinality, or
   complexity-controlled state;
4. no integer recovery or state-minimality result is proved; and
5. the exact formula has unit lattice density only because the theorem assumes
   the full standard lattice.  A general tangent lattice contributes its
   covolume.

Thus the draft succeeds as a clean analytic repair of the Vandermonde
exponent on a declared Morse fibre.  It does not yet meet the stronger claim
of a strict finite decorated response state closed on a nontrivial Morse
class.

## 1. The order estimate and the `n^(d/2)` count

For a fixed `z`, the proof of TM.1 is sound under its literal global
assumptions.  Feasibility is equivalent to `k/n in D_z`.  The upper quadratic
bound gives

```math
\sum_{k\in\mathbb Z^d}
 e^{-c\lVert k-nx_z\rVert^2/n}=O(n^{d/2})
```

uniformly in the shift `nx_z`: the sum factorizes into `d` shifted
one-dimensional Gaussian sums, each `O(sqrt(n))`.

For the lower bound, positive distance from `x_z` to `partial D_z` puts a
fixed Euclidean ball around `x_z` inside `D_z`.  Once `n^{-1/2}<eta`, every
lattice point with

```math
\lVert k-nx_z\rVert\le \sqrt n
```

is feasible.  A shifted radius-`sqrt(n)` ball contains
`Omega(n^(d/2))` integer points, uniformly in its centre, and (TM.4) bounds
each associated exponential below by `e^(nh(z)-C)`.  Multiplication by the
two input prefactors proves

```math
\log C_n(nz)
=nh(z)+(alpha+beta+d/2)\log n+O(1).
```

Consequently (TM.6) follows at speeds `(n,log n)`.

Three statement repairs are needed.

- TM.5 holds for **all sufficiently large** `n` such that `nz in Z^d`, not
  necessarily for every `n`.  Small admissible `n` can have an empty
  decomposition fibre even though interior lattice points exist eventually.
- A fixed `z` has infinitely many admissible integers only when its
  coordinates are rational.  A continuum response roof should instead be
  stated uniformly for grid queries `z=m/n` (or for rounded `z_n -> z`).
- “Uniformly Morse” at a single `z` is misleading terminology.  TM.4 is a
  two-sided quadratic-exposure condition.  Uniformity across a query set
  requires common `eta,c,C`, not merely that every query separately has a
  nondegenerate maximizer.

The comment after (TM.4) is correct with an omitted regularity qualification:
a unique interior maximum, a negative-definite Hessian and local `C^2`
regularity give the local quadratic bounds; compactness and uniqueness turn
the upper bound into a global one.  A bare “unique nondegenerate maximum on a
compact smooth domain” is enough only when nondegeneracy and smoothness refer
to the objective near that maximum.  Smoothness of the domain itself is not
what supplies the inequalities.

## 2. Global versus local assumptions

TM.1 deliberately assumes (TM.1) over the entire feasible fibre, because its
upper proof uses that envelope for every summand.  The phrase “throughout the
portions of the domains used below” is therefore too indefinite.  One of the
following precise alternatives is needed:

1. impose (TM.1) globally on every feasible lattice point; or
2. impose the two-sided `n^alpha` and `n^beta` asymptotics only in a fixed
   saddle neighbourhood, together with a cruder global upper bound and a
   uniform exponential gap outside that neighbourhood.

The second formulation is the one needed for the advertised counting models.
For example, with `H` the binary entropy,

```math
{n\choose 0}=1,
\qquad n^{-1/2}e^{nH(0)}=n^{-1/2}.
```

No fixed `K` makes the upper half of (TM.1) hold at `k=0` with
`alpha=-1/2`.  The same face-dimension change occurs for multinomial types at
the simplex boundary.  Restricting `D_A,D_B` to an interior compact set makes
TM.1 applicable only to a **truncated** convolution, not immediately to the
full Vandermonde or type convolution.

The standard repair is short.  Use uniform Stirling asymptotics in a fixed
interior neighbourhood of the saddle.  Away from that neighbourhood, strict
concavity supplies a uniform entropy gap, while the elementary global bound
`binomial(n,k) <= e^(nH(k/n))` (and its method-of-types analogue) makes the
discarded polynomially many terms exponentially negligible.  This
localization step should be stated before applying TM.1/TM.2 to the full
arrays.

## 3. Exact tangent amplitude and determinant

In the stated standard coordinates the constant in TM.14 is correct.  Put

```math
y={k-nx_z\over\sqrt n}.
```

The rescaled grid is a translate of `n^(-1/2)Z^d`, so its point density is
`n^(d/2)`.  Taylor expansion gives

```math
n(F_z(k/n)-h(z))=-\tfrac12 y^T J_z y+o(1),
```

and hence

```math
n^{-d/2}e^{-nh(z)}n^{-alpha-beta}C_n(nz)
\longrightarrow
a(x_z)b(z-x_z)
\int_{\mathbb R^d}e^{-y^TJ_zy/2}\,dy.
```

The integral is `(2pi)^(d/2)/sqrt(det J_z)`.  There is no missing factor of
`n`, `2`, or `2pi`.

For a general full-rank tangent lattice `A Z^d`, however, the constant is

```math
{(2\pi)^{d/2}a(x_z)b(z-x_z)
 \over |\det A|\sqrt{\det J_z}}.
```

Thus TM.14 is not coordinate-free: it has silently fixed lattice covolume
one.  This matters when types are described in the ambient simplex
hyperplane, or when parity/sub-lattice restrictions are allowed.  Eliminating
one multinomial count gives a standard `Z^(q-1)` coordinate lattice and makes
TM.14 literal.  In an orthonormal basis of the sum-zero hyperplane, the same
answer is obtained only after including the tangent lattice covolume.  A true
“tangent-mass decoration” should either retain this lattice datum or declare
the standard full-grid class permanently.

The local amplitude assumption is otherwise adequate because TM.2 inherits
the global bound from TM.1.  It should say explicitly that both `o(1)` errors
are uniform relative errors on a fixed saddle neighbourhood.  The proof's
“growing bounded window” also needs a rate: under `C^3`, one may take
`R_n -> infinity` with `R_n^3/sqrt(n) -> 0`, or avoid the phrase by first
working on fixed rescaled balls and then sending their radius to infinity.

For finitely many maxima, summing the displayed constants is correct provided
the maxima are uniformly separated and each has a uniform positive Hessian
gap.  This is an extension, not a consequence of the unique-maximizer theorem
as currently stated.  Coalescing saddles and changing saddle counts destroy
the required uniformity.

## 4. Uniform queries and repeated composition

The order-level uniform claim after TM.1 can be made correct.  For a compact
query set `Q`, a sufficient statement quantifies over all pairs `(n,z)` with
`z in Q`, `nz in Z^d`, and
requires:

- one global input constant `K`;
- common `eta,c,C` in (TM.4);
- a compact grid-query set whose saddles stay uniformly interior; and
- the same dimensions and full tangent lattice at every query.

The same lattice-sum proof then gives one output constant `K_C` and therefore
a uniform TM.1-type envelope for `C_n` on `Q`.

TM.2 does **not** yet prove analogous reusable exact-amplitude closure.  It
only proves `C_n(nz) ~ ...` for one fixed `z`.  To use `C_n` as an input to a
second convolution, one needs the output expansion uniformly for `z` in a
neighbourhood of every future saddle.  That requires uniform input
remainders, uniform derivative bounds, Hessian eigenvalues bounded away from
zero, and a smooth uniformly interior saddle map `z -> x_z`.  Under such
assumptions the implicit-function theorem and a parameter-uniform Laplace
argument do give the advertised output field `c(z)`, but that theorem is not
written here.

Nor is an actual closed class exhibited.  “Any class closed under these
conditions” is formally true but does not prove that a useful class exists.
A defensible example would be a finite-parameter family of uniformly strongly
concave `C^3` profiles on convex domains, with compact query sets bounded away
from every Minkowski-sum boundary and amplitudes bounded above and below.
The draft should either prove preservation for such a class or call TM.15 a
one-step composition formula rather than a reusable algebra.

Associativity of the formal Gaussian rule is plausible and classical: in a
three-factor convolution the two sequential Hessian determinants combine by
a Schur-complement identity into the determinant of the full `2d`-dimensional
saddle.  It is not checked in the draft, and analytic associativity still
depends on the missing uniform closure hypotheses.

## 5. Vandermonde and multinomial checks

### Vandermonde

The exponent calculation is correct.  On a subsequence with `np` integral,

```math
{n\choose np}
\sim {e^{nH(p)}\over\sqrt{2\pi n p(1-p)}}.
```

For the convolution target `2np`, the saddle is `x_z=p`,

```math
a(p)={1\over\sqrt{2\pi p(1-p)}},
\qquad
J_z={2\over p(1-p)}.
```

TM.14 therefore gives the order-one coefficient

```math
{\sqrt{2\pi}\,a(p)^2\over\sqrt{J_z}}
={1\over 2\sqrt{\pi p(1-p)}},
```

which agrees exactly with Stirling for `{2n choose 2np}` when the power is
written as `n^(-1/2)`.  At `p=1/2` this is `1/sqrt(pi)`, the constant checked
by the experiment.

The draft must add `np in Z` (or specify rounded types and evaluate entropy at
the rounded type).  It must also insert the localization argument from
Section 2; uniform Stirling on a compact interior interval alone does not
verify global (TM.1).

### Multinomial types

The intrinsic convolution dimension is indeed

```math
d=q-1,
```

because one of the `q` counts is fixed by their sum.  In the first `q-1`
count coordinates,

```math
{n!\over\prod_{i=1}^q (np_i)!}
\sim
(2\pi n)^{-d/2}
\left(\prod_{i=1}^q p_i\right)^{-1/2}e^{nH(p)}.
```

Two child powers `-d/2` plus a `d`-dimensional saddle mass `+d/2` return
`-d/2`.  Thus lines 190--201 have the correct dimension and exponent.

Again, all `np_i` must be integral, the query types must remain in a compact
simplex interior, and boundary terms require localization rather than the
literal global prefactor assumption.  Any exact determinant statement must
also specify whether the Hessian is taken in first-`q-1` integer coordinates
or in an orthonormal tangent basis together with the corresponding lattice
covolume.

## 6. Quartic falsifier

The quartic example is correct.  Since

```math
-n(k/n)^4=-k^4/n^3,
```

the active window is `|k|=Theta(n^(3/4))`, and the Riemann sum gives the
sharper asymptotic

```math
n^{-3/4}\sum_{|k|\le n}e^{-k^4/n^3}
\longrightarrow
\int_{\mathbb R}e^{-u^4}\,du
={\Gamma(1/4)\over2}.
```

It therefore genuinely falsifies a universal `d/2` correction outside the
quadratic-exposure class.  A flat interval similarly has `Theta(n)` lattice
mass.

This example does not threaten TM.1: `F(x)=-x^4` fails its quadratic upper
exposure near zero.  It does show that the proposed “finite stratification by
saddle type” is a research direction, not a consequence.  Without a declared
finite menu of degeneracies, arbitrary orders of flatness and non-polynomial
tangent profiles need not admit a finite stratification.

## 7. What the verification script establishes

The script runs successfully and reports 15 checks.  Its formulas and finite
calculations are internally correct.

- The log-gamma Vandermonde sum agrees with `{2n choose n}` to the stated
  tolerance for `n <= 640`.
- The difference between the full convolution and its largest summand trends
  toward `(1/2)log n`.
- The central-binomial output is within two percent of
  `4^n/sqrt(pi n)` on the tested range.
- The quartic raw slopes trend toward `3/4` and are already above `0.70`.

These are smoke tests, not a verification of the general theorems.  In
particular the script does not test:

- `d>1` or the multinomial `q-1` dimension;
- the determinant assembled from the two input amplitudes and Hessian;
- an off-centre or uniform family of queries;
- a second composition using the first output as an input;
- multiple saddles, boundary localization, or a non-unit tangent lattice;
- any integer recovery or finite-state encoding claim.

The quartic statistic `log(sum)/log(n)` converges slowly because it retains
the additive log of the limiting integral; its last printed value is about
`0.809`, not especially close to `0.75`.  A sharper regression would test the
normalized sum against `Gamma(1/4)/2`, or use dyadic log differences.  For
TM.14, a direct test should form a nontrivial Gaussian convolution (and a
two-dimensional example), compute `J_z`, `a`, and `b`, and compare to the
predicted constant.  The present central-binomial amplitude check uses the
known output Stirling formula and does not independently exercise the
determinant law.

## 8. Is this a strict finite decorated state?

No, not under the ordinary meanings of finite state or finite-dimensional
carrier.

At logarithmic precision the proposed object is `(f,alpha)` on a continuum
domain; at order-one precision it is `(f,alpha,a)`.  Both `f` and `a` may be
arbitrary continuous functions.  The amplitude field can carry unbounded
information, and no metric-entropy, parameter-count, finite alphabet, or
approximation theorem is supplied.  It is smaller than the entire
scale-dependent atomic array `(A_n)` because it discards finite-`n`
fluctuations, but that fact alone is not a strict finite-state reduction.

Moreover, the fixed full lattice is itself structural information.  Once
descriptor grids, sublattices, or face lattices vary, tangent density/coset
data must also be retained.  Omitting it does not make the carrier smaller; it
restricts the model class to the unit-density grid.

The theorem assumes positive real arrays.  It gives no recovery construction
from `(f,alpha,a)` to integer multiplicities, much less recovery inside codes,
graphs, or another constrained model.  An abstract recovery after adding a
large common leading normalization and flooring is likely easy on a compact
full grid, but it would still not establish constrained realization or a
finite encoding.

The Vandermonde example proves that a bare pointwise lexicographic maximum is
insufficient.  It does **not** prove that `(f,alpha,a)` is minimal, that every
sufficient carrier must store an amplitude field, or that no different
finite statistic works on a narrower family.  Accordingly “minimal
compositional repair” and “finite tangent-mass decoration” are currently
stronger than the results.

There are two honest ways to repair the scope.

1. Call this a **finite-depth functional decoration**: a finite list of
   limiting fields on a fixed standard lattice, with no finite-state claim.
2. For a strict finite carrier, restrict `f` and `a` to an explicitly
   finite-parameter family preserved by convolution (quadratic/Gaussian
   profiles are the cleanest starting point), or restrict to a fixed finite
   query alphabet, and then prove an encoding/recovery bound.

## 9. Classical versus new

The following ingredients are classical.

- Gaussian upper/lower estimates for lattice sums at a nondegenerate saddle;
- the `n^(d/2)` local lattice mass and the
  `(2pi)^(d/2)/sqrt(det J)` constant;
- the covolume correction for non-unit lattices;
- summing contributions from finitely many separated saddles;
- Stirling, Vandermonde, multinomial type asymptotics; and
- the change from `n^(1/2)` to `n^(3/4)` at a quartic saddle.

The cited [Moran paper](https://doi.org/10.2307/3213083) is a genuine 1979
coefficient/local-lattice-CLT reference and supports the claimed historical
boundary.  A recent explicit formulation is Hughes--Helton--Schlosser,
[The discrete Laplace asymptotic method](https://arxiv.org/abs/2509.16420),
whose standard theorem also displays the Hessian and lattice-covolume factors.

What is specific to this project is the **interpretation** of `d/2` as the
defect in the earlier pointwise two-speed roof and the packaging of the
classical formula as the candidate rule (TM.15).  That is a useful conceptual
bridge.  The draft does not establish a new analytic Laplace theorem, a
minimality theorem, a strict finite-state theorem, or a closed/recoverable
Morse carrier.  Novelty should therefore be claimed for the response-carrier
framing and scoped diagnosis, not for the asymptotic or for proven
minimal/finite closure.

## 10. Recommended scope repairs

1. Change the status to “rigorous fixed-query discrete-Laplace theorem;
   uniform functional-carrier extension pending,” or prove the missing
   parameter-uniform theorem.
2. Quantify TM.5 over sufficiently large admissible `n`, and formulate
   continuum queries through grid points `z=m/n` or rounded sequences.
3. Replace the vague domain phrase in TM.1 by explicit global assumptions, or
   state a local-saddle plus global-exponential-tail version.
4. Add the localization lemma before invoking binomial or multinomial
   examples; add all integrality conditions.
5. State that TM.14 uses the unit-density standard lattice.  Add a covolume
   factor, or include a fixed lattice in the declared class.
6. Prove a parameter-uniform TM.2 with a smooth saddle map before calling
   `(f,alpha,a)` a reusable algebra.  Treat multiple saddles in a separate
   theorem with uniform separation.
7. Exhibit one genuinely closed class under repeated convolution and check
   associativity/uniform output errors.
8. Define “finite.”  If it means only a finite tuple of continuum fields, say
   so.  If strict finite-dimensional/state reduction is intended, impose a
   finite-parameter family and give a recovery or complexity result.
9. Remove “minimal” unless accompanied by a lower-bound/equivalence theorem.
10. Describe finite saddle stratification as a proposed extension, not the
    uniquely legitimate next step.
11. Extend the experiment with an off-centre determinant check, a
    `q=3` multinomial check, a normalized quartic constant, and a lattice
    covolume example.

With these repairs, TM.1 and TM.2 are strong, clean supporting lemmas for the
multi-speed program.  Without them, the analytic core passes but the promoted
finite-state and reusable-closure interpretation does not.
