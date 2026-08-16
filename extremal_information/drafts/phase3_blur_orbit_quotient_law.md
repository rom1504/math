# Power-orbit and bounded-repair laws for reusable response charts

**Status.** Axiomatizer theorem draft.  The algebraic identities below are
self-contained.  They isolate one common mechanism behind the tropical
subgroup chart and the fixed-flat matroid chart, and also explain why the
formally analogous Minkowski blur cannot be reused at arbitrary depth.  The
scope is deliberately narrower than all congruences: the theorem
characterizes **fixed-element blur charts** and supplies a geometric,
fiber-free test for their response distortion.

## 1. Metric-monoid setup

Let `(M,star,1)` be a commutative monoid equipped with an extended
pseudometric `d` for which every translation is nonexpansive:

```math
d(x\star z,y\star z)\le d(x,y).                 \tag{BOQ.1}
```

Fix `b in M`, write `b^(star k)` for its `k`-th power, and define the blur

```math
P_bx=x\star b.                                   \tag{BOQ.2}
```

The metric may be infinite on some pairs.  This is useful for min-plus
indicator kernels.  All conclusions below remain valid with value
`+infinity`.

### Theorem BOQ.1 (exact power-orbit law)

For every `m>=1`,

```math
\sup_{x_1,\ldots,x_m\in M}
d\left(
 (P_bx_1)\star\cdots\star(P_bx_m),
 P_b(x_1\star\cdots\star x_m)
\right)
=d(b^{\star m},b).                              \tag{BOQ.3}
```

More generally, a composition program which multiplies exact inputs and
applies `P_b` a total of `k>=1` times, at arbitrary leaves or internal
nodes, outputs

```math
(x_1\star\cdots\star x_m)\star b^{\star k}.     \tag{BOQ.4}
```

Its worst-case distance from the once-blurred exact product is therefore
exactly `d(b^(star k),b)`.  Consequently the sharp arbitrary-depth algebra
defect of the chart is

```math
\boxed{D_{\rm orb}(b)=\sup_{k\ge1}d(b^{\star k},b).}             \tag{BOQ.5}
```

The following are equivalent:

1. `P_b(x star y)=P_bx star P_by` for every `x,y`;
2. `P_b` is idempotent;
3. `b star b=b`.

When these conditions hold, `P_b(M)` is a submonoid with identity `b`,
`P_b` is a multiplicative retraction onto it, and its fibers form a monoid
congruence.  Thus an idempotent blur has zero algebra defect at every depth.

#### Proof

Commutativity and associativity give

```math
(P_bx_1)\star\cdots\star(P_bx_m)
=(x_1\star\cdots\star x_m)\star b^{\star m}.    \tag{BOQ.6}
```

Translation nonexpansiveness bounds the distance in (BOQ.3) by
`d(b^(star m),b)`.  Equality is attained by taking every `x_i=1`.  Every
application of `P_b` contributes one central factor `b`, proving
(BOQ.4)--(BOQ.5) in the same way.

If `b star b=b`, then

```math
P_b(x\star y)=x\star y\star b
=(x\star b)\star(y\star b),                    \tag{BOQ.7}
```

and `P_b(P_bx)=P_bx`.  Conversely, either multiplicativity applied to
`x=y=1`, or idempotence applied to `1`, gives `b star b=b`.  The remaining
claims are immediate. `square`

This is not merely a repeated triangle inequality.  A generic estimate
would pay one local error per use of the chart.  Equation (BOQ.3) says that
all uses collapse to the single algebraic orbit of the discarded element
`b`, and it is sharp before applying any scalar observable.

## 2. A geometric test for observable distortion

The preceding theorem decides algebraic reuse.  A separate issue is whether
the quotient retains an extremal observable.  The next result replaces a
pairwise search over homomorphism fibers by a bounded **repair distance** to
their canonical saturated representatives.

Let `C subseteq M` be closed under `star`, suppose `b star b=b`, and let
`F` be real-valued on `C union P_b(C)`.  Assume:

1. `F` is `L`-Lipschitz for `d` on these points;
2. blur is one-sided for the observable,

   ```math
   F(P_bx)\le F(x)\qquad(x\in C);                \tag{BOQ.8}
   ```

3. the external geometric repair radius is bounded,

   ```math
   \sup_{x\in C}d(x,P_bx)\le D.                 \tag{BOQ.9}
   ```

The metric in (BOQ.9) is part of the model geometry, not the future-response
metric defined from `F`.  In particular, `P_bx`, its distance from `x`, and
the decoder base value `F(P_bx)` must all be well-defined even when the
saturated state is not itself an admissible exact object in `C`.

### Theorem BOQ.2 (bounded-repair idempotent chart)

Under (BOQ.8)--(BOQ.9), the state

```math
\sigma_b(x)=P_bx                                  \tag{BOQ.10}
```

is an exact, associative, summary-only composition state.  Every one of its
fibers has `F`-oscillation at most `LD`, uniformly after every future
context.  More precisely, if
`P_bx=P_by=s`, then

```math
F(x),F(y)\in[F(s),F(s)+LD].                     \tag{BOQ.11}
```

Therefore the decoder

```math
\widehat F(s)=F(s)+{LD\over2}                   \tag{BOQ.12}
```

has uniform error at most `LD/2`, after an arbitrary number of summarized
compositions and with no accumulation in depth.

The constant is sharp **for this fixed chart** whenever one fiber contains two
admissible objects `x_- ,x_+` with a common saturated representative
`P_bx_-=P_bx_+=s` and

```math
F(x_-)=F(s),\qquad
F(x_+)=F(s)+LD.                                  \tag{BOQ.13}
```

#### Proof

Theorem BOQ.1 makes `sigma_b` multiplicative.  If `s=P_bx`, one-sidedness,
Lipschitz continuity, and (BOQ.9) give

```math
0\le F(x)-F(s)\le Ld(x,s)\le LD.               \tag{BOQ.14}
```

This proves (BOQ.11), and midpoint decoding proves (BOQ.12).  For exact
inputs `x_1,...,x_m in C`, their summary-only product is exactly
`P_b(x_1 star ... star x_m)` by BOQ.1.  Since the exact product remains in
`C`, the same one-shot interval (BOQ.14) applies once at the root,
independently of `m` and of the evaluation tree.  Finally, the two values in
(BOQ.13) differ by `LD` while sharing one state, so no scalar decoder can
have worst-case error below `LD/2`. `square`

In particular, if `P_bx=P_by` and `c in C`, then
`P_b(x star c)=P_b(y star c)`; applying (BOQ.11) to `x star c` and
`y star c` bounds their future-response difference by `LD` as well.

The certificate is uniform over observables: once `(M,d,b,C,D)` is fixed,
the same quotient works simultaneously for the entire cone of one-sided
`L`-Lipschitz functions, with only the scalar decoder changing.  This is a
useful formal sense in which the repair geometry is not manufactured from a
single response table.

If the orientation (BOQ.8) is omitted but the same Lipschitz and repair
bounds hold, the argument only places `F(x)` in
`[F(s)-LD,F(s)+LD]`; fiber oscillation may then be `2LD` and the best
guaranteed decoder error is `LD`.  Thus the one-sided extremal order is what
earns the factor two in BOQ.2.  Idempotence and bounded repair have similarly
separate roles: BOQ.1 shows that dropping idempotence exposes the whole power
orbit, while dropping bounded repair leaves CSC.1 with no geometric bound on
fiber oscillation.

### Why BOQ.2 is stronger than invoking CSC.1

CSC.1 says that a *given* congruence is accurate precisely when its fibers
have small observable oscillation.  That is an exact characterization, but
checking it can require comparing the entire response landscape inside
every fiber.

BOQ.2 gives a different, checkable certificate:

```math
\text{idempotent algebraic blur}
+\text{ bounded geometric repair to }P_bx
+\text{ one-sided Lipschitz observable}.        \tag{BOQ.15}
```

Every fiber is automatically star-shaped around the same canonical
saturation `s=P_bx` in the only sense needed by the proof: all of its
observable values lie in the one-sided interval based at `F(s)`.  No
pairwise fiber comparison and no enumeration of future responses is used.
The conclusion is therefore not the assertion "homomorphism plus small
oscillation" with different notation; bounded repair is a sufficient
geometric mechanism which *proves* small oscillation.

The hypotheses are not universal.  CRL.2 shows that arbitrary response
tables can live in monotone union semilattices, so some additional structure
beyond monotonicity and idempotence is indispensable.  A common saturation
with bounded repair is one sufficient form, not a necessary classification.

## 3. Tropical subgroup/coset profiles

Let `G=F_2^w`.  Profiles `f:G->[0,infinity]` compose by min-plus
convolution.  Use the standard extended sup metric

```math
d_\infty(f,g)=\inf\{a\ge0:f\le g+a\text{ and }g\le f+a\},       \tag{BOQ.15a}
```

whose restriction to finite profiles is the usual sup norm.  This avoids
an undefined `infinity-infinity` when indicator kernels occur.  Translation
by a common profile is nonexpansive.  For a subgroup `H<=G`, take

```math
b=\iota_H,
\qquad
b(h)=0\ (h\in H),
\quad b(x)=+\infty\ (x\notin H).                \tag{BOQ.16}
```

Then `b star b=b`, and

```math
P_bf(x)=\min_{h\in H}f(x-h)                    \tag{BOQ.17}
```

is exactly the coset-minimum profile, equivalently a profile on `G/H`.

Fix a coordinate basis and let `H` be spanned by `d` basis vectors.  On the
convolution-closed class of word profiles whose supports contain that
basis, every `f` is one-Lipschitz for Hamming distance.  Hence

```math
0\le f(x)-P_bf(x)\le d                           \tag{BOQ.18}
```

for every `x`: choose a minimizing `h` in (BOQ.17) and use `|h|<=d`.
Thus `||f-P_bf||_infinity<=d`.  With `F(f)=max_x f(x)`, conditions
(BOQ.8)--(BOQ.9) hold with `L=1,D=d`.  Theorem BOQ.2 gives an exact
summary-only convolution algebra on quotient profiles and radius error at
most `d/2` at arbitrary depth.  This recovers TDS.4 from the common
bounded-repair principle.

For a general finite penalty kernel `b`, Theorem BOQ.1 instead gives

```math
D_{\rm orb}(b)=\sup_k\|b^{\star k}-b\|_\infty.
                                                               \tag{BOQ.19}
```

The min-plus powers converge to the shortest-path closure `b_*`, so this is
`||b-b_*||_infinity`.  Thus TDS.2 is the min-plus specialization of the
power-orbit law, while subgroup indicators supply a particularly useful
family of exact idempotent charts.  Other subadditive kernels can also be
min-plus idempotents.

## 4. Fixed-flat matroid quotients

Let `mathcal M` be a finite matroid of rank `R`, let `L(mathcal M)` be its
flat join-semilattice, and put

```math
F(X)=R-r(X).
```

Fix a flat `W` of rank `d` and take `b=W`.  Then

```math
P_bX=X\vee W,                                    \tag{BOQ.20}
```

and `W vee W=W`.  Equip the flat lattice with the explicit rank-increment
metric

```math
d_\vee(X,Y)=\max\{
 r(X\vee Y)-r(X),\ r(X\vee Y)-r(Y)
\}.                                             \tag{BOQ.21}
```

Submodularity shows that common join is nonexpansive for `d_vee`: the
marginal rank gained by adjoining `Y` can only decrease after the base is
enlarged.  The same diminishing-returns inequality proves the triangle
inequality, so (BOQ.21) is an intrinsic metric of the flat geometry rather
than a response table supplied as an assumption.

Indeed, with
`delta(X,Y)=r(X vee Y)-r(X)`, submodularity gives

```math
\delta(X\vee U,Y\vee U)\le\delta(X,Y),
\qquad
\delta(X,Z)\le\delta(X,Y)+\delta(Y,Z).          \tag{BOQ.21a}
```

The first inequality is the diminishing-returns form of submodularity.  For
the second, submodularity applied to `X vee Y` and `Y vee Z` gives

```math
r(X\vee Y)+r(Y\vee Z)
\ge r(Y)+r(X\vee Y\vee Z)
\ge r(Y)+r(X\vee Z),                            \tag{BOQ.21b}
```

which rearranges to the displayed directed triangle inequality.  Taking the
maximum of the two directed quantities proves the triangle inequality and
common-join contraction for `d_vee`; distinct flats have positive distance,
so it is a metric.  Moreover

```math
|F(X)-F(Y)|=|r(X)-r(Y)|\le d_\vee(X,Y),          \tag{BOQ.21c}
```

so the observable is one-Lipschitz.  Finally,

```math
d_\vee(X,X\vee W)=r(X\vee W)-r(X)\le d,         \tag{BOQ.22}
```

and, exactly,

```math
F(X)-F(P_bX)=r(X\vee W)-r(X).                   \tag{BOQ.23}
```

Thus BOQ.2 yields the fixed-flat quotient with arbitrary-depth error at
most `d/2`.  Its image is the interval of flats containing `W`, or
equivalently the flat lattice of the contraction `mathcal M/W`.  The bound
is sharp: the bottom flat and `W` have the same saturation and their
residual ranks differ by `d` (with the usual harmless replacement of the
empty set by the loop flat).  For the projective geometry matroid this is
exactly PMQ.1 under the identification `(X+W)/W`.

This application is not a coding relabeling.  Its geometric repair cost is
matroid rank, whereas the tropical application uses the diameter of a
subgroup in a word metric.  The common theorem predicts both because each
model has a bounded idempotent saturation invisible to all later
composition.

## 5. A no-go from the same law: convex uncertainty

Let `M` be the compact convex subsets of a finite-dimensional normed vector
space under Minkowski addition, with Hausdorff metric.  Translation by a
common summand is nonexpansive.  For a compact convex `K` containing zero,
take `b=K`.  Then

```math
b^{\star m}=mK,
\qquad
d_H(mK,K)=(m-1)\max_{x\in K}\|x\|.             \tag{BOQ.24}
```

Unless `K={0}`, the power-orbit radius is infinite.  BOQ.1 therefore proves
that a fixed Minkowski uncertainty blur has a linearly growing, sharply
attained algebra defect and cannot define an arbitrary-depth chart at fixed
error.  This recovers Proposition TDS.2b and shows that "nonexpansive
composition" alone is not the mechanism: the discarded element must have a
bounded algebraic power orbit, and exact quotient charts require it to be
idempotent.

## 6. Theory judgment

The surviving object is not an arbitrary low-oscillation congruence.  It is
a **bounded-repair algebraic saturation**:

1. the power orbit of the discarded element measures the exact cost of
   reusing the chart;
2. idempotence makes the chart a true quotient at all depths; and
3. external model geometry bounds the response width of every quotient
   fiber through one canonical representative.

This gives one law covering a tropical code profile and a matroid response
quotient, while correctly rejecting the analogous convex blur.  It is a
strict specialization of the universal CSC.1 characterization, but is more
generative: to find a reusable extremal state one may now search for a
central idempotent saturation with small image complexity and bounded repair
radius, instead of attempting to inspect every future-response fiber.

The theorem does not say that every useful congruence is generated by an
idempotent element.  Prime cyclic groups in CSC.2 and zero-separating
congruences in PMQ.2 lie outside this chart class.  The next structural
question is whether useful bounded-repair retractions which are not
principal (`x -> x star b`) admit an analogous power-orbit invariant, or
whether principal saturations already exhaust the natural strict quotients
in the current examples.
