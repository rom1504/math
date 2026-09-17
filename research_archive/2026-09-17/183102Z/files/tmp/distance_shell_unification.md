# Distance shells, tropical projection, and resource-complete response classes

**Status.** Independent theorem draft for the benchmark campaign.  All
statements below are finite and deterministic.  The tropical identities are
classical in character; the substantive project-level conclusion is the
resource-completeness and optimal-distortion theorem obtained when they are
combined with a private profile compiler and a matching sensitivity bound.

## 1. The distance kernel and its projectors

Let `(X,d)` be a nonempty finite metric space and let `lambda>=0`.  All
statements extend verbatim to a finite pseudometric after quotienting its
zero-distance classes; this covers resource vectors with zero coordinates.
For a function `f:X->R`, define

```math
(P_lambda f)(x)=max_(y in X){f(y)-lambda d(x,y)},       \tag{DS.1}
```

and its lower companion

```math
(Q_lambda f)(x)=min_(y in X){f(y)+lambda d(x,y)}.       \tag{DS.2}
```

Write `Lip_lambda(X,d)` for the functions satisfying

```math
|f(x)-f(y)|<=lambda d(x,y).                            \tag{DS.3}
```

### Theorem DS.1 (distance-shell projection)

For every finite metric space and every `lambda>=0`:

1. `P_lambda f` and `Q_lambda f` are `lambda`-Lipschitz;
2. `P_lambda f>=f>=Q_lambda f`;
3. `P_lambda f=f` if and only if `f in Lip_lambda(X,d)`, and the same
   equivalence holds for `Q_lambda`;
4. `P_lambda f` is the least `lambda`-Lipschitz majorant of `f`, while
   `Q_lambda f` is the greatest `lambda`-Lipschitz minorant;
5. the projectors obey the bottleneck law

   ```math
   P_mu P_lambda=P_(min{lambda,mu}),
   \qquad Q_mu Q_lambda=Q_(min{lambda,mu}).          \tag{DS.4}
   ```

#### Proof

For fixed `y`, the function `x -> f(y)-lambda d(x,y)` is
`lambda`-Lipschitz.  A pointwise maximum of functions with the same
Lipschitz constant is again `lambda`-Lipschitz.  Taking `y=x` in (DS.1)
gives `P_lambda f>=f`.  If `f` is `lambda`-Lipschitz, then

```math
f(y)-lambda d(x,y)<=f(x),
```

so equality holds.  Conversely, equality makes `f` Lipschitz because
`P_lambda f` is Lipschitz.  If `h>=f` is `lambda`-Lipschitz, then

```math
h(x)>=h(y)-lambda d(x,y)>=f(y)-lambda d(x,y)
```

for every `y`, hence `h>=P_lambda f`.  This proves the majorant claim.
The statements for `Q_lambda` follow by sign reversal.

Finally,

```math
P_mu P_lambda f(x)
=max_(y,z){f(z)-lambda d(y,z)-mu d(x,y)}.
```

The triangle inequality bounds the two-distance cost below by
`min(lambda,mu)d(x,z)`, and choosing `y=x` or `y=z` attains the smaller
endpoint cost.  This proves (DS.4). `square`

In max-plus kernel notation put

```math
B_lambda(x,y)=-lambda d(x,y),
\qquad
(A odot B)(x,z)=max_y{A(x,y)+B(y,z)}.              \tag{DS.5}
```

Then (DS.4) is the kernel identity

```math
\boxed{B_lambda odot B_mu=B_(min{lambda,mu}).}     \tag{DS.6}
```

In particular `B_lambda` is max-plus idempotent.  In min-plus convention,
`D_lambda=lambda d` obeys the identical law with `min_y` in place of
`max_y`.

## 2. Isometric twists and holonomy

Let `g` be a bijective self-isometry of `(X,d)` and define the oriented shell

```math
B_(lambda,g)(a,t)=-lambda d(t,g(a)).                \tag{DS.7}
```

It acts on a profile by

```math
(T_(lambda,g)f)(t)=max_a{f(a)+B_(lambda,g)(a,t)}.  \tag{DS.8}
```

### Theorem DS.2 (bottleneck-holonomy semigroup)

For nonnegative strengths and self-isometries,

```math
\boxed{
B_(lambda,g) odot B_(mu,h)
=B_(min{lambda,mu},h circ g).}                     \tag{DS.9}
```

Consequently a chain has

```math
T_(lambda_T,g_T) cdots T_(lambda_1,g_1)f
=T_(lambda_*,G)f,
\qquad
lambda_*=min_i lambda_i,
\quad G=g_T circ cdots circ g_1.                  \tag{DS.10}
```

Equivalently,

```math
T_(lambda_*,G)f
=P_(lambda_*) (f circ G^(-1)).                    \tag{DS.11}
```

If `f` is `lambda_*`-Lipschitz, the chain therefore transmits it exactly,
up to the endpoint relabelling:

```math
T_(lambda_T,g_T) cdots T_(lambda_1,g_1)f(t)
=f(G^(-1)t).                                      \tag{DS.12}
```

#### Proof

Move `t` by `h^(-1)`.  For every intermediate `u`,

```math
lambda d(u,g(a))+mu d(u,h^(-1)t)
>=min(lambda,mu)d(g(a),h^(-1)t).
```

Choosing either endpoint attains equality.  Negating this min-plus identity
gives (DS.9), and induction proves (DS.10).  Changing variables from
`a` to `G(a)` proves (DS.11); Theorem DS.1 then gives (DS.12). `square`

This is the exact point of contact with the finite-metric directed-response
law.  The same kernel that projects arbitrary profiles onto a Lipschitz
class carries the bottleneck scalar and isometry holonomy under serial
composition.

### Corollary DS.2a (anisotropic Hamming resource lattice)

The scalar bottleneck has a coordinatewise strengthening on the Hamming
interface.  For `ell,m in R_+^w`, let

```math
d_ell(x,y)=sum_i ell_i 1{x_i ne y_i},
```

either on the full binary cube, or on the antipodal quotient with

```math
d_ell([x],[y])=min{d_ell(x,y),d_ell(x,-y)}.
```

Put `B_ell=-d_ell` and `(ell wedge m)_i=min(ell_i,m_i)`.  Then

```math
\boxed{B_ell odot B_m=B_(ell wedge m)}             \tag{DS.12a}
```

on both spaces.  Hence a chain of aligned anisotropic shells has effective
resource vector

```math
ell_*=ell^(1) wedge cdots wedge ell^(T),            \tag{DS.12b}
```

and acts once by `P_(d_(ell_*))`.  If `g,h` also permute coordinates, define
`h_*ell` by

```math
d_(h_*ell)(hx,hy)=d_ell(x,y).
```

For the twisted kernel `B_(ell,g)(a,t)=-d_ell(t,g(a))`, the exact semidirect
law is

```math
\boxed{
B_(ell,g) odot B_(m,h)
=B_(m wedge h_*ell,h circ g).}                    \tag{DS.12c}
```

Thus the earlier resource vector is transported before the coordinatewise
minimum is taken, while the permutation product is the holonomy.

For the full cube, minimization over the intermediate word separates by
coordinate: when the endpoint bits differ, that coordinate is changed in
whichever layer has smaller price.  This proves (DS.12a).  On the antipodal
quotient, expand each quotient distance as a minimum over one global choice
of endpoint orientation.  For each pair of orientation choices the same
coordinatewise calculation applies, and the remaining minimum is exactly
the antipodal distance for `ell wedge m`. `square`

Thus anisotropic Max-Cut boundary loads form an idempotent resource lattice,
not merely a one-parameter bottleneck family.

### Robust form

If an actual shell `K_i` satisfies

```math
||K_i-B_(lambda_i,g_i)||_infinity<=eta_i,
```

then max-plus composition gives

```math
\left\|K_1 odot cdots odot K_T-B_(lambda_*,G)\right\|_infinity
<=E:=sum_i eta_i.                                  \tag{DS.13}
```

Hence every transmitted profile is within `E` in uniform norm of (DS.11),
and every directed row gap is within `2E` of its ideal value.  The error is
the sum of actual kernel perturbations; there is no additional per-layer
projection loss.

## 3. An abstract compiler--bridge theorem

Consider a compositional optimization language with a labelled interface
`X`.  Its conditional response is a real function on `X`, considered either
literally or modulo an independently stored additive constant.  Assume the
following.

1. **Private compiler.**  For every target profile `f:X->R` (or, minimally,
   for every `lambda`-Lipschitz target), there is a component whose response
   at a private inner copy of `X` is `c_f+f`.
2. **Distance bridge.**  At resource `lambda`, the language realizes a
   component from the inner interface `y` to the outer interface `x` whose
   conditional response is

   ```math
   c_lambda-lambda d(x,y).                         \tag{DS.14}
   ```

   Only this bridge is charged to the exposed-interface resource.
3. **Screening.**  After the inner assignment is fixed, the compiler and
   bridge have disjoint private variables.  Eliminating the inner interface
   therefore applies the max-plus operator `P_lambda`.
4. **Sensitivity converse.**  Every admissible component of exposed
   resource at most `lambda` has a `lambda`-Lipschitz response.  A sufficient
   finite hypothesis is that for every fixed private configuration `z`,

   ```math
   |H(x,z)-H(x',z)|<=lambda d(x,x').               \tag{DS.15}
   ```

### Theorem DS.3 (resource-complete Lipschitz response class)

Under assumptions 1--4, the response shapes at exposed resource `lambda`
are exactly

```math
\boxed{
Resp_lambda/mathbb R=Lip_lambda(X,d)/mathbb R.}    \tag{DS.16}
```

Moreover distance-shell future contexts of the same resource already expose
the complete state.  If `f,g` are `lambda`-Lipschitz, then

```math
max_x{f(x)-lambda d(x,t)}=f(t),                    \tag{DS.17}
```

and consequently

```math
sup_(t in X)|max_x(f(x)-lambda d(x,t))
              -max_x(g(x)-lambda d(x,t))|
=||f-g||_infinity.                                \tag{DS.18}
```

Modulo a common calibration, the operational shape metric is

```math
d_sh([f],[g])={1\over2}osc(f-g).                   \tag{DS.19}
```

#### Proof

The sensitivity hypothesis proves the inclusion from left to right in
(DS.16).  Conversely compile a `lambda`-Lipschitz `f`, attach one distance
bridge, and eliminate the private interface.  The result is
`c_f+c_lambda+P_lambda f`, which equals `f` modulo its stored constant by
Theorem DS.1.  Equation (DS.17) is the same fixed-point identity with the
roles of `x,t` renamed.  Thus every coordinate of both tables is exposed by
a resource-admissible context, proving (DS.18); calibrating the range of
`f-g` proves (DS.19). `square`

The theorem deliberately charges only exposed-interface resource.  It says
nothing about the compiler's internal size, total internal weight, or
description length.  If those are also constrained, assumption 1 may fail
and the conclusion must be replaced by the metric entropy of the actually
compilable subclass.

## 4. Exact approximation at a resource bottleneck

The distance shell gives more than a fixed-point characterization.  It is an
optimal repair of an arbitrary response profile when additive constants are
stored separately.

Define the directed Lipschitz defect

```math
Delta_lambda(f)
=max_(x,y in X){f(y)-f(x)-lambda d(x,y)}>=0.       \tag{DS.20}
```

### Theorem DS.4 (exact resource-distortion formula)

For every finite metric profile,

```math
||P_lambda f-f||_infinity
=||f-Q_lambda f||_infinity
=Delta_lambda(f).                                 \tag{DS.21}
```

The exact distance to the resource-admissible class is

```math
\boxed{
inf_(h in Lip_lambda)d_sh([f],[h])
={Delta_lambda(f)\over2}
=d_sh([f],[P_lambda f]).}                         \tag{DS.22}
```

For literal, uncalibrated sup norm,

```math
\boxed{
inf_(h in Lip_lambda)||f-h||_infinity
={Delta_lambda(f)\over2},}                        \tag{DS.23}
```

attained by

```math
h={P_lambda f+Q_lambda f\over2}.                  \tag{DS.24}
```

Thus, in any compiler--bridge language satisfying Theorem DS.3, the one-shell
output is a nearest admissible **shape** to an arbitrary privately compiled
profile.  By shifting the compiler's free additive constant by
`-Delta_lambda(f)/2`, it is also a nearest literal response.

#### Proof

For each `x`,

```math
P_lambda f(x)-f(x)
=max_y{f(y)-f(x)-lambda d(x,y)}.
```

Maximizing over `x` gives the first equality in (DS.21); the second follows
by swapping `x,y`.  Both errors are nonnegative.  At a global maximum of
`f`, `P_lambda f=f`, so `P_lambda f-f` ranges from zero to `Delta_lambda`.
This proves the last equality in (DS.22).

For any `lambda`-Lipschitz `h`, choose `x,y` attaining (DS.20).  Since
`h(y)-h(x)<=lambda d(x,y)`,

```math
(f-h)(y)-(f-h)(x)>=Delta_lambda(f).                \tag{DS.25}
```

Hence `osc(f-h)>=Delta_lambda`, which proves the lower bound in (DS.22),
and `||f-h||_infinity>=Delta_lambda/2`, which proves the lower bound in
(DS.23).  The functions `P_lambda f` and `Q_lambda f` are Lipschitz, so their
midpoint is Lipschitz.  Since both one-sided errors lie in
`[0,Delta_lambda]`, its sup error from `f` is at most
`Delta_lambda/2`. `square`

### Corollary DS.4a (no cumulative bottleneck distortion)

For a twisted chain, compare the output with the holonomy-relabelled input.
Equations (DS.10)--(DS.11) and (DS.22) give

```math
d_sh\left(
 [T_(lambda_T,g_T)cdots T_(lambda_1,g_1)f],
 [f circ G^(-1)]
 \right)
={Delta_(lambda_*)(f)\over2}.                     \tag{DS.26}
```

Thus repeated interacting composition creates exactly the distortion of the
single weakest resource layer.  It does not pay that distortion once per
layer.  With perturbed shells, the right side increases by at most `E` from
(DS.13).

This is the strongest general corollary of the synthesis: it is simultaneously
an exact approximation theorem, a repeated-composition theorem, and an
operational response statement.

## 5. The unit-load Max-Cut shell as an instance

For the pure Max-Cut benchmark take

```math
X={+1,-1}^w/{s~-s}
```

with the anisotropic projective Hamming metric

```math
d_ell([s],[t])
=min\left\{
 sum_(i:s_i ne t_i)ell_i,
 sum_(i:s_i=t_i)ell_i
 \right\}.                                        \tag{DS.27}
```

The unrestricted pure-Max-Cut lookup gadget is the private compiler.  Join
outer spin `s_i` to inner spin `y_i` by a fresh two-edge path of weight
`ell_i`.  After eliminating its middle spin, the path contributes a constant
minus `ell_i 1{s_i ne y_i}`.  Maximizing over the common inner orientation
therefore realizes the projective distance bridge (DS.14).  Only the first
edge touches the true boundary, so its exposed load is exactly `ell_i`.

Conversely, fixing all private spins and changing the outer assignment alters
the cut score by at most the total weight incident to the changed boundary
coordinates.  This is (DS.15).  Theorem DS.3 consequently gives

```math
Resp(G_w(ell))/mathbb R
=Lip_1(X,d_ell)/mathbb R.                          \tag{DS.28}
```

For unit load, `ell_i=1`, this is the entire one-Lipschitz ball on the
projective Hamming cube.  Distance-pin Max-Cut contexts have the same unit
load and expose every profile coordinate by (DS.17).

There is also an exact resource gauge.  Put

```math
Delta_i(f)=max_s|f(s)-f(s^(i))|.                  \tag{DS.29}
```

Every realization must have boundary load at least `Delta_i(f)`.  Conversely
telescoping coordinate flips makes `f` Lipschitz for the weighted Hamming
metric with weights `Delta_i(f)`, and the shell realizes it.  Thus the
coordinatewise minimum exposed load vector is exactly `(Delta_i(f))`.

Theorem DS.4 adds a new operational interpretation.  If an arbitrary private
profile `f` must pass through outer load `ell`, its unavoidable normalized
response distortion is exactly half of its weighted-Hamming Lipschitz defect,
and the shell attains that optimum.  Scalar metric shells pay only the weakest
strength.  More strongly, aligned anisotropic Hamming shells pay the
coordinatewise minimum load vector from Corollary DS.2a, not the number of
separators.  For unrelated metrics outside this lattice no such claim is
made.

The same compiler--bridge theorem applies to universal binary CSP/Ising
interfaces whenever arbitrary private tables and metric distance factors are
available.  It need not apply to a fixed finite grammar or a polynomial-size
component class.

## 6. What is unified, and what is not

### Genuine unification

The following three facts are now one theorem rather than parallel
observations.

1. **Resource-normalized separator responses.**  A distance bridge projects
   the private separator profile onto the complete Lipschitz response class,
   and matching sensitivity proves that no other profiles occur.
2. **Query completeness.**  The same distance kernels are admissible future
   contexts and expose every coordinate without increasing the resource.
3. **Metric holonomy.**  Under serial composition, these kernels retain only
   the weakest strength and the product of their isometry labels.  An
   arbitrary Lipschitz profile is transported through the whole chain with
   no cumulative loss.

This explains why the unit-load Max-Cut shell and the permutation/metric
bottleneck calculation fit together: they are the object-completion and
serial-algebra faces of the same distance kernel.

### Classical core

The identities `P_lambda^2=P_lambda`, fixed points equal Lipschitz functions,
and distance-kernel idempotence are classical tropical/mathematical-morphology
facts.  On their own they do not constitute a new theory.  The abstract
compiler--bridge theorem is also a short consequence once its hypotheses are
stated.

The nontrivial project-level content is therefore narrower:

* the pure-Max-Cut language satisfies both the compiler and the sharp
  exposed-load sensitivity hypotheses;
* the same normalized query language remains coordinate-exposing;
* the exact distortion formula (DS.22) converts a resource deficit into an
  operational response error; and
* (DS.26) shows that this error is bottleneck-controlled, not accumulated,
  under repeated interacting composition.

### It does not collapse separator DP to the holonomy state

For an arbitrary Lipschitz profile, the separator state is still the whole
function `f`; its metric entropy can remain exponential in the separator
width.  The tiny `(lambda,G)` holonomy state describes the **restricted
shell semigroup**, not every object propagated through it.  Compression to
holonomy alone is valid only when the future query asks about the kernels or
the transported profile is already known modulo the declared relabelling.

Thus the synthesis unifies the composition operation and its resource
geometry, but it does not identify the two state-complexity regimes.  Any
claim that separator DP itself now has `O(log|Isom(X)|)` state would be
false.

## 7. Falsifiers and scope boundaries

1. If the private compiler realizes only a strict subclass, the output class
   is `P_lambda` of that subclass, not the full Lipschitz ball.
2. If a continuation may bypass the exposed interface and touch private
   compiler variables, the response profile is no longer sufficient.
3. If the bridge kernel differs from a distance cone by an uncontrolled
   amount, idempotence and the bottleneck law need not hold.
4. If the resource does not imply pointwise sensitivity such as (DS.15), the
   converse inclusion in (DS.16) can fail.
5. Internal size is uncontrolled.  The Max-Cut compiler can be exponential
   in `w`, so (DS.28) is not a polynomial-size realization theorem.
6. The isometric composition law requires a common metric and bijective
   self-isometries.  General maps or changing metrics do not obey (DS.9).
7. Exact holonomy information is visible only at a positive bottleneck scale;
   below the approximation radius, different isometries may merge.

## 8. Director judgment

The synthesis is more than a vocabulary change because Theorem DS.4 and
Corollary DS.4a formulate an exact resource--distortion law that neither the
separator table nor the bare holonomy calculation states alone.  It applies
to at least pure Max-Cut and universal pairwise-CSP/Ising interfaces.

It is not, however, a fundamentally new algebra: its engine is classical
tropical distance projection.  Its value to extremal information theory is
as a benchmark-quality **closure theorem**:

> a universal private compiler plus one resource-priced idempotent bridge
> turns a local sensitivity upper bound into an exact response-class and
> exact distortion characterization.

The next meaningful theorem would replace the exact metric bridge by a
structured approximate bridge family and characterize when its repeated
closure still has a one-bottleneck distortion law rather than linear error
accumulation.
