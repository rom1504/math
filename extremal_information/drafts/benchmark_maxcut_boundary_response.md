# Boundary responses for partial Max-Cut: exact contextual state and coarse metric entropy

**Status.** Solution-hidden benchmark report.  Sections 1--7 were derived
before any file under `extremal_information/` was opened.  Section 8 records
the subsequent collision audit.  The finite checks are in
[`verify_maxcut_boundary_response.py`](../experiments/verify_maxcut_boundary_response.py).

## 1. Experiment and conventions

Fix a labelled boundary `B` of size `w`.  A finite binary pairwise instance
`H` has boundary variables `x in X={0,1}^B`, private variables `y`, and score

```math
S_H(x,y)=c_H+
 \sum_v \theta_v(z_v)+\sum_{uv}\theta_{uv}(z_u,z_v).
```

All scores are finite in the main statements.  Hard constraints can instead
be represented by `-infinity`, together with the evident feasibility mask.
A continuation `C` has the same labelled boundary, a disjoint private
variable set, and no access to the private variables of `H`.  Gluing
identifies the two copies of `B` and adds scores.  Thus boundary factors
present on both sides are deliberately counted twice.

The partial weighted Max-Cut specialization has

```math
S_H(z)=\sum_{uv\in E(H)}a_{uv}{\bf 1}\{z_u\ne z_v\},
\qquad a_{uv}\ge0.
```

The declared observable is the optimum after an arbitrary continuation:

```math
\operatorname{Val}(H\oplus C)=\max_{x,y,y'}
  \{S_H(x,y)+S_C(x,y')\}.                         \tag{MC.1}
```

Two partial instances are **contextually equivalent** when (MC.1) is equal
for every continuation.  Equality is literal, not equality up to an
additive constant.

## 2. The response discovered by the experiment

Define the conditional optimum, or boundary response,

```math
F_H(x)=\max_y S_H(x,y),\qquad x\in X.              \tag{MC.2}
```

### Theorem MC.1 (exact gluing law)

For partial instances whose private variables are disjoint,

```math
\boxed{F_{H\oplus C}(x)=F_H(x)+F_C(x)}             \tag{MC.3}
```

pointwise.  If the common boundary is then eliminated,

```math
\boxed{\operatorname{Val}(H\oplus C)
       =\max_{x\in X}\{F_H(x)+F_C(x)\}.}          \tag{MC.4}
```

#### Proof

After `x` is fixed, `y` and `y'` occur in separate summands and can be
maximized independently.  This gives (MC.3); maximizing the result over `x`
gives (MC.4). `square`

Parallel gluing therefore uses ordinary pointwise addition.  The final
boundary elimination is the max pairing.  This differs from serial
two-boundary gluing, where eliminating a middle boundary gives max-plus
matrix multiplication.

## 3. Exposure and exact minimality

Assume first that arbitrary unary boundary rewards are valid pairwise-CSP
continuations.  For `a in X` and `M>0`, the pinning continuation has response

```math
P_{a,M}(x)=-M d_H(x,a).                             \tag{MC.5}
```

If `M>osc(F):=max F-min F`, then `a` is the unique maximizer of
`F+P_{a,M}` and

```math
\max_x\{F(x)+P_{a,M}(x)\}=F(a).                   \tag{MC.6}
```

Consequently every coordinate of (MC.2) is operationally exposed.

### Theorem MC.2 (coarsest exact state and exact response metric)

For finite binary pairwise instances,

```math
H\equiv G\quad\Longleftrightarrow\quad F_H=F_G.   \tag{MC.7}
```

Moreover, if

```math
d_{ctx}(F,G)=\sup_C
 \left|\max_x(F(x)+F_C(x))-\max_x(G(x)+F_C(x))\right|,
```

then

```math
\boxed{d_{ctx}(F,G)=\|F-G\|_\infty.}              \tag{MC.8}
```

Thus the response table, up to injective recoding, is the coarsest exact
deterministic state for this continuation experiment.

#### Proof

For every continuation response `Q`, the elementary maximum inequality gives

```math
|\max(F+Q)-\max(G+Q)|\le\|F-G\|_\infty.
```

Choose `a` at which the right side is attained and take `M` larger than both
response oscillations.  The same pin (MC.5) makes `a` optimize both tables,
so the response difference is `|F(a)-G(a)|`.  This proves (MC.8), and (MC.7)
follows.  Any exact summary assigning the same state to `H,G` would give the
same answer in every continuation, hence must separate precisely the classes
in (MC.7). `square`

### Pure Max-Cut exposure without fields

A cut response has the global-flip symmetry

```math
F_H(x)=F_H({\bf1}-x).                               \tag{MC.9}
```

Hence its natural domain is the set `\bar X=X/<x~1-x>`, of size `2^(w-1)`
for `w>=1`.  This is not a loss of distinguishability: a positive-weight
Max-Cut continuation can expose every orbit.

Fix a target `a`.  Add private vertices `r,s`, an edge `rs` of weight `A`,
and, for each `b in B`, an edge of weight `M` from `b` to `r` when `a_b=1`
and from `b` to `s` when `a_b=0`.  If `A>Mw`, every optimum has `r` and `s`
opposite, and the continuation response is

```math
Q_{a,M}(x)=A+M\bigl(w-d_{orb}(x,a)\bigr),
\quad
d_{orb}(x,a)=\min\{d_H(x,a),d_H(x,{\bf1}-a)\}.     \tag{MC.10}
```

Indeed, opposite `r,s` earn `A` and the better orientation satisfies
`w-d_orb(x,a)` terminal edges.  Equal `r,s` earn at most `Mw<A`.  If
`M>osc(F)`, (MC.10) exposes exactly the target orbit.  Therefore Theorem
MC.2 holds for pure weighted Max-Cut after replacing `X` by `\bar X`.

## 4. Additive normalization is part of the state

Adding an internal component independent of `B` translates every response
entry by its optimum.  Such a translation changes every future optimum by
the same amount.  It is therefore **not** contextual equivalence under the
literal definition (MC.1).

A fixed-anchor gauge is the simplest exact representation.  Choose
`x_* in X` (an orbit representative for Max-Cut) and store

```math
\alpha_H=F_H(x_*),\qquad f_H=F_H-\alpha_H,
\qquad f_H(x_*)=0.                                 \tag{MC.11}
```

Then gluing is simply

```math
(\alpha,f)\oplus(\beta,g)=(\alpha+\beta,f+g),      \tag{MC.12}
```

and closing the boundary returns

```math
\alpha+\beta+\max_x(f(x)+g(x)).                   \tag{MC.13}
```

Max normalization is often numerically preferable.  Write

```math
\alpha=\max F,\qquad f=F-\alpha,\qquad\max f=0.
```

For two normalized states put `u=f+g` and `m=max u`.  The exact normalized
gluing law is

```math
\boxed{(\alpha,f)\otimes(\beta,g)
 =(\alpha+\beta+m,\ u-m).}                         \tag{MC.14}
```

Discarding `alpha` retains only the response **shape** and changes the
equivalence relation to equality modulo a constant.

The operational metric on such shapes has a closed form.  Define

```math
d_{sh}([F],[G])
=\inf_{c\in\mathbb R}\sup_C
 \left|\operatorname{Val}(F,C)-
       \operatorname{Val}(G,C)-c\right|.           \tag{MC.15}
```

Then

```math
\boxed{d_{sh}([F],[G])
 =\inf_c\|F-G-c{\bf1}\|_\infty
 ={1\over2}\operatorname{osc}(F-G).}              \tag{MC.16}
```

To prove it, put `D=F-G`.  Every contextual difference lies in
`[min D,max D]`, while the pins expose every value `D(x)`.  The best single
calibration is the midpoint of this interval.  Formula (MC.16) is the exact
reason that a normalized table must not silently be compared in an
uncalibrated sup norm.

## 5. Sharp generic covering and packing modulo constants

Without a score bound there is no finite cover at any fixed additive error:
translations alone are unbounded.  Even after quotienting translations, a
scale bound is necessary because response oscillations can be scaled
arbitrarily.

Let `Y` be an `n`-point boundary state space and let

```math
\mathcal R_{n,W}=
 \{[f]:f:Y\to\mathbb R,\ \operatorname{osc}(f)\le W\}
```

with metric (MC.16).  Write `Cov_delta` for the minimum size of an internal
radius-`delta` cover and `Pack_delta` for the maximum size of a set whose
distinct members have distance at least `delta`.

### Theorem MC.3 (response-cube entropy)

For `0<delta<=W/2`,

```math
\boxed{
\left(\left\lfloor{W\over2\delta}\right\rfloor+1\right)^{n-1}
\le \operatorname{Pack}_\delta(\mathcal R_{n,W}),} \tag{MC.17}
```

and

```math
\boxed{
\operatorname{Cov}_\delta(\mathcal R_{n,W})
\le n\left(\left\lfloor{W\over\delta}\right\rfloor+1\right)^{n-1}.} \tag{MC.18}
```

In addition, when `0<delta<=W/6`,

```math
\left(\left\lfloor{W\over6\delta}\right\rfloor+1\right)^{n-1}
\le \operatorname{Cov}_\delta(\mathcal R_{n,W}),
```

and

```math
\operatorname{Pack}_\delta(\mathcal R_{n,W})
\le n\left(\left\lfloor{3W\over\delta}\right\rfloor+1\right)^{n-1}.
```

In particular, for `delta=epsilon W`,

```math
\log \operatorname{Cov}_{\epsilon W},\quad
\log \operatorname{Pack}_{\epsilon W}
=\Theta\bigl(n\log(1/\epsilon)\bigr)               \tag{MC.19}
```

up to universal constant changes in radius, throughout the nontrivial
small-`epsilon` regime.

#### Proof

For the cover, represent each class by `min f=0`, record the least minimizer,
and round every other coordinate down to the grid
`{0,delta,...,floor(W/delta)delta}`.  The rounded table still has oscillation
at most `W`, and its sup error is less than `delta`; hence its shape error is
less than `delta`.  This gives (MC.18).

For the packing, fix one coordinate `y_0` to zero and let all other
coordinates range independently over

```math
\{0,2\delta,\ldots,
 2\delta\lfloor W/(2\delta)\rfloor\}.
```

For two distinct tables, their difference is zero at `y_0` and has absolute
value at least `2delta` somewhere else.  Its oscillation is therefore at
least `2delta`, so (MC.16) gives distance at least `delta`.

For the additional covering lower bound, use the same anchored construction
with grid spacing `6delta`.  Its shape distance is at least `3delta`, so a
radius-`delta` ball contains at most one packing point.  For the packing upper
bound, apply (MC.18) at radius `delta/3`: each such cover ball has diameter at
most `2delta/3` and hence contains at most one member of a `delta`-packing.
This also justifies (MC.19) for each quantity separately, rather than only up
to an informal packing/covering comparison. `square`

For binary pairwise CSP, `n=2^w`.  For pure Max-Cut, symmetry gives the upper
bound with `n=2^(w-1)`; (MC.17) is not asserted for the restricted class of
graph-realizable cut responses.

The packing lower bound is, however, attained inside unrestricted binary
pairwise CSP rather than merely in an abstract table cube.

### Lemma MC.4 (every finite table is a pairwise-CSP response)

For every `f:{0,1}^w->[0,W]` there is a finite binary pairwise instance with
boundary response exactly `f`.

#### Proof

Introduce a selector bit `z_a` for every `a in {0,1}^w`.  Choose
`L>W`, `C>Lw`, and `P>2(C+W)`, and use the score

```math
-C+\sum_a z_a\bigl(C+f(a)-L d_H(x,a)\bigr)
 -P\sum_{a<b}z_az_b.                               \tag{MC.20}
```

Every term is unary or pairwise: `z_a d_H(x,a)` is a sum of pair potentials
between `z_a` and the boundary bits.  No active selector scores `-C<0`.
Two or more active selectors have score below zero by the choice of `P`.
One active selector `a` scores `f(a)-L d_H(x,a)`; at boundary value `x`, the
unique nonnegative maximum is obtained by `a=x` and equals `f(x)`. `square`

Taking `W=w` in (MC.17)--(MC.19) gives the requested additive
`epsilon*w` scale.  One concrete justification is a unit boundary-sensitivity
normalization: if the total variation of factors incident to each boundary
variable is at most one, then changing a set `S` of boundary bits changes
the response by at most `|S|`, and hence `osc(F)<=w`.  For Max-Cut it is
enough that every boundary vertex have weighted degree at most one.

The offset in (MC.11) must either be stored exactly or separately quantized.
If it is known to lie in an interval of length `A`, quantizing it to error
`eta` multiplies the cover size by at most `ceil(A/eta)+1`; without such a
bound no finite full-response cover exists.

## 6. A Hamming-Lipschitz refinement at error `epsilon*w`

The generic cube bound treats nearby boundary assignments as unrelated.  A
locally normalized instance has more geometry.  Let `Lip_w` be the set of
one-Lipschitz functions on the Hamming cube `{0,1}^w`, modulo constants, with
metric (MC.16).  Put

```math
V(w,r)=\sum_{j=0}^r {w\choose j}.
```

### Theorem MC.5 (coarse Lipschitz response entropy)

Let `0<eta<delta<=w`, and let `r` be a nonnegative integer with
`eta+2r<=delta`.  There is a radius-`delta` internal cover of `Lip_w` with at
most

```math
\left({2w\over\eta}+3\right)^{s+1},
\qquad
s=\left\lceil {2^w\over V(w,r)}(w\log2+1)\right\rceil.          \tag{MC.21}
```

Conversely, put `h=ceil(delta)` and assume `2h<=w`.  There is a
`delta`-packing of cardinality at least

```math
2^m,
\qquad
m=\left\lfloor
 {2^w-V(w,h-1)\over V(w,2h-1)}
 \right\rfloor.                                   \tag{MC.22}
```

For fixed `0<epsilon<1/4` and `delta=epsilon w`, these imply

```math
\log_2\operatorname{Cov}_{\epsilon w}(Lip_w)
\le 2^{(1-H_2(\epsilon/2)+o(1))w},                \tag{MC.23}
```

```math
\log_2\operatorname{Pack}_{\epsilon w}(Lip_w)
\ge 2^{(1-H_2(2\epsilon)+o(1))w}.                 \tag{MC.24}
```

The constants are not claimed sharp.  The point is that a linear additive
error reduces the logarithmic entropy from the generic `Theta(2^w)` scale to
a Hamming covering-code scale, while it remains exponential in `w`.

#### Proof of the cover

A random-centre argument gives an `r`-cover `S` of the cube with `|S|<=s`:
`k` uniform random centres leave expected uncovered size at most
`2^w exp(-kV(w,r)/2^w)`, which is below one for the displayed `k`.
Adjoin a fixed anchor `x_0`.

Normalize `f(x_0)=0` and round every sampled value upward at mesh `eta`:

```math
q_s=\eta\lceil f(s)/\eta\rceil.
```

There are at most `2w/eta+3` choices per sampled value.  Define the
McShane-type extension

```math
g(x)=\min_{s\in S\cup\{x_0\}}\{q_s+d_H(x,s)\}.    \tag{MC.25}
```

It is one-Lipschitz.  Since `q_s>=f(s)`, every term in (MC.25) is at least
`f(x)`.  If `s` is within `r` of `x`, then

```math
g(x)<f(s)+eta+r<=f(x)+eta+2r<=f(x)+delta.
```

Thus `||f-g||_infinity<=delta`.  The anchor term and Lipschitz consistency
also give `g(x_0)=0`, so the centres lie in `Lip_w`.  This proves (MC.21).

#### Proof of the packing

Fix `x_0`.  Greedily choose a set `C` outside the radius-`h-1` ball of `x_0`
whose mutual distances are at least `2h`.  Each choice deletes at most
`V(w,2h-1)` candidates, giving `|C|>=m`.

For every sign vector `sigma in {-1,+1}^C`, prescribe

```math
u_sigma(x_0)=0,\qquad u_sigma(c)=h\sigma_c.
```

These data are one-Lipschitz: a centre is at least `h` from the anchor and
oppositely signed centres are at least `2h` apart.  Extend them to the cube
by (MC.25), without quantization.  Distinct sign vectors differ by `2h` at
some centre and agree at `x_0`; their difference has oscillation at least
`2h`, so their shape distance is at least `h>=delta`.  This proves (MC.22).
For the covering estimate take, for example, `eta=1` and
`r=floor((epsilon w-1)/2)`; the extra `O(log w)` bits per landmark do not
alter the displayed double-exponential rate.  The packing estimate and both
entropy estimates use
`V(w,alpha w)=2^{(H_2(alpha)+o(1))w}`. `square`

Every function used in the lower packing is realizable by Lemma MC.4.  The
upper bound applies in particular to responses of pairwise instances, or of
Max-Cut instances, satisfying the unit terminal-sensitivity condition.  The
lower bound is only claimed for pairwise CSP under an **effective response**
Lipschitz promise; the selector realization need not obey a syntactic bound
on the sum of all incident factor magnitudes.

## 7. Finite enumeration checks

The accompanying dependency-free experiment performs the following exact
checks.

1. It enumerates every unweighted graph on three boundary and two private
   vertices, computes (MC.2), and verifies global-flip symmetry.
2. It tests (MC.3) against direct optimization after gluing sampled graph
   pairs.
3. For every enumerated graph and every boundary orbit, it constructs
   (MC.10) with `M>osc(F)` and verifies that exactly the desired orbit is
   exposed.
4. It checks the selector construction (MC.20) on deterministic random
   tables.
5. It enumerates the binary anchored packing used in (MC.17), evaluates
   (MC.16) directly, and checks its minimum separation.
6. It checks the rounded McShane extension and its advertised error bound on
   a small Hamming cube.

These computations are falsifiers for formulas and constants, not proofs.

## 8. Post-derivation collision audit and independent content

Only after Sections 1--7 had been fixed were repository files inspected.
The following parts collide with existing work.

- `phase2_feature_growth.md`, especially FG.1--FG.4, already proves that a
  conditional separator maximum is the coarsest endpoint response, that
  endpoint pinning makes the response metric a sup norm, that serial gluing
  is max-plus composition, and that a universal bounded kernel cube has
  sharp grid entropy.  MC.1--MC.3 are the one-boundary/parallel-gluing
  specialization of that general principle.
- `theorems.md`, Proposition 4.6 and Section 5, already emphasize that full
  Boolean pinning recovers a complete landscape and that packing/covering is
  the correct deterministic response complexity.
- `phase3_contextual_response_law.md` already states the general syntactic
  future-response pseudometric and coarsest congruence principle.

The following items were genuinely obtained independently in this benchmark
and were not found in those comparison passages.

1. the positive-edge two-reference Max-Cut gadget (MC.10), including its
   reduction to global-flip orbits;
2. the calibrated shape metric `osc(F-G)/2` and the exact normalized gluing
   formula (MC.14), which isolate the additive-offset issue;
3. the explicit pairwise selector realization (MC.20) used to transfer the
   full quotient-cube packing to actual binary pairwise instances; and
4. the Hamming-Lipschitz covering/packing theorem MC.5 at additive
   `epsilon*w` scale.

The benchmark therefore reproduces the expected separator response from the
continuation experiment alone, while the orbit exposure, gauge audit, and
coarse Lipschitz entropy are the non-colliding outputs.
