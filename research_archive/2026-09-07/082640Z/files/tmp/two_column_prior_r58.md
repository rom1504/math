# Wave 58: two-column selected priors and the bare-column switching barrier

## Status and scope

The finite LP/KKT statements and switching identities below are exact.  The
finite calculations exhaust all selectors and oriented projective cuts of the
stored exact minimizers `A_8` and `A_9`.  They use the full bare favorable
event from (10.1352), including the retained parent-deficit term; they are not
calculations on the older hard complement subcolumn.

The main conclusion is negative but scoped.  Optimizing the two-column prior
is not a new asymptotic obligation: with a constant relaxation of the row
cap, the selected-prior LP always rounds to one column with only a constant
loss in incidence.  At a strict finite row cap, two columns can improve the
answer, but exact minimality and the retained-deficit switching identities do
not force a monotone local descent.  An exact `A_8` column is a nonglobal local
minimum for both the logarithmic scalar and the correct selected-prior scalar,
even after orientation reversal is allowed and even though all one-spin
neighbors are nonempty.  This does not rule out a global block-switching or
orbit-Harnack theorem.

## 1. The two different convex hulls

For the bare event `F_t`, put

```math
u_d=U_m(F_t^d),\qquad h_d=\log(1/u_d),\qquad a_d=e^{h_d}=1/u_d,
\qquad R_d=x^{\mathsf T}A^2x.
```

If `pi` is the captured output law of a selector-independent prior, then the
prior itself and its incidence are

```math
Z(\pi)=\frac1{\sum_d\pi_da_d},
\qquad \nu_d=Z(\pi)\pi_da_d.
\tag{R58.1}
```

Thus the selected-prior problem at row cap `R_0` is the finite LP

```math
M_{\rm sel}(R_0):=Z_{\rm sel}(R_0)^{-1}
=\min_{\pi\in\Delta}\left\{\sum_d\pi_da_d:
\sum_d\pi_dR_d\le R_0\right\}.
\tag{R58.2}
```

Its exact KKT dual is

```math
\boxed{
M_{\rm sel}(R_0)
=\sup_{\gamma\ge0}
\left[\min_d\{a_d+\gamma R_d\}-\gamma R_0\right].}
\tag{R58.3}
```

Consequently an active pair lies on a supporting line in the **`(R,a)`**
plane, not the `(R,h)` plane.  If `R_-<R_0<R_+` and the constraint binds,
then

```math
\pi_-={R_+-R_0\over R_+-R_-},\qquad
\pi_+={R_0-R_-\over R_+-R_-},
\tag{R58.4}
```

and the two active columns obey

```math
a_-+\gamma R_-=a_++\gamma R_+
=\min_d(a_d+\gamma R_d).
\tag{R58.5}
```

By contrast, (10.1368) uses the hull of `(R,h)` and the cost
`E_pi h`.  A supporting pair there minimizes `h_d+lambda R_d`.  The two
costs satisfy only

```math
\boxed{
\mathbb E_\pi h_D\le
\log\mathbb E_\pi e^{h_D}=-\log Z(\pi),}
\tag{R58.6}
```

with equality only when `h_D` is constant on the captured support.  A local
descent theorem for `h+lambda R` would therefore concern the reference-free
information LP; it is not, without an additional spread bound, a KKT theorem
for the selected-prior objective.

## 2. Constant-slack rounding closes the proposed asymptotic subproblem

Let

```math
U_*(L)=\max_{d:R_d\le L}u_d.
```

For every `b>1` and `R_0>0`, one has the exact sandwich

```math
\boxed{
U_*(R_0)\le Z_{\rm sel}(R_0)
\le {b\over b-1}\,U_*(bR_0).}
\tag{R58.7}
```

Indeed, the first inequality uses a point-mass prior.  For the second, let
`pi` be any feasible captured law and `M=E_pi e^h`.  Markov gives
`pi{R<=bR_0}>=1-1/b`.  Conditional averaging on that set produces a column
with

```math
R_d\le bR_0,
\qquad
e^{h_d}\le {M\over1-1/b},
\qquad
u_d\ge (1-1/b)Z(\pi).
\tag{R58.8}
```

This is the row-only specialization of selected-prior extraction (10.1353),
restated here because it collapses the Wave-58 proposal.  In particular,

```math
Z_{\rm sel}(R_0)\ge e^{-O(H)},\quad R_0=O(n^{9/4-c})
```

already gives one column with `u_d>=e^{-O(H)}` and
`R_d=O(n^{9/4-c})`.  Optimizing or constructing the exact two-column support
can improve strict constants, but cannot supply a new exponent-scale route to
(10.795).

## 3. Exact orientation and switching calculus for the retained deficit

Write

```math
e(d)=\sigma x^{\mathsf T}Ax,\qquad
e_S(d)=\sigma x_S^{\mathsf T}A[S]x_S,
```

and retain

```math
\Delta(d)=q-e(d),\qquad
\delta_S(d)=Q(A[S])-e_S(d),
```

so the pointwise bare loss is

```math
g_d(S)=\widehat\ell(S,d)
=\delta_S(d)-p_2\Delta(d)-B,
\qquad B=(p^{3/2}-p_2)q.
\tag{R58.9}
```

Orientation reversal has the exact sum

```math
\boxed{
g_d(S)+g_{-d}(S)
=2\{Q(A[S])-p^{3/2}q\}.}
\tag{R58.10}
```

It leaves `R_d` unchanged.  Thus an orientation pair cannot convexify row;
it merely chooses between two columns on the same vertical row line.  Also,
using (R58.10) to lower-bound either fibre requires a lower tail for the
endpoint excess `Q(A[S])-p^(3/2)q`, which is the unresolved restriction
problem rather than an independent orientation theorem.

For a vertex block `U`, flip the physical spins in `U` and put

```math
W_U=\sum_{i\in U,j\notin U}a_{ij}x_ix_j,
\qquad
W_{U,S}=\sum_{i\in U\cap S,j\in S\setminus U}a_{ij}x_ix_j.
```

Direct expansion gives the exact full-bare-column increment

```math
\boxed{
g_{d^U}(S)-g_d(S)
=4\sigma\{W_{U,S}-p_2W_U\},}
\tag{R58.11}
```

and the row increment

```math
\boxed{
R(d^U)-R(d)
=-4\sum_{i\in U}x_i(A^2x)_i
+4\lVert A[:,U]x_U\rVert_2^2.}
\tag{R58.12}
```

For singletons,

```math
R(d^i)-R(d)=4\{n-1-x_i(A^2x)_i\},
\qquad
\sum_i[R(d^i)-R(d)]=4\{n(n-1)-R(d)\}.
\tag{R58.13}
```

Equations (R58.11)--(R58.13) are the bare-event counterparts of the older
hard-complement switching formulas (10.1157)--(10.1158).  They control the
pointwise threshold displacement but not the logarithm of the number of
selectors that cross the threshold.

## 4. An exact full-bare-column local barrier on `A_8`

Take the stored exact order-eight minimizer `A_8`, `m=5`, `t=0`; then
`q=20` and there are `N=56` selectors.  Consider

```math
x=(1,-1,-1,-1,-1,-1,1,1),\qquad \sigma=-1.
\tag{R58.14}
```

Exhaustive evaluation of the exact event (R58.9) gives

```math
|F^d|=32,\qquad R(d)=72,\qquad \Delta(d)=20.
\tag{R58.15}
```

The reversed orientation has the same three values.  Every one-spin neighbor
has row `64`; four have `( |F|,Delta)=(19,8)` and four have
`( |F|,Delta)=(15,32)`.  In particular all neighbors are nonempty, and the
best row-descending neighbors also improve the parent deficit from `20` to
`8`.

At the exact price `1/20`, define both relevant scalar objectives

```math
F(d)=h_d+R_d/20,\qquad G(d)=e^{h_d}+R_d/20.
\tag{R58.16}
```

For a best one-spin neighbor, their increments are

```math
F(d^i)-F(d)=\log(32/19)-2/5>0,
```

and

```math
G(d^i)-G(d)=56/19-7/4-2/5>0.
\tag{R58.17}
```

The four remaining neighbors have `19` replaced by `15` and hence larger
increments.  Orientation reversal has zero increment.  Thus `d` is a local
minimum (strict modulo its orientation mate) for both the `(h,R)` scalar and
the correct `(e^h,R)` selected-prior scalar.

It is very far from global.  A column with `( |F|,R,Delta)=(36,8,20)` has

```math
F=\log(56/36)+8/20<\log(56/32)+72/20,
```

and

```math
G=56/36+8/20<56/32+72/20.
\tag{R58.18}
```

The best fixed-orientation block switches at projective distances
`k=0,1,2,3,4` have respectively

```text
k:             0          1          2          3          4
(|F|,R,Delta): (32,72,20) (19,64,8)  (20,40,4)  (39,32,16) (36,8,20).
```

So a two-spin block already descends both scalars, but every sequential
one-spin path must first go uphill.  This is a new scoped wall beyond the
empty-neighbor examples after (10.1159): it uses the complete bare event,
retains the deficit subsidy exactly, all neighbors are nonempty, and it
blocks both scalarizations.  It disproves a claim that exact signing
minimality plus (R58.10)--(R58.13) forces a monotone orientation/one-spin
descent.  It does **not** disprove a global block choice, orbit-average
Harnack estimate, or a migration theorem that controls the whole favorable
fibre.

## 5. What strict-cap two-column gains actually look like

The full `A_8` selected-prior optimum at the strict cap `R_0=40` is supported
on

```math
(u_+,R_+,\Delta_+)=(40/56,64,40),\qquad
(u_-,R_-,\Delta_-)=(39/56,32,16),
\tag{R58.19}
```

with captured weights `pi_+=1/4`, `pi_-=3/4`.  Hence

```math
M={371\over260},\qquad Z={260\over371},\qquad
(\nu_+,\nu_-)=({91\over371},{280\over371}).
\tag{R58.20}
```

The exact selected-prior supporting price is `gamma=7/6240`, since every
column satisfies

```math
a_d+\gamma R_d\ge {287\over195},
```

with equality on (R58.19).  The best single column under the strict cap has
coverage `39/56`, so the mixture improves it, but only by the factor

```math
{260/371\over39/56}=1.006289\ldots.
```

The information cost and selected-prior surprise are, distinctly,

```math
\mathbb E h
=\tfrac14\log(7/5)+\tfrac34\log(56/39),
\qquad
-\log Z=\log(371/260),
\tag{R58.21}
```

and the latter is strictly larger.  The high-row column in this strict gain
is the opposite-sector maximizer `Delta=2q`; restricting to `Delta<=q`
collapses the `A_8,R_0=40` optimum to the single `(39/56,32,16)` column.

This last sector feature is not universal.  On `A_9,m=6`, after restricting
to `Delta<=q=24`, the strict cap `R_0=64` has an exact adjacent one-spin pair

```math
(u_-,R_-,\Delta_-)=(30/84,56,24),\qquad
(u_+,R_+,\Delta_+)=(32/84,72,24),
\tag{R58.22}
```

with captured weights `1/2,1/2`, selected incidence `80/217`, and supporting
price `gamma=7/640`.  The best sector-restricted single column under the cap
has incidence `30/84`.  This is positive finite evidence that a switching
edge can straddle the strict row cap, but no exact-minimality inequality
forces such an edge at asymptotic scale.  In the unrestricted `A_9` LP these
columns are not optimal because opposite-sector columns have larger bare
fibres.

## 6. Resulting judgment

The direct two-column construction should be retired as an independent
asymptotic target: (R58.7) turns any successful selected prior into the single
column already required by (10.795), with harmless constant changes.  If a
two-column law arises naturally, it remains a valid certificate and (R58.3)--
(R58.5) give its exact audit.

The proposed local-minimum contradiction does not follow from exact signing
minimality, orientation pairing, or the one-spin Euler sum.  The full-bare
`A_8` barrier (R58.14)--(R58.18) is the precise finite obstruction.  A viable
continuation would need a genuinely global switching/Harnack theorem or a
fibre-stable migration result; merely finding a row- and deficit-descending
neighbor response does not control its column surprise.

Reproducible audit: `tmp/two_column_prior_r58_check.py`.
