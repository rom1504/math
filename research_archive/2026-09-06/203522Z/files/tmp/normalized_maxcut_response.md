# Unit-boundary-load pure Max-Cut realizes the entire Lipschitz response class

**Status.** Rigorous theorem report for the normalized-boundary-load benchmark.
The argument is independent of the earlier spectral atom packing.  It uses the
already verified unrestricted projective-table compiler only *behind* a new
unit-load interface shell.  Thus it does not pay the lookup table separately
at the exposed boundary.

## 1. Setup and normalization

Let

```math
X_w=\{\pm1\}^w/\{s\sim-s\},
\qquad
d_{\rm pr}([s],[t])=\min\{d_H(s,t),d_H(s,-t)\}.
```

For a nonnegatively weighted graph `G` with labelled boundary `B=[w]`, put

```math
h_G([s])=\max_z \operatorname{Cut}_G(s,z).
```

The syntactic load of boundary vertex `i` is

```math
\ell_i(G)=\sum_{e\ni i}c_e.
```

An edge with both endpoints in the boundary contributes to both endpoint
loads.  Let `G_w(ell)` be the class with `ell_i(G)<=ell_i` for every `i`.
Constants in a response are stored separately, since private components can
add an arbitrary common offset.  Thus the natural shape metric is

```math
d_{\rm sh}([f],[g])
=\inf_c\|f-g-c\mathbf1\|_\infty
=\tfrac12\operatorname{osc}(f-g).
```

## 2. Exact boundary-load realization theorem

### Theorem 1 (load vector equals the anisotropic Lipschitz envelope)

For every vector `ell=(ell_1,...,ell_w)>=0`, the response shapes of pure
weighted Max-Cut components in `G_w(ell)` are exactly

```math
\boxed{
\operatorname{Resp}(\mathcal G_w(\ell))/\mathbb R
=\operatorname{Lip}_1(X_w,d_\ell)/\mathbb R,}
\tag{1}
```

where

```math
d_\ell([s],[t])
=\min\left\{\sum_{i:s_i\ne t_i}\ell_i,
             \sum_{i:s_i=t_i}\ell_i\right\}.
\tag{2}
```

In particular, if every exposed boundary vertex has load at most one, the
realizable response shapes are **all** one-Lipschitz functions on the
projective Hamming cube, not a smaller model-dependent subclass.

#### Necessity

Fix two oriented boundary words `s,t` and any private spin assignment `z`.
Changing the boundary coordinates in
`D={i:s_i!=t_i}` changes the cut contribution by at most
`sum_(i in D) ell_i(G)`: compare using the same private assignment, and note
that double-counting a boundary-boundary edge only weakens this bound.
Taking maxima over `z` in both directions gives

```math
|h_G(s)-h_G(t)|\le\sum_{i\in D}\ell_i.
\tag{3}
```

The profile is invariant under global spin flip, so applying (3) to `t` and
to `-t` proves the `d_ell`-Lipschitz bound.

#### Sufficiency: the unit-load McShane shell

Let `f:X_w->R` be `d_ell`-Lipschitz.  Add a common constant so that `f>=0`.
The projective lookup theorem gives a pure weighted Max-Cut component `H_f`
with an *inner* labelled interface `Y={y_1,...,y_w}` and

```math
h_{H_f}([y])=C_f+f([y]).
\tag{4}
```

The loads of the inner interface are irrelevant because these vertices now
become private.  Introduce the true outer boundary spins `s_i`.  Between
`s_i` and `y_i`, put a fresh two-edge path

```text
s_i --(ell_i)-- p_i --(ell_i)-- y_i.
```

After maximizing `p_i`, this path scores `2 ell_i` if `s_i=y_i` and
`ell_i` otherwise.  Therefore the outer response is

```math
\begin{aligned}
h_G([s])
 &=C_f+2\sum_i\ell_i
   +\max_y\left\{f([y])-\sum_{i:s_i\ne y_i}\ell_i\right\}\\
 &=C_f+2\sum_i\ell_i+f([s]).
\end{aligned}
\tag{5}
```

The second equality is precisely the max-plus McShane identity: `y=s`
attains `f([s])`, while Lipschitzness gives

```math
f([y])-\sum_{i:s_i\ne y_i}\ell_i\le f([s]).
```

Only the first edge of the `i`-th path meets the outer boundary, so
`ell_i(G)=ell_i`.  This proves (1).  Notice the point of the construction:
the exponentially large lookup is compiled into a private landscape and is
paid only once; a single anisotropic distance kernel couples it to the
separator.

### Corollary 1 (the exact minimum boundary load of a shape)

For a projective response shape `f`, define its coordinate oscillations

```math
\Delta_i(f)=\max_s|f(s)-f(s^{(i)})|.
\tag{6}
```

Then

```math
\boxed{
\inf_{G:\,[h_G]=[f]}\ell_i(G)=\Delta_i(f)
\quad\text{coordinatewise},}
\tag{7}
```

and consequently

```math
\boxed{
\inf_{G:\,[h_G]=[f]}\sum_i\ell_i(G)
=\sum_i\Delta_i(f).}
\tag{8}
```

Indeed, (3) gives the lower bounds.  Conversely, telescoping coordinate
flips shows that `f` is Lipschitz for the weighted Hamming metric with
weights `ell_i=Delta_i(f)`, and Theorem 1 realizes it with exactly those
loads.  Thus the promise `sum_i ell_i(G)<=w` is characterized exactly by
`sum_i Delta_i(f)<=w`; it is not merely bounded by that seminorm.

This is a useful structural falsifier: normalized *syntactic boundary load
alone* cannot force a profile algebra smaller than the corresponding
Lipschitz ball.

### Abstract form: the max-plus Lipschitz projector

The mechanism is not peculiar to the lookup formula.  For any finite metric
space `(X,d)`, define

```math
(P_df)(x)=\max_y\{f(y)-d(x,y)\}.
\tag{8a}
```

Then `P_df` is one-Lipschitz, `P_df>=f`, and

```math
\boxed{P_df=f\quad\Longleftrightarrow\quad f\in\operatorname{Lip}_1(X,d).}
\tag{8b}
```

Moreover the distance kernel is max-plus idempotent:

```math
\max_y\{-d(x,y)-d(y,z)\}=-d(x,z).
\tag{8c}
```

Thus a language with (i) a private compiler for arbitrary inner profiles and
(ii) a realizable distance bridge automatically realizes its entire
Lipschitz response ball at the resource cost of one bridge.  If the same
resource implies the matching Lipschitz upper bound, the characterization is
exact.  Theorem 1 is this abstract projector instantiated in pure Max-Cut.
The same statement applies to universal binary-CSP/Ising interfaces.  It
also identifies this shell with the metric-kernel holonomy algebra appearing
in the directed-response campaign: repeated shells collapse exactly rather
than accumulating loss.

## 3. Unit-load future contexts still expose every coordinate

One might worry that (1) uses arbitrary future attachments to distinguish
its tables.  It does not.

### Theorem 2 (restricted-context isometry)

Let `f,g` be responses of components with boundary load at most one per
vertex, and let future contexts also have boundary load at most one per
vertex.  Then

```math
\boxed{
\sup_C\left|
 \max_x(f(x)+h_C(x))-\max_x(g(x)+h_C(x))
\right|=\|f-g\|_\infty.}
\tag{9}
```

For every target `t in X_w`, the standard pure-Max-Cut pin has unit exposed
load and response

```math
q_t(x)=K-d_{\rm pr}(x,t).
\tag{10}
```

Since `f` is one-Lipschitz,

```math
\max_x\{f(x)-d_{\rm pr}(x,t)\}=f(t),
\tag{11}
```

and likewise for `g`.  Hence (10) exposes `f(t)-g(t)` without amplification.
Choosing a coordinate attaining the sup norm proves the lower bound in (9);
the elementary maximum inequality proves the upper bound.

Modulo a separately stored offset, the same argument exposes the whole
range of `f-g`, so the restricted operational shape metric is exactly
`d_sh`.  Thus neither the object class nor the declared query class needs an
unnormalized boundary load.

## 4. Macroscopic response rate-distortion

Let `L_w=Lip_1(X_w,d_pr)/R`; by Theorem 1 this is exactly the unit-load
pure-Max-Cut response-shape class.  Write `Cov_delta(L_w)` for its internal
covering number in `d_sh`.  Let

```math
V(w,r)=\sum_{j=0}^{\lfloor r\rfloor}\binom wj.
```

For radii below `w/2`, a ball in `X_w` has exactly `V(w,r)` elements.

### Theorem 3 (double-exponential macroscopic response complexity)

For every fixed `0<epsilon<1/4`,

```math
\boxed{
2^{(1-H_2(2\epsilon)+o(1))w}
\le \log_2\operatorname{Cov}_{\epsilon w}(L_w)
\le 2^{(1-H_2(\epsilon)+o(1))w}.}
\tag{12}
```

Thus a deterministic state answering all normalized pure-Max-Cut futures to
additive error `epsilon w` requires a number of bits exponential in `w`;
equivalently its number of possible states is double exponential in `w`.
The construction establishing the lower bound has total exposed boundary
load exactly `w` in every member.

#### Upper bound

Let `D=floor(w/2)` and choose an `r`-net `S` of `X_w`.  Normalize
`min f=0`, so `0<=f<=D`, and round each `f(s)`, `s in S`, upward to a
multiple `q_s` of `eta`.  Define

```math
u(x)=\min_{s\in S}\{q_s+d_{\rm pr}(x,s)\},
\quad
l(x)=\max_{s\in S}\{q_s-d_{\rm pr}(x,s)\},
\quad
g=(u+l)/2.
\tag{13}
```

Both envelopes and `g` are one-Lipschitz.  A landmark within distance `r`
and the global Lipschitz inequalities give

```math
-r\le g-f<r+\eta.
```

Consequently `d_sh([f],[g])<r+eta/2`.  The number of possible centers is at
most

```math
(D/\eta+2)^{|S|}.
\tag{14}
```

Random centers give

```math
|S|\le \left\lceil{2^{w-1}\over V(w,r)}
 (1+(w-1)\log2)\right\rceil.
\tag{15}
```

Take `eta=o(w)` and `r=(epsilon-o(1))w`; (14)--(15) and Hamming-ball
asymptotics give the upper half of (12).  By Theorem 1 every center `g` is
itself a unit-load pure-Max-Cut response shape, so this is an internal cover,
not an ambient relaxation.

#### Lower bound

Put `h=(epsilon+o(1))w`, with the positive `o(1)` chosen so that `2h` is an
integer and `2h>2 epsilon w`.  A greedy code gives a set `C subset X_w` with
pairwise distance at least `2h` and

```math
|C|\ge {2^{w-1}\over V(w,2h-1)}.
\tag{16}
```

For every subset `U subset C` of cardinality `floor(|C|/2)`, prescribe on
`C` the values

```math
a_U(c)=\begin{cases}+h,&c\in U,\\-h,&c\notin U.\end{cases}
\tag{17}
```

These data are one-Lipschitz and therefore admit a one-Lipschitz McShane
extension `f_U` to all of `X_w`.  If `U!=V` have equal cardinality, there is
one point of `U\setminus V` and one of `V\setminus U`; hence `f_U-f_V`
takes both values `+2h` and `-2h` on `C`.  Therefore

```math
d_sh([f_U],[f_V])\ge2h>2 epsilon w.
\tag{18}
```

No ball of radius `epsilon w` contains two members, and so

```math
\operatorname{Cov}_{\epsilon w}(L_w)
\ge\binom{|C|}{\lfloor|C|/2\rfloor}.
\tag{19}
```

Equations (16)--(19) give the lower half of (12).  This constant-weight
labeling is important: arbitrary sign labels only force shape separation
`h`, losing an unnecessary factor two in the rate-distortion lower bound.

For a literal, rather than shape-normalized, response packing, choose the
lower McShane extensions

```math
f_U(x)=\min_{c\in C}\{a_U(c)+d_{\rm pr}(x,c)\}.
```

All lie in the same fixed interval `[-h,D+h]`.  Shift them all by `h`, use
the common-offset padding in the inner lookup compiler, and apply the unit
shell.  The resulting pure-Max-Cut components have a common literal offset,
unit load at every true boundary vertex, and the separations (18).

## 5. Interpretation and remaining boundary

This gives the requested alternative (i): a double-exponential packing at
macroscopic error with *linear total / unit per-coordinate exposed boundary
load*.  It also gives an exact description of the whole normalized response
class, not just a packing.

The result is not obtained by rescaling the earlier lookup gadget.  Directly
placing that gadget at the exposed separator costs
`4 sum_t F(t)` at every boundary coordinate.  Here the lookup interface is
private, and the distance shell performs a single joint max over its hidden
word.  Cancellation/selection occurs before the outer response is paid.

What the theorem does **not** show is that polynomial-size graphs have the
same entropy.  The private compiler remains exponential in `w`.  Therefore
the next genuinely different question is:

```math
\boxed{\text{What is the response entropy under both unit boundary load
and }|V(G)|,|E(G)|\le\operatorname{poly}(w)?}
```

The present result decisively rules out normalized boundary sensitivity by
itself as a compression mechanism.  Any smaller state theorem must use a
global resource promise (size, total internal weight/bit complexity, or a
restricted construction grammar), not merely separator load.

The companion script
`/home/math/quadra/tmp/verify_normalized_maxcut_shell.py` exhausts every
oriented outer word in 1,000 random anisotropic projective trials for
`w=2,...,6`.  It checks the Lipschitz inequalities, the exact shell identity
(5), and the realization metric generated by the coordinate oscillations.
It prints:

```text
PASS: 1000 anisotropic projective shell trials (w=2..6)
```
