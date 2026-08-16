# The Hamming Grassmannian already contains classical coding rate

**Status.** The finite statements below are proved.  The exhaustive check in
[`verify_phase3_hamming_grassmannian_coding_barrier.py`](../experiments/verify_phase3_hamming_grassmannian_coding_barrier.py)
verifies the line metric, the exact ball formula, the coding sandwich, and the
systematic-chart inequality for all tested small parameters.

This note studies

```math
\Pi_{D,k}(t)=\operatorname{Pack}
\bigl(\operatorname{Gr}_k(\mathbb F_2^D),d_{\rm Hs},t\bigr),
```

where `d_Hs` is Hausdorff distance induced by Hamming distance and a packing
has pairwise distance strictly greater than the integer `t`.  For a real
threshold `Delta`, take `t=floor(Delta)`.  Thus the `P_(D,k)` notation in the
project's open-question file is `log_2 Pi_(D,k)` here.

The principal conclusion is negative but sharp: already at `k=1`, this
Grassmannian packing number agrees up to one point with the unrestricted
binary coding number.  Thus determining its asymptotic exponent in full
generality contains the classical unknown binary coding-rate problem.  At
the same time, the Hamming sphere-packing bound proves that the anticode
puncturing quotient retains exponentially too many states.

## 1. Exact line metric

Write `L_v=span(v)={0,v}` for `0 != v in F_2^D`, and put

```math
a=\operatorname{wt}(v),\qquad
b=\operatorname{wt}(w),\qquad
c=\operatorname{wt}(v+w).
```

### Lemma HG.1 (line Hausdorff formula)

For distinct nonzero `v,w`,

```math
d_{\rm Hs}(L_v,L_w)
=\max\{\min(a,c),\min(b,c)\}.                 \tag{HG.1}
```

Consequently,

```math
d_{\rm Hs}(L_v,L_w)>t
\quad\Longleftrightarrow\quad
c>t\ \text{ and }\ (a>t\ \text{or}\ b>t).   \tag{HG.2}
```

#### Proof

The directed distance from `L_v` to `L_w` is

```math
\max_{x\in\{0,v\}}d(x,L_w)
=d(v,L_w)=\min\{\operatorname{wt}(v),
                  \operatorname{wt}(v+w)\}.
```

The reverse directed distance is `min(b,c)`.  Their maximum is (HG.1),
and (HG.2) follows immediately. `\square`

The low-weight exception in (HG.2) is essential: line Hausdorff distance is
not literally Hamming distance on all nonzero representatives.

## 2. Exact reduction to unrestricted binary codes

Let `A_2(D,d)` be the largest cardinality of a binary length-`D` code of
minimum Hamming distance at least `d`; a singleton code is allowed.

### Theorem HG.2 (binary coding sandwich)

For every `D>=1` and integer `0<=t<D`,

```math
\boxed{
A_2(D,t+1)-1\le \Pi_{D,1}(t)\le A_2(D,t+1).
}                                               \tag{HG.3}
```

In particular the two quantities have exactly the same exponential rate in
every asymptotic regime in which `A_2(D,t+1)` grows exponentially.

#### Proof

For the lower bound, take an optimal binary code, translate it so that it
contains zero, and discard zero.  Every remaining word has weight greater
than `t`, and every pair has Hamming distance greater than `t`.  By (HG.2),
their one-dimensional spans form a line packing.  Over `F_2`, distinct
nonzero vectors define distinct lines.

Conversely, represent a line packing by its unique nonzero vectors `V`.
Equation (HG.2) says that all distinct members of `V` have Hamming distance
greater than `t`, and that at most one member of `V` can have weight at most
`t`.  If there is such a low-weight member, remove it and insert zero; the
result is a binary code of cardinality `|V|` and distance at least `t+1`.
If there is no low-weight member, insert zero without removing anything and
obtain a code of cardinality `|V|+1`.  Either case gives
`|V|<=A_2(D,t+1)`. `\square`

This is not merely a comparison with linear codes.  The right object is the
unrestricted, possibly nonlinear coding number.

### Corollary HG.3 (asymptotic packing barrier)

Fix `0<delta<1/2` and let `t_D=floor(delta D)`.  Then

```math
1-H_2(\delta)
\le \liminf_{D\to\infty}{1\over D}\log_2\Pi_{D,1}(t_D),
                                                        \tag{HG.4}
```

```math
\limsup_{D\to\infty}{1\over D}\log_2\Pi_{D,1}(t_D)
\le1-H_2(\delta/2).                              \tag{HG.5}
```

More precisely, the limsup and liminf are the corresponding asymptotic
rates of `A_2(D,t_D+1)`.  Hence finding their exact value is at least the
classical binary asymptotic coding-rate problem.

Puncturing `t_D` coordinates gives a response quotient with exactly
`2^(D-t_D)` possible projected subspaces of dimension zero or one.  Its state
count exceeds the maximum same-scale packing count by at least

```math
(D-t_D)-\log_2\Pi_{D,1}(t_D)
\ge
\bigl(H_2(\delta/2)-\delta-o(1)\bigr)D.          \tag{HG.6}
```

The coefficient is positive because `H_2(x)>2x` for `0<x<1/2`.

#### Proof

The greedy Gilbert bound and disjoint-ball Hamming bound are

```math
A_2(D,t+1)\ge {2^D\over\sum_{i=0}^{t}{D\choose i}},
\qquad
A_2(D,t+1)\le
{2^D\over\sum_{i=0}^{\lfloor t/2\rfloor}{D\choose i}}.
```

Combine them with (HG.3) and the standard binomial-ball entropy asymptotic.
For the quotient claim, if two lines have the same projection after deleting
`t` coordinates, every point of either line can be matched to a point of the
other outside those coordinates, so their Hausdorff distance is at most
`t`.  Thus projection is injective on a `>t` packing.  There are one zero
subspace and `2^(D-t)-1` lines in the quotient.  Finally use (HG.5). `\square`

This packing comparison alone does not prove that a smaller covering summary
exists; Proposition HG.4a below supplies that operational conclusion.  The
common-separated-host lower bound at `k=1` stores the nonzero words of a
linear code of minimum distance `t+1`.  The actual invariant in (HG.3) is the
optimal nonlinear code.  The anticode quotient is strictly too large.  Thus
the missing middle object is real even in the smallest channel dimension.

### Corollary HG.3a (presented-response transfer)

Give every line `L_v` any presentation cost
`pi_v:L_v->[0,p]`, and let

```math
F_v(x)=\min_{c\in L_v}
 \{d_{\rm Ham}(x,c)+\pi_v(c)\}.
```

If `R_D(s)` is the largest subfamily with pairwise uniform response distance
greater than `s`, then, for `s>p`,

```math
\Pi_{D,1}(s+p)\le R_D(s)\le\Pi_{D,1}(s-p),       \tag{HG.6a}
```

where a real threshold in `Pi` is replaced by its floor.  In particular,

```math
A_2(D,\lfloor s+p\rfloor+1)-1
\le R_D(s)
\le A_2(D,\lfloor s-p\rfloor+1).                \tag{HG.6b}
```

#### Proof

The presented-carrier capacity law gives

```math
\big|\|F_v-F_w\|_\infty-d_{\rm Hs}(L_v,L_w)\big|\le p.
```

Packing transfer gives (HG.6a), and HG.2 gives (HG.6b). `\square`

Thus a bounded presentation toll does not evade the coding barrier; it only
shifts the required integer distance by a bounded amount.

## 3. Exact line-ball volume

Put

```math
V_D(t)=\sum_{i=0}^t{D\choose i}
```

and, for `b=wt(w)`,

```math
I_D(t,b)=
\sum_{i=0}^b {b\choose i}
\sum_{j=0}^{D-b}{D-b\choose j}
\mathbf1\{i+j\le t,\ b-i+j\le t\}.             \tag{HG.7}
```

### Proposition HG.4 (line-ball formula)

The number of lines in the closed radius-`t` Hausdorff ball about `L_w` is

```math
|B_{\rm Hs}(L_w,t)|=
\begin{cases}
V_D(t),&\operatorname{wt}(w)>t,\\
2V_D(t)-I_D(t,\operatorname{wt}(w))-1,
   &\operatorname{wt}(w)\le t.
\end{cases}                                      \tag{HG.8}
```

#### Proof

By (HG.1), a nonzero `v` lies in the ball exactly when

```math
\operatorname{wt}(v+w)\le t,
```

or both `wt(v)<=t` and `wt(w)<=t`.  If `wt(w)>t`, this is the Hamming ball
about `w`, which does not contain zero and has size `V_D(t)`.  If
`wt(w)<=t`, it is the union of the Hamming balls about zero and `w`, with
zero deleted.  Their intersection size is (HG.7), obtained by recording the
number `i` of ones of `v` on `supp(w)` and the number `j` off that support.
Inclusion--exclusion gives (HG.8). `\square`

Thus even the finite ball geometry at `k=1` is ordinary Hamming ball geometry
up to the single low-weight exceptional region.

### Proposition HG.4a (a strict covering improvement over puncturing)

Let `N_(D,1)(t)` be the minimum number of closed radius-`t` Hausdorff balls
covering all binary lines.  Then

```math
N_{D,1}(t)
\le{(2^D-1)(\log(2^D-1)+1)\over V_D(t)}+1.                  \tag{HG.8a}
```

Consequently, for `t=floor(delta D)` and fixed `0<delta<1/2`,

```math
\log_2N_{D,1}(t)\le(1-H_2(\delta)+o(1))D.                  \tag{HG.8b}
```

Thus puncturing uses at least

```math
(H_2(\delta)-\delta-o(1))D                                  \tag{HG.8c}
```

more bits than an existing same-scale carrier summary.

#### Proof

Proposition HG.4 shows that every radius-`t` line ball has at least `V_D(t)`
members.  Select each of the `2^D-1` centers independently with probability
`p=(log(2^D-1)+1)/V_D(t)`; if `p>1`, use all centers.  Symmetry of the metric
implies that every line is covered by at least `V_D(t)` candidate centers,
so its probability of remaining uncovered is below `1/(2^D-1)`.  Add all
uncovered lines.  The expected total is at most the right side of (HG.8a),
hence such a cover exists.  Hamming-ball entropy proves (HG.8b), and
comparison with the `D-t` puncturing bits proves (HG.8c). `square`

Any response map that is one-Lipschitz in line Hausdorff distance inherits
this cover.  A bounded presentation toll changes the radius by only `O(1)`
and leaves the exponent unchanged.

## 4. A general systematic-chart upper bound

The coding obstruction persists for higher `k`, although no exact reduction
is claimed.  Let `A_q(L,d)` denote the largest length-`L`, alphabet-size-`q`
code of minimum symbol distance at least `d`.

### Theorem HG.5 (growing-alphabet shadow)

For all `1<=k<D` and `0<=t<D`,

```math
\boxed{
\Pi_{D,k}(t)
\le {D\choose k}\,A_{2^k}(D-k,t+1).
}                                               \tag{HG.9}
```

Consequently, when `t<=D-k`,

```math
\log_2\Pi_{D,k}(t)
\le k(D-k-t)+\log_2{D\choose k}.               \tag{HG.10}
```

The more explicit sphere-packing bound is

```math
\Pi_{D,k}(t)
\le {D\choose k}
{2^{k(D-k)}\over
 \displaystyle\sum_{i=0}^{\lfloor t/2\rfloor}
 {D-k\choose i}(2^k-1)^i}.                     \tag{HG.11}
```

#### Proof

Every `k`-subspace has a coordinate information set `I` of size `k`.
Partition a packing by one chosen information set.  In the chart belonging
to `I`, a subspace has the unique systematic generator

```math
[I_k\mid X],\qquad X\in\mathbb F_2^{k\times(D-k)}.
```

Regard the columns of `X` as a word of length `D-k` over an alphabet of size
`2^k`.  If `X,Y` differ in `s` columns, then matching
`(u,uX)` to `(u,uY)` shows, in both directed senses,

```math
d_{\rm Hs}(C_X,C_Y)
\le\max_u\operatorname{wt}(u(X-Y))\le s.       \tag{HG.12}
```

Thus a Hausdorff `>t` packing inside one chart maps injectively to a
`2^k`-ary code of distance at least `t+1`.  Summing over the at most
`binom(D,k)` charts proves (HG.9).  The Singleton and Hamming bounds for the
alphabet code give (HG.10) and (HG.11). `\square`

### Proposition HG.6 (the chart shadow saturates Singleton)

Put `L=D-k`.  If `2^k>=L` and `0<=t<=L`, there is a `2^k`-ary code of length
`L`, distance `t+1`, and size exactly

```math
(2^k)^{L-t}.                                    \tag{HG.13}
```

Hence no argument using only the column-support consequence in (HG.12) can
improve the leading puncturing exponent in the regime `k=Theta(D)`.

#### Proof

For `t=L`, take the singleton zero code.  Otherwise `L-t>=1`.
Choose `L` distinct elements of `F_(2^k)` and evaluate every polynomial of
degree less than `L-t` on them.  Two different polynomials agree at at most
`L-t-1` points, so their evaluation words differ at at least `t+1` points.
There are `(2^k)^(L-t)` polynomials, meeting Singleton. `\square`

The MDS words in HG.6 need not form a Hausdorff packing: changing the common
input `u` in (HG.12) can recouple the two graph subspaces and reduce their
actual distance.  That recoupling is exactly the information discarded by
the growing-alphabet shadow.

## 5. Exponent regimes and the remaining theorem

Let `t=delta D` and `k=k_D`.

* For `k=1`, the exact exponent is the unrestricted binary coding exponent,
  by HG.2.  It is strictly below the puncturing quotient upper bound by HG.3.
* For fixed `k`, HG.9 imports constant-alphabet coding upper bounds.  Exact
  determination remains at least as delicate as ordinary coding theory.
* For `k=o(D)` with `k->infinity`, the natural description scale is `kD`.
  A common linear host of asymptotic rate `R_lin(delta)` gives exponent
  `k(R_lin(delta)D-k)`, while puncturing gives at most
  `k((1-delta)D-k)+O(D)`.
* For `k_D/D->kappa`, `t_D/D->delta`, and
  `0<kappa<(1-delta)/2`, suppose that for every sufficiently large `D`
  there is a dimension-`s_D` binary linear code of minimum distance greater
  than `t_D`, with `s_D/D->R>kappa`.  Then the proved leading interval is

  ```math
  \kappa(R-\kappa)
  \ \le\ \liminf {\log_2\Pi_{D,k}(t)\over D^2}
  \ \le\ \limsup {\log_2\Pi_{D,k}(t)\over D^2}
  \ \le\ \kappa(1-\delta-\kappa),              \tag{HG.14}
  ```

  The elementary linear Gilbert estimate permits every
  `R<1-H_2(delta)`, provided `kappa<R`.  For completeness, the lower bound
  takes all `k_D`-subspaces of the dimension-`s_D` host.  If two differ,
  a vector in one but not the other is at distance greater than `t_D` from
  every vector of the other, because every nonzero host word has weight
  greater than `t_D`.  The Gaussian-binomial count has logarithm
  `k_D(s_D-k_D)+o(D^2)`.

At linear `k`, the growing alphabet is large enough that its coding shadow
meets Singleton exactly.  The unresolved issue is therefore not another
ordinary ball-volume estimate.  One needs an inequality that charges the
same-switch recoupling in (HG.12), or a construction showing that this loss
can be avoided by a quadratic-bit family.  HG.2 warns that any proposed
universal closed formula must specialize to the nonlinear binary coding
rate at `k=1`.

## 6. Director interpretation

The result narrows the Hamming Grassmannian target in three ways.

1. **The anticode quotient is rigorously nonminimal at one channel.**
   Proposition HG.4a constructs a same-scale line cover whose rate is
   strictly below puncturing; this is stronger than comparing packing counts.
2. **The middle invariant is classical coding entropy at one channel.**  An
   exact general answer cannot reasonably be expected without parameter
   restrictions, because it contains `A_2(D,d)`.
3. **The linear-channel regime has a different bottleneck.**  Ordinary
   growing-alphabet coding bounds saturate the quotient exponent.  Progress
   there must use the coherent same-input minimization defining Hausdorff
   distance between graph subspaces, not just column support.

This is a theorem-level negative resolution of the unrestricted version of
the target, not a solution of the linear-`k` Grassmannian exponent.
