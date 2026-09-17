# Wave 49C: Bellman slack-cover identity and the exact full-slack rounding wall

## Status

- **Verified theorem:** a probability law on maximal selectors has an exact
  all-state Bellman identity, `(R49C.5)`.  It produces a positive fractional
  flip margin whenever its selected blocks are bad, and unlike the earlier
  ground-face pressure it explicitly retains every positive completion-slack
  layer.
- **Verified exact obstruction:** on the exact minimizer `A8`, the two bad
  maximal four-blocks give canonical fractional margin `8`, while all three
  integral threshold realizations have minimum margin `0`; their adaptive
  witnesses occur at full slacks `0`, `32`, and `40`.
- **Verified stronger local wall:** on the four pressure edges from
  `(10.1227)`, the exact full-state fractional game has value `4`, but its
  Boolean game has value `0`.  A half ground/half slack-eight law is an exact
  dual certificate.  Thus retaining full slack does not by itself supply a
  fractional-to-Boolean theorem.
- **Research judgment:** the identity is a useful reduction, but any integral
  successor must exploit the hypothesis that *all* maximal selectors are bad.
  A rounding theorem based on an arbitrary family of bad blocks, on threshold
  rounding of selector loads, or on a universal additive integrality gap is
  already false on `A8`.

The reproducible exact checker and output are
`tmp/full_slack_crossblock_r49_check.py` and
`tmp/full_slack_crossblock_r49.out`.

## 1. Full-slack notation

Let `d` range over all oriented parent states and write

```math
E_d=2\sum_e M_{d,e},\qquad s(d)=q_n-E_d.
\tag{R49C.1}
```

Fix a selector size `m`, put

```math
q_* = \max_{|S|=m}Q(A[S]),
```

and let `F_m` be the selectors attaining `q_*`.  For `S in F_m` and an
oriented state `y` on `S`, use the Bellman completion slack from `(10.1225)`
and the child deficit

```math
s_S(y)=q_n-\max_{d:d[S]=y}E_d,
\qquad
g_S(y)=q_*-c_S(y).
\tag{R49C.2}
```

The quantity which detects a genuinely tight lift is

```math
j_S(y):=s_S(y)+g_S(y),\qquad
\delta_S:=\min_y j_S(y).
\tag{R49C.3}
```

All four quantities in `(R49C.2)--(R49C.3)` are nonnegative, and the energy
lattice makes `j_S` and `delta_S` multiples of four.  The selector `S` has a
parent-ground/child-ground pair exactly when `delta_S=0`; a bad maximal
selector has `delta_S>=4`.

The Bellman identity gives, for every full state `d`,

```math
s_S(d[S])\le s(d).
\tag{R49C.4}
```

This inequality is the positive-slack information that was absent from the
ground-face calculation in Wave 48.

## 2. Exact weighted Bellman slack-cover identity

Let `lambda` be any probability law on `F_m`, and define its edge-load vector

```math
w_e=\mathbb P_{S\sim\lambda}\{e\subset S\},
\qquad
p_e=\frac{1-w_e}{4}.
```

For a full state `d`, abbreviate

```math
\bar s(d)=\mathbb E_\lambda s_S(d[S]),
\qquad
\bar j(d)=\mathbb E_\lambda j_S(d[S]).
```

Then the fractional flip `p` has the **exact all-state margin identity**

```math
\boxed{
s(d)+4\sum_ep_eM_{d,e}
=\frac12\{q_n-q_*+s(d)-\bar s(d)+\bar j(d)\}.}
\tag{R49C.5}
```

Indeed, `c_S(d[S])=q_*+s_S(d[S])-j_S(d[S])`, so averaging over selectors
gives

```math
2\sum_e w_eM_{d,e}=q_*+\bar s(d)-\bar j(d).
```

Subtract this from `2 sum_e M_(d,e)=q_n-s(d)` and use
`4p_e=1-w_e`.  This proves `(R49C.5)` without an inequality or a discarded
slack term.

By `(R49C.4)`, `s(d)-bar s(d)>=0`.  Therefore

```math
\boxed{
\min_d\left\{s(d)+4\sum_ep_eM_{d,e}\right\}
\ge \frac12\left(q_n-q_*+\mathbb E_\lambda\delta_S\right).}
\tag{R49C.6}
```

In particular, if `lambda` is supported on bad maximal selectors, the right
side is at least `(q_n-q_*+4)/2`.  If tight decomposition fails at size `m`,
*every* maximal selector is bad, so `(R49C.6)` applies to any law on the full
maximal-selector family.  This is the cleanest positive output of the wave:
Bellman composition and weighted selector loads do force a strictly positive
margin against every state simultaneously, but at a fractional edge point.
Here `q_n>=q_*`: averaging the external terms over all completions of any
oriented child state leaves its child energy, so at least one completion has
parent energy at least that large.

More generally, if a `lambda`-mass `beta` of its selectors is bad, then the
right side is at least `(q_n-q_*+4 beta)/2`.  Thus the identity also measures
partial failure rather than only its zero layer.

## 3. Why the identity does not round: exact `A8` selector-load wall

For the exact minimizer `A8`, at `m=4`,

```math
q_8=20,\qquad q_*=12.
```

The two bad maximal selectors are

```math
S_0=\{0,3,4,5\},\qquad S_1=\{1,2,6,7\},
```

and exact profile enumeration gives

```math
\delta_{S_0}=\delta_{S_1}=4.
\tag{R49C.7}
```

Take `lambda_(S_0)=lambda_(S_1)=1/2`.  Their internal edge sets are disjoint,
so `w_e=1/2` on the twelve internal edges and `w_e=0` on the sixteen cross
edges.  The canonical point in `(R49C.5)` is consequently

```math
p_e=\begin{cases}
1/8,&e\in E(S_0)\cup E(S_1),\\
1/4,&e\notin E(S_0)\cup E(S_1).
\end{cases}
\tag{R49C.8}
```

The theorem gives margin at least `6`; exhaustive exact evaluation of all
`256` distinct oriented states gives the stronger value

```math
\boxed{\min_d\{s(d)+4p\mathbin\cdot M_d\}=8.}
\tag{R49C.9}
```

Nevertheless, `(R49C.8)` has a natural threshold realization supported on
only three integral flips:

```math
\mathbb P(F=\varnothing)=3/4,\quad
\mathbb P(F=E(K_8)\setminus(E(S_0)\cup E(S_1)))=1/8,\quad
\mathbb P(F=E(K_8))=1/8.
\tag{R49C.10}
```

Every edge then has exactly the marginal probability `(R49C.8)`.  But the
minimum integral margin is zero for each outcome in `(R49C.10)`.  The exact
witness layers are

| integral flip | witness full slack | number of witnesses |
|:--|--:|--:|
| empty | `0` | `8` |
| sixteen cross edges | `32` | `8` |
| all twenty-eight edges | `40` | `8` |

Thus each fixed state has positive expected margin, but the maximizing state
migrates through widely separated positive-slack layers after the integral
flip is revealed.  This is a full-profile obstruction, not the old omission
of positive slack.

## 4. Exact four-edge full-slack integrality gap

The failure is already visible on the four pressure edges

```math
H=\{05,34,16,27\}.
```

Define the full-state fractional and Boolean games

```math
\Gamma_{\rm frac}(H)
=\max_{p\in[0,1]^H}\min_d
\left\{s(d)+4\sum_{e\in H}p_eM_{d,e}\right\},
```

```math
\Gamma_{\rm bool}(H)
=\max_{F\subset H}\min_d
\left\{s(d)+4\sum_{e\in F}M_{d,e}\right\}.
\tag{R49C.11}
```

At `p_e=1/2`, exact enumeration gives

```math
s(d)+2\sum_{e\in H}M_{d,e}\ge4
\qquad\text{for every }d.
\tag{R49C.12}
```

For the reverse inequality, there are two actual oriented states with

```math
(s,M_H)=(0,(-1,1,1,1)),\qquad
(s,M_H)=(8,(1,-1,-1,-1)).
```

Their half--half law has expected slack `4` and zero mean on every edge of
`H`.  It is therefore an exact dual certificate for the standard minimax
formula, proving

```math
\boxed{\Gamma_{\rm frac}(H)=4.}
\tag{R49C.13}
```

Exhausting the sixteen subsets of `H` gives ten with negative minimum margin
and six with zero minimum margin.  Hence

```math
\boxed{\Gamma_{\rm bool}(H)=0.}
\tag{R49C.14}
```

There is also a short profile explanation of the upper bound.  On the ground
layer, the four `H`-patterns with exactly one negative edge all occur.  They
defeat flip sizes zero, one, and two.  On the slack-eight layer, all patterns
with three negative edges occur; they defeat flip sizes three and four.
This is precisely the old-ground/slack-eight migration, now expressed as an
exact full-state integrality gap with matching primal and dual certificates.

The full-edge version has an even larger generic gap: `p_e=1/2` gives margin
`q_n` identically, while exact minimality forces every integral signing to
have minimum margin at most zero.  Consequently a universal additive
fractional-to-Boolean theorem cannot be the missing step; it fails both in
the global half-flip game and in the four-edge signing-specific `A8` game.

## 5. Surviving exact target

Identity `(R49C.5)` reduces a hypothetical failure of tight decomposition to
a positive, state-dependent fractional certificate.  The remaining target
can be stated without suppressing any slack layer.  If every maximal
selector is bad, seek a law `lambda` on the *entire* maximal-selector family
and an integral `F` such that, for every oriented state `d`,

```math
4\sum_e(\mathbf1_F(e)-p_e)M_{d,e}
> -\frac12\{q_n-q_*+s(d)-\bar s(d)+\bar j(d)\},
\qquad p_e=(1-w_e)/4.
\tag{R49C.15}
```

Together with `(R49C.5)`, this is exactly the strict all-state margin
`(R48C.19)` and would contradict exact minimality, proving tight
decomposition.  It is falsifiable on any proposed minimizer/profile by a
finite weighted discrepancy problem.

The `A8` calculation imposes essential hypotheses on any attempt at
`(R49C.15)`:

1. it must use that **all** maximal selectors are bad, not merely an arbitrary
   jointly realizable subfamily of bad selectors;
2. it must couple the identity's state-dependent budgets, rather than replace
   them by their positive uniform lower bound;
3. it cannot use selector-load threshold rounding or claim a universal
   additive loss below four; and
4. it must control adaptive witnesses on high slack layers as well as the
   ground and first positive layer.

This leaves a narrower but genuine route: a global all-maximal-selector
anti-migration theorem.  No such theorem is proved here.  The result does not
prove tight decomposition or convergence, but it turns “use full slack” into
the exact identity `(R49C.5)` and isolates the precise integral obstruction.
