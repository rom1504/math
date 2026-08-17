# A bounded-cap exact-sign contextual metric compiler

**Status.** Rigorous task-local theorem.  This answers the weakest
bounded-cap compiler question positively for a linear-rate family of
regular-Hadamard switchings.  It does **not** reconstruct every declared
same-support response coordinate.  Instead, it embeds the resulting
contextual metric, with constant bi-Lipschitz gain, into caps of disjoint
complete-sign parents.

The construction is deliberately outside the scope of UP.1.  A query does
not make one old configuration optimal for every child.  The old optimizer
is allowed to switch with the child--query pair.

## 1. A linear Rayleigh-separated switching code

Let `n=q^2`, where `q=2^m`.  Index coordinates by
`(u,v) in F_2^m times F_2^m`, let

```math
W_{(a,b),(u,v)}=(-1)^{a\cdot u+b\cdot v},
\qquad b(u,v)=(-1)^{u\cdot v},
\qquad \mathcal H=D_bWD_b.                            \tag{BCX.0}
```

The direct Walsh transform of the self-dual bent vector `b` gives

```math
\mathcal H^2=nI,
\qquad \mathcal H\mathbf1=q\mathbf1,
\qquad \operatorname{tr}\mathcal H=0.                \tag{BCX.1}
```

Put

```math
A=\mathcal H-\operatorname{diag}(\mathcal H),
\qquad A_s=D_sAD_s\quad(s\in\{\pm1\}^n).             \tag{BCX.2}
```

On Boolean vectors,

```math
H_(A_s)(x)={1\over2}(s\odot x)^T\mathcal H(s\odot x),
\qquad Q(A_s)={1\over2}qn.                            \tag{BCX.3}
```

### Lemma BCX.1 (linear-size two-sided Rayleigh code)

There is an absolute `gamma>0` and, for every sufficiently large `n=q^2`,
a set

```math
\mathcal S\subset\{\pm1\}^n,
\qquad |\mathcal S|\ge\exp(\gamma n),                \tag{BCX.4}
```

such that, for all distinct `s,t in mathcal S`, writing `w=s odot t`,

```math
|w^T\mathcal Hw|\le {1\over4}qn.                     \tag{BCX.5}
```

#### Proof

For a uniform Boolean `W`, the Rademacher Hanson--Wright inequality and
(BCX.1) give

```math
\Pr\{|W^T\mathcal HW|>qn/4\}
\le2\exp\left[-c\min\left\{
 {q^2n^2\over16||\mathcal H||_F^2},
 {qn\over4||\mathcal H||_op}\right\}\right]
\le2e^{-c_1n},                                      \tag{BCX.6}
```

because `||mathcal H||_F=n` and `||mathcal H||_op=q`.
Choose `M=ceil(exp(gamma n))` independent uniform switches.  Every pair
product is uniform, so a union bound makes all pairs satisfy (BCX.5) with
positive probability when `2gamma<c_1`.  Condition (BCX.5) also excludes
equal and antipodal switches, since those have normalized Rayleigh
coordinate one in absolute value. `square`

This use of Hanson--Wright is only a code-selection lemma.  Everything below
is deterministic once `mathcal S` is fixed.

## 2. The exact-sign anti-pin

For a declared query `t in mathcal S`, append only `q=sqrt(n)` new spins.
Use the query-owned complete old--new block

```math
B_t=t\mathbf1_q^T\in\{\pm1\}^{n\times q}             \tag{BCX.7}
```

and put the public positive clique `C=J_q-I_q` on the new shore.  For child
`s`, the resulting complete hollow signing of order

```math
N=n+q=n+\sqrt n                                      \tag{BCX.8}
```

has energy

```math
P_(s,t)(x,y)
=H_(A_s)(x)+(t\mathbin\cdot x)(\mathbf1\mathbin\cdot y)
  +H_C(y).                                           \tag{BCX.9}
```

Every edge is owned by exactly one source: the old block by `s`, the cross
block by `t`, and the new clique publicly.  There is no child--query joint
edge coefficient.

Write

```math
F_s(t)=Q(P_(s,t)),
\qquad E_q={q\choose2}.                               \tag{BCX.10}
```

### Theorem BCX.2 (diagonal/off-diagonal bounded-cap gap)

For every `s in mathcal S`,

```math
F_s(s)={3\over2}qn+E_q.                              \tag{BCX.11}
```

For distinct `s,t in mathcal S`,

```math
F_s(t)\le {11\over8}qn+E_q.                          \tag{BCX.12}
```

In particular every parent is an exact complete signing with

```math
Q(P_(s,t))\le {3\over2}n^{3/2}+O(n)
=O(N^{3/2}),                                         \tag{BCX.13}
```

and the query indexed by a child separates it from every other child by at
least `n^(3/2)/8`.

#### Proof

First omit the clique term.  For fixed `x`, maximizing the absolute value of

```math
H_(A_s)(x)+(t\mathbin\cdot x)(\mathbf1\mathbin\cdot y)
```

over `y` uses one of the endpoints `1 dot y=+-q`.  Since the child landscape
is even, the resulting cap is exactly

```math
Q_0(s,t)=\max\{R_+(w),R_-(w)\},                     \tag{BCX.14}
```

where `w=s odot t` and

```math
R_\pm(w)=\max_{u\in\{\pm1\}^n}
 \left\{\mathord\pm{1\over2}u^T\mathcal Hu
                 +q w^Tu\right\}.                  \tag{BCX.15}
```

For the positive channel put `K_-=2qI-mathcal H`.  It is positive definite,
and

```math
K_-^{-1}={2qI+\mathcal H\over3q^2}.                 \tag{BCX.16}
```

Completing the square, using `||u||_2^2=n`, gives

```math
R_+(w)
\le qn+{1\over2}(qw)^TK_-^{-1}(qw)
=qn\left(1+{2+\rho(w)\over6}\right),               \tag{BCX.17}
```

where

```math
\rho(w)={w^T\mathcal Hw\over qn}.                    \tag{BCX.18}
```

Likewise `K_+=2qI+mathcal H` has inverse
`(2qI-mathcal H)/(3q^2)`, so

```math
R_-(w)
\le qn\left(1+{2-\rho(w)\over6}\right).             \tag{BCX.19}
```

For an off-diagonal code pair, `|rho(w)|<=1/4`, and both channels are at
most `11qn/8`.  Restoring the clique changes the cap by at most

```math
Q(C)=E_q,                                             \tag{BCX.20}
```

which proves (BCX.12).

On the diagonal `w=1`.  The child cap and the cross-term bound give
`Q_0(s,s)<=qn/2+qn=3qn/2`, while `x=s` attains equality.  At this same old
state, `y=mathbf1` also attains `H_C(y)=E_q`.  Hence the lower bound in
(BCX.11) meets the universal upper bound `Q_0+Q(C)`, proving equality.
`square`

The small appended clique is allowed to have quadratic cap in its own order:
its order is only `sqrt(n)`, so its `Theta(q^2)=Theta(n)` contribution is
lower order at the parent `n^(3/2)` scale.  No auxiliary optimizer is assumed
to equal a prescribed vector off the diagonal.

## 3. Constant-gain embedding of the same-support metric

Declare the original same-support response language

```math
G_s(r)=Q(A_s-A_r)\qquad(r\in\mathcal S).             \tag{BCX.21a}
```

Its projective contextual metric is exactly

```math
d_0(s,t)=Q(A_s-A_t).                                  \tag{BCX.21}
```

Indeed, the reverse triangle inequality gives
`|G_s(r)-G_t(r)|<=d_0(s,t)` for every `r`, while the two queries `r=s,t`
give differences `-d_0(s,t)` and `+d_0(s,t)`.  Thus half the oscillation is
exactly `d_0`.

Let the compiled contextual metric be the projective response distance

```math
d_C(s,t)={1\over2}\operatorname{osc}_{r\in\mathcal S}
              (F_s(r)-F_t(r)).                       \tag{BCX.22}
```

### Theorem BCX.3 (bounded-cap metric compilation)

For all distinct `s,t in mathcal S`,

```math
{1\over8}d_0(s,t)\le d_C(s,t)\le d_0(s,t).          \tag{BCX.23}
```

Thus a family with `Omega(n)` same-support response bits has a disjoint,
linear-order, complete-sign, bounded-cap contextual metric compiler with
constant gain and distortion at most eight.

#### Proof

At query `r=s`, (BCX.11)--(BCX.12) give

```math
F_s(s)-F_t(s)\ge qn/8.
```

At query `r=t` the reverse difference is at least the same amount.  Hence
`d_C(s,t)>=qn/8`.

Changing the child from `A_s` to `A_t` changes no other parent edge.  The
cap functional is therefore Lipschitz:

```math
|F_s(r)-F_t(r)|\le Q(A_s-A_t)=d_0(s,t)               \tag{BCX.24}
```

for every `r`, and so `d_C<=d_0`.

Finally, the triangle inequality and (BCX.3) give

```math
d_0(s,t)\le Q(A_s)+Q(A_t)=qn.                        \tag{BCX.25}
```

Therefore `qn/8>=d_0(s,t)/8`, proving the other half of (BCX.23). `square`

There is also a direct lower bound

```math
d_0(s,t)\ge{1-\rho(s\odot t)\over2}qn\ge{3\over8}qn
                                                                  \tag{BCX.26}
```

obtained by evaluating `H_(A_s)-H_(A_t)` at `x=s`.  Hence both the original
and compiled metrics have the natural `n^(3/2)` scale on this code.

## 4. The `sqrt(n)` auxiliary width is sharp for repeated-query bridges

The construction uses exactly `q=sqrt(n)` repeated query columns.  That
order is forced within this architecture, independently of any cap bound.

### Proposition BCX.4 (rank-one query-width ceiling)

Let a response family have parents

```math
P_(s,t)(x,y)=H_s(x)+(t\mathbin\cdot x)
                         (\mathbf1_m\mathbin\cdot y)+H_(C_t)(y),
                                                                  \tag{BCX.27}
```

where `x in {+-1}^n`, `y in {+-1}^m`, and `t in {+-1}^n`.  Put
`F_s(t)=Q(P_(s,t))`.  If the auxiliary block is public, so `C_t=C`, then
every pair of child response profiles satisfies

```math
{1\over2}\operatorname{osc}_t(F_s(t)-F_u(t))
\le2nm.                                                \tag{BCX.28}
```

If `C_t` may be an arbitrary query-owned complete signing, the right side
can be replaced by

```math
2nm+m(m-1).                                           \tag{BCX.29}
```

Consequently a projective contextual gap `c n^(3/2)` with fixed `c>0`
requires `m=Omega(sqrt(n))`.  In the public-block architecture of BCX.2,
its proved gap `n^(3/2)/8` could not be obtained with
`m<sqrt(n)/16`.

#### Proof

For two queries `t,r`, the two cross Hamiltonians differ pointwise by at
most

```math
|(t-r)\mathbin\cdot x|\,|\mathbf1_m\mathbin\cdot y|
\le2nm.                                                \tag{BCX.30}
```

The cap functional is one-Lipschitz in uniform Hamiltonian norm, so
`|F_s(t)-F_s(r)|<=2nm` for every child `s` when the auxiliary block is
public.  Hence, for `Delta(t)=F_s(t)-F_u(t)`,

```math
|\Delta(t)-\Delta(r)|\le4nm.
```

Taking half the oscillation proves (BCX.28).  If the auxiliary signing also
changes, its two energies differ pointwise by at most
`2 binom(m,2)=m(m-1)`.  Repeating the same argument gives (BCX.29).
Finally, both right sides are `o(n^(3/2))` when `m=o(sqrt(n))`; the explicit
constant in the public case follows directly from (BCX.28). `square`

Thus BCX is not merely linear-order in the dominant child size.  Its
`sqrt(n)` query shore is asymptotically minimal among repeated rank-one
query interfaces capable of carrying a macroscopic projective response.

## 5. What this solves and what it does not

1. **Positive bounded-cap closure at the weakest metric level.**  The
   response vector `t -> F_s(t)` is an anti-pin: its own coordinate is high
   rather than zero.  Contextual metrics retain the information because two
   response vectors have opposite signed gaps at their two coordinates.
   This is not a pointwise approximation of
   `t -> Q(A_s-A_t)`.
2. **Why UP.1 is evaded.**  The repeated query column does not force `x=t`
   for every child.  It only gives a field of magnitude `q=sqrt(n)`, exactly
   the natural zero-temperature scale.  Off-diagonal old and auxiliary
   optimizers may depend arbitrarily on `(s,t)`.
3. **Why the rank-one bridge is now affordable.**  EL.1 used `n` appended
   copies and therefore a `Theta(n^2)` lock.  Here only `sqrt(n)` copies are
   used.  The cross cap is exactly `n^(3/2)`, while the completion cost on
   the new shore is only `O(n)`.
4. **Scope.**  The construction applies to a linear-rate subfamily of one
   regular-Hadamard switching orbit.  It does not compile the non-switching
   alternating-form Gram family of Theorem 21.26, arbitrary bounded-cap
   signings, or exact/near minimizers below cap `1/2`.
5. **New remaining question.**  The metric obstruction is no longer exact
   signs, disjointness, bounded cap, or child-dependent switching by
   themselves.  It is whether a comparable two-sided Rayleigh/synchronizing
   code exists inside a family relevant to lower-cap near-minimizers, or
   whether pointwise response-coordinate preservation is genuinely needed
   for composition.
6. **Explicit versus existential data.**  The matrix `mathcal H`, every
   child attached to a supplied switch, and every parent coefficient are
   given by the explicit formulas (BCX.0), (BCX.2), and (BCX.7)--(BCX.9).
   Only the choice of the exponentially large pair-Rayleigh code
   `mathcal S` is probabilistic/existential.  Exhaustive finite search makes
   that choice uniform but not polynomial-time; no efficiency claim is
   used by the theorem.

The exact finite regression is
[`../experiments/verify_bounded_cap_contextual_metric_compiler.py`](../experiments/verify_bounded_cap_contextual_metric_compiler.py).
