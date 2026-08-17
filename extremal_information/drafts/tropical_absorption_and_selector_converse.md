# Tropical absorption and the robust selector converse

**Status.** Proof source for the depth-uniform continuation theorem.  Exact
finite checks are in
[`../experiments/verify_tropical_absorption.py`](../experiments/verify_tropical_absorption.py).

The result separates two perturbation models that cannot share one converse.
A fresh adversarial residual at every layer is controlled exactly by tangent
resets after endpoint gauges are removed.  Two fixed coherent continuation
families can instead remain close because exact semigroup relations absorb
the apparent local defect.

## 1. Bounded normal forms

### Theorem 1 (finite-semigroup absorption)

Let a finite alphabet generate a finite semigroup (S), and let

```math
L=\max_{s\in S}\min\{|w|:w\text{ represents }s\}.
```

Suppose (F) and (G) are two exact actions of (S) by nonexpansive maps
on the same metric space.  If corresponding generators satisfy

```math
\sup_x d(F_a x,G_a x)\le\varepsilon,
```

then every word (w), of arbitrary written length, satisfies

```math
\sup_x d(F_wx,G_wx)\le L\varepsilon.                \tag{1.1}
```

#### Proof

Replace (w) by a word (v) of length at most (L) representing the same
semigroup element.  Exact factorization gives (F_w=F_v) and (G_w=G_v).
Replace the factors of (F_v) by those of (G_v) one at a time; every
common suffix is nonexpansive, so each replacement costs at most
(\varepsilon).  `square`

This is neither a small-image reset nor an additive coboundary.  The
relation can be nonlinear, such as idempotence.

## 2. An all-finite max-plus counterexample

Use the projective coordinate (z=u_2-u_1) on
(\mathbb R^2/\mathbb R\mathbf1), with
(d_H(z,z')=|z-z'|/2).  In the column-output convention set

```math
S_0=\begin{pmatrix}0&0\\-1&0\end{pmatrix},
\qquad
S_\delta=\begin{pmatrix}0&\delta\\-1&0\end{pmatrix},
\qquad0<\delta<1.
```

Their projective maps are

```math
P_0(z)=\operatorname {clip}(z,0,1),
\qquad
P_\delta(z)=\operatorname {clip}(z,\delta,1).        \tag{2.1}
```

Both are idempotent.  Consequently

```math
\sup_zd_H(P_0^tz,P_\delta^tz)=\delta/2
\qquad(t\ge1).                                      \tag{2.2}
```

Yet their kernel difference has rectangular circulation (-\delta), so it
is not a row-plus-column endpoint gauge.  Their projective image diameters
are (1/2) and ((1-\delta)/2), unchanged by positive powers; there is no
(O(\delta))-image reset.  Thus gauge plus small **full-image** reset is not
complete for coherent fixed perturbations.

The distinction is sharp.  The equally close matrix

```math
\widehat S_\delta
=\begin{pmatrix}0&0\\-1+\delta&\delta\end{pmatrix}
```

induces (\widehat P_\delta(z)=\operatorname {clip}(z+\delta,0,1)), so from
(z=0) its error reaches order one after (\Theta(1/\delta)) iterations.
One-step distance and within-kernel rectangle defects therefore cannot
recognize coherent depth stability; the relations of the transition
semigroup can.

## 3. A sharp converse for arbitrary selector residuals

Let (V=\mathbb R^r/\mathbb R\mathbf1) with
(\|[v]\|_H=\operatorname {osc}(v)/2).  For
(\sigma:[r]\to[r]), let

```math
(P_\sigma v)_j=v_{\sigma(j)}.
```

A selector product is a **tangent reset** when its composite map is constant,
equivalently when it is zero on (V).  Consider a factorial language of
selector words and the disturbed recursion

```math
e_t=P_t e_{t-1}+\eta_t,
\qquad e_0=0,
\qquad\|\eta_t\|_H\le\varepsilon.                  \tag{3.1}
```

### Theorem 2 (syndetic tangent-reset criterion)

If every allowed word of length (L) contains a tangent-reset factor, then

```math
\|e_T\|_H\le L\varepsilon                           \tag{3.2}
```

for all (T\ge L).  Conversely, if there is a reset-free allowed word of
length (T), one can choose residuals in (3.1) such that

```math
\|e_T\|_H
\ge\left\lfloor{T\over r(r-1)}\right\rfloor
       \varepsilon.                                 \tag{3.3}
```

Hence arbitrary-residual stability (\|e_T\|_H\le C\varepsilon) forces
every reset-free word to have length below ((C+1)r(r-1)).

#### Proof

Unroll (3.1):

```math
e_T=\sum_{s=1}^T P_T\cdots P_{s+1}\eta_s.           \tag{3.4}
```

A reset in the final length-(L) window kills every older summand; at most
(L) nonexpansively transported residuals remain.

For the converse, every suffix product of a reset-free word is a nonconstant
selector.  For each insertion time choose an ordered output pair sent to two
different input coordinates.  One of the (r(r-1)) ordered pairs occurs at
least (\lfloor T/[r(r-1)]\rfloor) times.  At each of those times put
(+\varepsilon) and (-\varepsilon) on the two distinct input coordinates
selected by that output pair.  All
chosen contributions have the same sign at the common output pair, yielding
(3.3).  `square`

If the error has endpoint-gauge form

```math
e_t=P_te_{t-1}+h_t-P_th_{t-1}+\eta_t,               \tag{3.5}
```

then (e_t-h_t) obeys (3.1).  Thus endpoint gauge plus tangent resets is
necessary and sufficient in this adversarial linearized model.  A tangent
reset is weaker than a small full state image: only transported error
directions must die.

## 4. Twisted cycle holonomy

### Theorem 3 (affine-selector cycle criterion)

Let (A(e)=P_\sigma e+b) act on (V).  Its projective iterates are bounded
for every starting point if and only if all directed cycles (C) of the
functional graph of (\sigma) have one common mean

```math
\beta_C={1\over|C|}\sum_{j\in C}b_j.                \tag{4.1}
```

Equivalently, for some (p\in\mathbb R^r) and (\beta\in\mathbb R),

```math
b=p-P_\sigma p+\beta\mathbf1.                       \tag{4.2}
```

Then

```math
A^k(e)=p+P_\sigma^k(e-p)+k\beta\mathbf1.            \tag{4.3}
```

If two cycles have different means, projective separation grows linearly at
their mean difference.

#### Proof

Equation (4.2) telescopes to (4.3).  Conversely, (4.2) requires the sum of
(b_j-\beta) around every cycle to vanish.  If all cycle means agree, define
(p_j-p_{\sigma(j)}=b_j-\beta) consistently on every cycle and then along
the trees feeding it.  Unequal cycle means produce coordinates with
different linear drifts.  `square`

The correct holonomy is twisted by the selector transport.  Ordinary sums
of untransported edge labels can therefore miss drift or report false drift.

## 5. Consequence for the theory

Depth-uniform reuse has three proven finite mechanisms:

1. endpoint or twisted gauges telescope;
2. full-image or tangent resets erase old information; and
3. exact finite-semigroup relations give bounded normal forms and absorb a
   coherent defect.

Only the first two are robust against fresh adversarial residuals.  The third
is nevertheless essential for fixed continuation kernels.  No global
converse is claimed for arbitrary switching max-plus systems; the next
finite target is a paired-selector skew-product certificate combining
twisted cycle means, tangent resets, and exact semigroup relations.

## 6. A finite decision theorem for regular selector languages

The robust reset premise itself is decidable by a finite cycle test.  Let all
finite paths of a finite directed graph be legal trajectories and label every
edge by a selector on `[r]`; initialize at `(v,emptyset)` for each permitted
start vertex, and at every vertex when all graph paths are legal.  For a reset-free path prefix ending at graph vertex
`v`, store the set `S` of composite selectors of all its nonempty suffixes.
If an edge labelled `sigma` is appended, update

```math
S'={sigma} union {rho circ sigma:rho in S}.           \tag{6.1}
```

Reject the update if `S'` contains a constant map.  The lifted graph has at
most

```math
|V_graph| 2^(r^r-r)                                 \tag{6.2}
```

states, since `[r]^[r]` has `r^r` maps and `r` are constant.

### Theorem 4 (suffix-product cycle criterion)

There are reset-free legal paths of arbitrarily large length if and only if
a reachable directed cycle exists in the lifted graph `(v,S)`.  If no such
cycle exists and `H` is the largest edge length of a lifted path, every legal word
of length `H+1` contains a tangent-reset factor.  Hence all arbitrary-
residual trajectories satisfy

```math
||e_T||_H<=(H+1)epsilon                              \tag{6.3}
```

after the final-window condition applies.  If a reachable lifted cycle
exists, pumping it and applying Theorem 2 produces residuals with unbounded
linear response error.  The cycle is therefore an explicit instability
certificate.

#### Proof

Induction on path length shows that `S` in (6.1) is exactly the set of
products of all nonempty suffixes.  A factor ending at the newest edge is a
reset precisely when its product is constant; factors ending earlier were
checked at their own update.  Thus every reset-free legal path has a unique
lift, and every lifted path projects to one.  A finite directed graph has paths of
unbounded length exactly when a reachable directed cycle exists.  In the
acyclic case its longest path gives the claimed reset gap; in the cyclic case
repeat the cycle and invoke the adversarial lower bound.  `square`

This is a finite certificate for the robust tie-free selector model.  It
does not handle a perturbation which changes the active selector cell; that
requires the paired-cell skew product in the next theorem target.

## 7. Minimal reset quotient: kernel partitions

The suffix-set lift above is correct but nonminimal. Constants form a
two-sided ideal in the full transformation semigroup, so a word has a
constant factor exactly when its whole composite is constant. If `Pi` is
the kernel partition of the current composite, appending `sigma` updates it
by the pullback `sigma^(-1)Pi`; the one-block partition is the reset sink.
Thus the exact lift has only

```math
|V_graph|(Bell(r)-1)+1
```

states. At a fixed control vertex with the full selector alphabet, this
partition component is minimal: two different partitions contain a pair
joined in exactly one, and a continuation with image that pair resets exactly
one state. A restricted language may have a smaller residual quotient, and
distinct control vertices may also merge. The suffix construction remains a
valid upper presentation but should not be used as the benchmark complexity.

For coherent switching across different selectors, see
[`paired_selector_witness_cycles.md`](paired_selector_witness_cycles.md).
