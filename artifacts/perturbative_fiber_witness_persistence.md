# Perturbative persistence of finite fiber witnesses

Status: proved extension of the exact Boolean-channel theorem in
`fixed_fiber_boolean_channel_test.md`.  The earlier theorem handles exact
fixed-generator tensor algebras.  This note treats depth-dependent and
non-Kronecker repairs that remain quantitatively close to such a channel.

## 1. Perturbed, depth-dependent lift

Let `T_r` be a symmetric signing of order `N_r`.  At depth `r`, let `G_r` be
a symmetric matrix of order `q_r>=2` with a Boolean extremal channel

```math
G_rv_r=\sigma_r\sqrt{q_r}\,v_r,
\qquad v_r\in\{+1,-1\}^{q_r},
\qquad \sigma_r\in\{+1,-1\}.                         \tag{P1}
```

Assume `N_(r+1)=N_rq_r` and, after an allowed vertex permutation and
switching of `T_(r+1)`, write

```math
T_{r+1}=T_r\otimes G_r+E_r.                           \tag{P2}
```

The error `E_r` may be depth dependent, non-Kronecker, and dense.  It also
contains any diagonal correction needed to make the right side a signing.
Thus (P2) is strictly broader than exact closure in a fixed fiber algebra.

Start with any Boolean witness `z_(r_0)` and propagate it by

```math
z_{r+1}=z_r\otimes v_r,                               \tag{P3}
```

followed by the inverse switching/permutation used in (P2).  The initial
witness may already be fully entangled across all earlier macro and fiber
coordinates.

## 2. Perturbative persistence theorem

Define its normalized absolute energy

```math
c_r={|H_{T_r}(z_r)|\over N_r^{3/2}}.
```

Kronecker multiplication and (P1) give

```math
H_{T_r\otimes G_r}(z_r\otimes v_r)
=\sigma_r q_r^{3/2}H_{T_r}(z_r).                     \tag{P4}
```

Therefore the reverse triangle inequality proves

```math
c_{r+1}\ge c_r-\eta_r,
\qquad
\eta_r={|H_{E_r}(z_r\otimes v_r)|\over N_{r+1}^{3/2}}. \tag{P5}
```

**Theorem P (summable perturbative persistence).** If

```math
\sum_{r\ge r_0}\eta_r<\infty,                        \tag{P6}
```

then every descendant satisfies

```math
{\operatorname{cap}(T_s)\over N_s^{3/2}}
\ge c_{r_0}-\sum_{r=r_0}^{s-1}\eta_r.                \tag{P7}
```

In particular, if the complete future budget in (P6) is smaller than
`c_(r_0)-C_*`, then an all-order upper construction
`M_N<=(C_*+o(1))N^(3/2)` gives a linear landing gap along the descendant
orders:

```math
\operatorname{cap}(T_s)^{2/3}-M_{N_s}^{2/3}
\ge\left[
 \left(c_{r_0}-\sum_{r\ge r_0}\eta_r\right)^{2/3}
 -C_*^{2/3}-o(1)
\right]N_s.                                           \tag{P8}
```

The proof is the telescope of (P5).  It uses the actual propagated Boolean
state, not a spectral upper bound or a separable ansatz for the finite
witness.

## 3. Checkable norm and edit budgets

### 3.1 Operator norm

For every Boolean vector `w` of length `N_(r+1)`,

```math
|H_{E_r}(w)|
={1\over2}|w^TE_rw|
\le {N_{r+1}\over2}\|E_r\|_{op}.                    \tag{P9}
```

Hence the witness-independent sufficient condition is

```math
\eta_r\le {\|E_r\|_{op}\over2\sqrt{N_{r+1}}}.       \tag{P10}
```

The theorem therefore applies whenever

```math
\sum_r {\|E_r\|_{op}\over\sqrt{N_{r+1}}}<\infty.    \tag{P11}
```

This covers a genuinely non-Kronecker repair: `E_r` need not factor, have
small rank, or belong to the old algebra.  Only its operator size relative to
the new order is used.

If `q_r>=2` and, for some fixed `delta>0`,

```math
\|E_r\|_{op}\le C N_r^{1/2-\delta},                 \tag{P12}
```

then the future loss from a prefix of order `N_(r_0)` is at most

```math
{C\over2\sqrt2}
{N_{r_0}^{-\delta}\over1-2^{-\delta}}.              \tag{P13}
```

Thus any uniform power saving below the leading `sqrt(N)` operator scale is
summable.  A merely pointwise statement `||E_r||_op=o(sqrt(N_r))` is not
enough: its normalized errors can still form a divergent series.

### 3.2 Edge edits

Suppose `T_(r+1)` differs from the tensor model on at most `t_r` unordered
edges.  On each edited edge the difference is at most two, so for every
Boolean state

```math
|H_{E_r}(w)|\le2t_r,
\qquad
\eta_r\le {2t_r\over N_{r+1}^{3/2}}.                 \tag{P14}
```

Consequently any edit repair satisfying

```math
\sum_r {t_r\over N_{r+1}^{3/2}}<\infty              \tag{P15}
```

preserves a finite obstruction up to that explicit tail.  In particular,
`t_r=O(N_(r+1)^(3/2-delta))` has a geometrically vanishing tail when the
orders grow geometrically.  This rules out sparse refills, diagonal repairs,
and other subcritical edit modifications even when their placement is
arbitrary and depth dependent.

### 3.3 Bounded correction algebras

If

```math
E_r=\sum_{a=1}^d C_{a,r}\otimes R_{a,r},             \tag{P16}
```

then

```math
\|E_r\|_{op}
\le\sum_{a=1}^d\|C_{a,r}\|_{op}\|R_{a,r}\|_{op}.   \tag{P17}
```

For bounded fiber order, fixed `d`, bounded `R_(a,r)`, and uniformly
bounded old-coordinate factors `C_(a,r)`, (P10) is `O(N_r^(-1/2))` and is
geometrically summable.  This includes depth-dependent coefficients and
generators; exact membership in one fixed algebra is unnecessary.

The trace-zero diagonal correction in the balanced Hadamard example is the
zero-error endpoint of (P16): its energy on every Boolean old state vanishes
exactly.  The exact theorem in `fixed_fiber_boolean_channel_test.md` already
covers that endpoint, so it is not reproved here.

## 4. Quantitative necessary cost of a repair

Suppose a prefix witness has `c_(r_0)>C_*`.  Any continuation that drives its
propagated energy down to `C_*` must spend at least

```math
\sum_{r\ge r_0}
 {|H_{E_r}(z_r\otimes v_r)|\over N_{r+1}^{3/2}}
\ge c_{r_0}-C_*.                                     \tag{P18}
```

The witness-independent norm consequences are

```math
\sum_{r\ge r_0}{\|E_r\|_{op}\over2\sqrt{N_{r+1}}}
\ge c_{r_0}-C_*                                      \tag{P19}
```

or, for pure edge repairs,

```math
\sum_{r\ge r_0}{2t_r\over N_{r+1}^{3/2}}
\ge c_{r_0}-C_*.                                     \tag{P20}
```

These are necessary budget inequalities, not sufficient cap controls.  They
show that a repair of the order-56 excess `0.524977...-1/2` cannot be hidden
inside a summable small perturbation.  It must make a leading cumulative
change on the explicit propagated entangled state.

## 5. What is and is not ruled out

Theorem P newly rules out the following escapes from the exact fixed-algebra
obstruction:

1. depth-dependent generators that still have a compatible Boolean extremal
   vector at each depth and whose deviation satisfies (P11);
2. non-Kronecker dense corrections with a uniform operator-norm power saving
   as in (P12);
3. arbitrary sparse or subcritical edge refills satisfying (P15); and
4. bounded correction algebras whose coefficients change with depth but obey
   (P17) at summable scale.

It does **not** rule out:

1. a micro construction with no compatible Boolean extremal channel;
2. a correction with cumulative leading budget at least (P18), even if that
   correction can somehow be controlled by a new Boolean inequality;
3. a nonlinear order-changing map not close to any tensor model under the
   quadratic form, operator, or edit budgets above; or
4. a growing state that deliberately prevents the finite witness from having
   compatible descendants.

In particular, an abstract summable defect in the desired `b_N` recurrence
does not by itself imply (P11) or (P15).  The theorem applies only after an
exact mapping from the proposed construction to (P2) and one of the explicit
perturbation budgets.  A genuinely surviving non-tensor construction must
therefore expose a leading change on the propagated entangled witnesses and
still prove a summable all-order **cap** defect.  No current construction
does both.
