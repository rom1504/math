# Wave 26 memo: retained-deficit noise meets the budgets but not the gap

## Status

The conditional expectation identities, row-square identity, and asymptotic
budget calculation below are **proved**.  They are independently enumerated
on the exact order-6, order-8, and order-9 examples by
`tmp/retained_noise_r26.py`.  The conclusion is a scoped assessment of the
explicit independent-noise channel and its universal entropy and spectral
certificates.  It does not lower-bound the channel's actual mutual
information for quadratic minimizers, and it does not rule out
minimizer-specific shared structure.

Write

```math
r=\frac mn,
\qquad
p_2=\frac{(m)_2}{(n)_2},
\qquad
q=Q(A).
```

For each `m`-set `S`, choose deterministically an oriented child ground
`d_S^*`.  Independently flip its child spins with probability `s`, flip its
orientation with probability `u`, and fill every outside relative spin
uniformly.  Put

```math
\theta=1-2s,
\qquad
\eta=1-2u,
\qquad
a=\eta\theta^2.
```

This is the noisy-ground channel of (R21.1)--(R21.11), now evaluated against
the centered effective loss (10.792).

## Exact retained-deficit identity

Let `Q_S=Q(A[S])`.  The child calculation from (R21.1) gives

```math
\mathbb E[c_A(S,D)\mid S]=aQ_S.
```

Uniform independent outside spins kill every cross and outside--outside
term.  Therefore the *full parent energy* has the same conditional mean:

```math
\boxed{
\mathbb E[\langle A,D\rangle\mid S]=aQ_S.
}
\tag{R26.1}
```

It follows exactly that

```math
\begin{aligned}
\mathbb E[\ell(S,D)\mid S]&=(1-a)Q_S,\\
\mathbb E[\Delta_A(D)\mid S]&=q-aQ_S,
\end{aligned}
```

and hence

```math
\boxed{
\mathbb E[\widehat\ell(S,D)\mid S]
=\bigl[1-a(1-p_2)\bigr]Q_S-r^{3/2}q.
}
\tag{R26.2}
```

For `S` uniform, with
`\overline Q=\mathbb E_{U_m}Q(A[S])`, this becomes

```math
\boxed{
\mathbb E\widehat\ell
=\overline Q-r^{3/2}q-a(1-p_2)\overline Q.
}
\tag{R26.3}
```

Thus selector-dependent child information earns exactly the last negative
term.  This is the term that was hidden when the old raw child loss alone was
used.

## Exact row-square identity

Represent the chosen child ground by a spin vector `y_S`, extended by zero
outside `S`.  Orientation does not affect row squares.  Expanding each
squared row field and averaging the independent noise gives

```math
\boxed{
\mathbb E[R_2(D)\mid S]
=n(n-1)+\theta^2
\left(\lVert A[:,S]y_S\rVert_2^2-m(n-1)\right).
}
\tag{R26.4}
```

Indeed, diagonal terms contribute `n(n-1)`; two distinct selected spin
coordinates have correlation `theta^2`; and every correlation involving a
uniform outside coordinate vanishes.  The exact minimizer spectral bound
`\lVert A\rVert_{\rm op}^2\le2q` therefore implies

```math
\boxed{
\mathbb E[R_2(D)\mid S]
\le n(n-1)+2\theta^2mq.
}
\tag{R26.5}
```

For a prescribed retained coefficient `a`, the universal row bound is best
at `u=0`, `theta^2=a`.  This is also optimal to leading order for the
universal capacity certificate (R21.11).

There is a useful comparison with noise around one fixed full parent ground
`g`.  If every full spin is independently flipped, then

```math
\boxed{
\mathbb E R_2(D)
=n(n-1)+\theta^2\bigl(R_2(g)-n(n-1)\bigr).
}
\tag{R26.6}
```

But this output is independent of `S`.  Consequently (10.796), averaged over
the noise, says

```math
\mathbb E_{S,D}\widehat\ell(S,D)
=\overline Q-r^{3/2}q
```

for every noise level.  Global noise can regularize row square, but retained
deficit makes it incapable of changing the effective restriction mean.

## The small-retention budget audit

Use no orientation noise and set

```math
a=\theta^2=n^{-1/4-c},
\qquad 0<c<1/4.
\tag{R26.7}
```

The exact capacity certificate from (R21.10)--(R21.12) is

```math
I(S;D)\le C_m^*(a)
=mf(a)+O(1),
\qquad
f(a)=\frac a2+O(a^2),
```

so, at fixed density,

```math
I(S;D)=O(n^{3/4-c}).
\tag{R26.8}
```

Since `q=O(n^{3/2})`, (R26.5) simultaneously gives

```math
\mathbb E R_2(D)
=O(n^2+an^{5/2})
=O(n^{9/4-c}).
\tag{R26.9}
```

These are exactly the information and row-square scales required by
(10.794) at `lambda` of order `n^(-3/4)`.  The old conclusion that noisy
grounds necessarily have a linear universal capacity does not survive the
new retained-deficit interface: a fixed positive retained coefficient is no
longer forced merely by coefficient matching.

However, the gain in (R26.3) at (R26.7) is only

```math
a(1-p_2)\overline Q=O(n^{5/4-c}),
\tag{R26.10}
```

which is smaller by `n^(1/4)` than the allowed residual
`O(n^(3/2-c))`.  Thus the channel closes (10.794) only if

```math
\overline Q-r^{3/2}q=O(n^{3/2-c})
\tag{R26.11}
```

already holds.  But (R26.11) is the stronger mean restriction lemma, since
some `S` has `Q(A[S])<=overline Q`; it proves the desired restriction edge
without the channel.

Conversely, a positive leading mean gap

```math
\overline Q-r^{3/2}q\ge g n^{3/2}
```

can be cancelled through (R26.3) only by a retained coefficient bounded away
from zero (with constants depending on `g` and the fixed density).  The
universal certificates then revert to `I=O(n)` and
`E R_2=O(n^(5/2))`; their contributions to (10.794) at
`lambda` of order `n^(-3/4)` are `O(n^(7/4))`, outside the desired budget.
This is not an actual-information lower bound: exceptional overlap among the
selected child grounds could still defeat the universal capacity estimate.

## Consequence for route selection

Independent weak noise is a valid budget-compatible regularizer, but its
retained selector correlation is sub-residual.  Constant noise has enough
correlation to address a leading gap, but the available universal
information and row certificates are too large.  Therefore this explicit
channel does not advance the single-cut target (10.795) or the joint-law
bound (10.794).

The surviving noise-based possibility is necessarily minimizer-specific:
prove that child ground labels have substantially more shared overlap than
the alphabet certificate detects, while simultaneously improving (R26.4),
or construct low-row-square approximate child solutions directly.  Outside
randomization by itself preserves exact child payoff, but leaves precisely
this shared-label and row-square problem.
