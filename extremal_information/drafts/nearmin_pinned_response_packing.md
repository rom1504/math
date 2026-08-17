# A linear-rate pinned-response packing inside every near-minimizer halo

Date: 2026-08-17.

Status: proof draft.  This is a Level-5 negative theorem for the contextual
response metric with amplitude-`n` linear fields.  It does **not** use a
bounded-cap all-spins-free completion, and therefore does not settle the
balanced-interface version relevant to a low-cap parent.

## 1. Statement

For a hollow signing `a`, write

```math
H_a(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad Q(a)=\max_x|H_a(x)|,
```

and define its one-sided field response by

```math
R_a(g)=\max_{x\in\{+-1\}^n}\{H_a(x)+g\mathbin\cdot x\}.       \tag{PR.1}
```

### Theorem PR.1 (near-minimality does not compress pinned responses)

There are absolute constants `c,C,delta>0` such that the following holds.
Let

```math
{C\over\sqrt n}\le\kappa\le {1\over10}.                    \tag{PR.2}
```

For every exact order-`n` minimizer `a`, there are a set
`U subset {+-1}^n/{+-1}` and exact signings `(b^u)_(u in U)` such that

```math
|U|\ge e^{cn},
\qquad Q(b^u)\le M_n+2\kappa n^{3/2},                        \tag{PR.3}
```

and, for every distinct `u,v in U`,

```math
\boxed{
 R_(b^u)(nu)-R_(b^v)(nu)
 \ge\delta\kappa n^{3/2}.}                                  \tag{PR.4}
```

Consequently any summary which answers all queries `(nu)_(u in U)` to
additive error below `delta kappa n^(3/2)/2`, uniformly on the
`2kappa`-near-minimizer class, has at least `|U|` states and hence requires
`Omega(n)` bits.

The state family is a one-hot packing, not an independently writable
`Omega(n)`-bit cube.  Its information lower bound is nevertheless linear
because the message set itself has exponential cardinality.

## 2. A universal code of query spins

Fix `rho=1/2`.  For all sufficiently large `n`, the standard random-code
argument gives a set of projective classes
`U subset {+-1}^n/{+-1}` with

```math
|U|\ge e^{c_0n},
\qquad |u\mathbin\cdot v|\le\rho n\quad(u\ne v),              \tag{PR.5}
```

where one may take any fixed `c_0<1/32`.  Choose and fix one representative
`u in {+-1}^n` of each class; every later occurrence of `u`, `nu`, and
`b^u` uses these representatives.  Indeed, for independent uniform
Boolean vectors,

```math
\Pr\{|u\mathbin\cdot v|>n/2\}\le2e^{-n/8},
```

and a union bound over `exp(c_0n)` sampled vectors proves existence after
removing projective repetitions.

For `u in U`, put

```math
z_u=(u_iu_j)_(i<j),
\qquad g_u=a\odot z_u,
\qquad N_u=\{e:(g_u)_e=-1\}.                                  \tag{PR.6}
```

Independently for every `u`, include each edge of `N_u` in `F_u` with
probability

```math
q={\kappa\over\sqrt n},                                       \tag{PR.7}
```

and obtain `b^u` by flipping the edges in `F_u`.

## 3. The cancellation identity

For distinct `u,v`, direct expansion gives

```math
H_(b^u)(u)-H_(b^v)(u)
=2|F_u|+2\sum_(e\in F_v)(g_u)_e.                              \tag{PR.8}
```

The expectation of this random difference has no dependence on the base
signing.  Since

```math
1_(N_v)(e)={1-(g_v)_e\over2},
```

we get

```math
\begin{aligned}
\mathbb E[2|F_u|]
 &=q\left(E-H_a(u)\right),\\
\mathbb E\left[2\sum_(e\in F_v)(g_u)_e\right]
 &=q\left(H_a(u)-\langle z_u,z_v\rangle\right).
\end{aligned}                                                 \tag{PR.9}
```

Thus, with `E=binom(n,2)`,

```math
\mathbb E\left[H_(b^u)(u)-H_(b^v)(u)\right]
=q\left(E-\langle z_u,z_v\rangle\right)
={q\over2}\left(n^2-(u\mathbin\cdot v)^2\right)
\ge {3\over8}\kappa n^{3/2}.                                \tag{PR.10}
```

This identity is the point of the construction: the base landscape cancels
before an absolute value or a separate-channel estimate is taken.

## 4. Simultaneous concentration and the cap bound

The right side of (PR.8) is a sum of independent variables bounded by two,
with total variance at most

```math
4q(|N_u|+|N_v|)\le8qE\le4\kappa n^{3/2}.                     \tag{PR.11}
```

Bernstein's inequality therefore gives

```math
\Pr\left\{
 H_(b^u)(u)-H_(b^v)(u)<{3\over16}\kappa n^{3/2}
\right\}
\le\exp(-c_1\kappa n^{3/2}).                                 \tag{PR.12}
```

The same estimate, or a Chernoff bound, gives

```math
\Pr\{|F_u|>\kappa n^{3/2}\}
\le\exp(-c_2\kappa n^{3/2}),                                \tag{PR.13}
```

because `qE<=kappa n^(3/2)/2`.  If `kappa sqrt(n)>=C`, choose
`c<c_0` and then `C` large enough.  A union bound over all vertices and all
ordered pairs of an `e^{cn}`-element subcode makes (PR.12)--(PR.13) hold
simultaneously with positive probability.

Fix such a realization.  Edge Lipschitzness now gives

```math
Q(b^u)\le Q(a)+2|F_u|
        \le M_n+2\kappa n^{3/2},                              \tag{PR.14}
```

while (PR.12) gives the claimed oriented evaluation gap.

## 5. Pinning is an exact optimized query

For any two Boolean spins `u,x`, let `d` be their Hamming distance.  Then

```math
H_b(x)-H_b(u)\le2d(n-d)\le2nd.                                \tag{PR.15}
```

The field `nu` loses exactly `2nd` between `u` and `x`.  Hence `u` is a
maximizer in (PR.1) and

```math
\boxed{R_b(nu)=n^2+H_b(u).}                                   \tag{PR.16}
```

Substituting (PR.16) into the simultaneous evaluation gap proves (PR.4),
for example with `delta=3/16`.  If two landscapes shared one summary state,
one decoder value at query `nu` would have to approximate responses
separated by (PR.4), proving the information lower bound.

## 6. What this proves and what it does not

This theorem is a genuine all-order near-minimizer statement.  It needs no
conference, Walsh, spectral, or product structure, and it remains true for
`kappa` decreasing as fast as a sufficiently large constant times
`n^(-1/2)`.  It proves that **unrestricted amplitude-`n` future queries do
not become information-light merely because the cap is nearly minimal**.

It does not yet falsify the live balanced-interface route.  Equation
(PR.16) is full Boolean pinning, and an all-spins-free exact-sign compiler
which forces this field can carry a common `Theta(n^2)` energy baseline.
The theorem therefore does not show an `Omega(n)` response rate under the
more restrictive class of continuations whose entire parent cap remains
`O(n^(3/2))`.  Any use in the convergence program must either:

1. compile the same packing into such low-cap parents; or
2. prove that the intended composition theorem really must answer the
   amplitude-`n` pinned query family.

The correct frontier label is accordingly:

```text
PROVES AN ARROW for unrestricted contextual incompressibility at Level 5;
WEAKENS L_pack to the low-cap/balanced physical-context restriction;
does not prove incompressibility for the convergence-relevant interface.
```
