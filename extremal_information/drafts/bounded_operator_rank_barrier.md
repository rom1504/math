# Bounded-operator sign bridges have extensive spectral interface

Status: rigorous theorem draft.  This closes the rank-growth benchmark for
spectral truncation only; it is not an information lower bound against
nonlinear symmetry quotients.

## Theorem BR.1 (scale-visible rank barrier)

Let `R in {+-1}^{n by n}` and suppose

```math
||R||_(2 to 2)<=C sqrt(n),\qquad C>=1.                 \tag{BR.1}
```

For `0<=epsilon<1`, let

```math
r_epsilon(R)=#{j:sigma_j(R)>epsilon sqrt(n)}.          \tag{BR.2}
```

Then

```math
r_epsilon(R)>=
n{1-epsilon^2\over C^2-epsilon^2}.                    \tag{BR.3}
```

In particular, every fixed-error spectral approximation at the natural
`n^(3/2)` Boolean interaction scale retains `Omega_(C,epsilon)(n)` singular
features.  More generally, if `rank(R)<=r`, then

```math
||R||_(2 to 2)>= {n\over sqrt r}.                      \tag{BR.4}
```

Thus a dense sign bridge cannot simultaneously have sublinear algebraic
rank and bounded `sqrt(n)` operator scale.

### Proof

The sign constraint gives the exact Frobenius identity

```math
n^2=||R||_F^2=sum_(j=1)^n sigma_j(R)^2.                \tag{BR.5}
```

There are `r_epsilon` summands at most `C^2n` and the remaining summands are
at most `epsilon^2n`.  Therefore

```math
n^2<=r_epsilon C^2n+(n-r_epsilon)epsilon^2n,
```

which rearranges to (BR.3).  If only `r` singular values are nonzero, then
`n^2<=r||R||^2`, proving (BR.4). `square`

## Corollary BR.2 (the SVD bridge hierarchy has no intermediate regime)

In the balanced setting of Theorem 18.7, replacing `R` by its rank-`r`
truncation gives the certified worst-case Boolean response upper bound

```math
n sigma_(r+1)(R).                                      \tag{BR.6}
```

For bounded-operator sign bridges (BR.1), certifying error at most
`epsilon n^(3/2)` through this operator-norm/SVD bound requires retaining at
least the number of features in (BR.3).  Hence the regimes

```text
fixed rank -> slowly growing rank -> arbitrary bounded-op dense signs
```

do not interpolate through a subextensive SVD interface at fixed normalized
accuracy: the final regime jumps to linear visible rank.

This conclusion is deliberately scoped.  A common permutation orbit can
have full algebraic rank and still admit a polynomial response quotient by
rearrangement synchronization (Theorem 18.6).  The theorem therefore rules
out **rank truncation as the explanation**, not strict response compression
by other algebraic structure.

## Sharpness and falsifiers

1. A Hadamard matrix has `C=1`, all singular values `sqrt(n)`, and equality
   in (BR.3) for every `epsilon<1`.
2. The all-ones matrix has rank one and operator norm `n`, attaining (BR.4).
   It is nevertheless response-compressible by magnetization.  Thus replacing
   the conclusion by a general contextual-information lower bound is false.
3. For arbitrary real weighted bridges the Frobenius identity (BR.5) is
   absent.  The sign alphabet is essential.

## General lesson

At the original interaction scale, spectral tail control and response-state
compression are different resources.  Bounded operator norm prevents a few
singular directions from carrying the dense sign mass; only synchronization,
symmetry, or another nonlinear congruence could then compress a linear-rank
interface.
