# A universal affine near-top algebra from low local fields

Date: 2026-08-17.

Status: proof draft.  This gives a genuine growing-interface response state
for every bounded-cap signing.  It is a one-step contextual theorem, not a
cross-order recurrence: its most unbalanced physical endpoints can add a
leading `n^(3/2)` amount.

## 1. The near-top affine coset

Let `A` be a hollow symmetric signing of order `n`.  Choose an absolute
ground state `x` and an independent sign `rho in {+-1}` so that

```math
\rho H_A(x)=Q(A)=Q>0.
```

Replace `A` by `rho A` for the proof and rename the oriented matrix `A`.
This global negation preserves the absolute cap and the absolute trust
response below.  (Replacing `x` by `-x` would not orient a quadratic form.)

Put

```math
\ell_i=x_i(Ax)_i.
```

One-spin optimality and the quadratic identity give

```math
0\le\ell_i\le Q,
\qquad \sum_i\ell_i=2Q.                                      \tag{LA.1}
```

For `S subset [n]`, write `x^S` for the spin obtained by flipping exactly
the coordinates in `S`.

### Theorem LA.1 (low-field affine shell algebra)

For every integer `1<=k<=n-1`, there is a `k`-element vertex set `I` such
that the projective affine coset

```math
\mathcal C_I=\{x^S:S\subseteq I\}                             \tag{LA.2}
```

has `2^k` members and every member obeys

```math
\boxed{
Q-|H_A(x^S)|
\le\Delta_k:={4kQ\over n}+2k(k-1).}                          \tag{LA.3}
```

Moreover `mathcal C_I` is closed under every odd coordinatewise product:
the product of `2r+1` members is another member.

#### Proof

Take `I` to be the indices of the `k` smallest local fields.  Then

```math
\sum_(i\in I)\ell_i\le {2kQ\over n}.                         \tag{LA.4}
```

If `s_ij=a_ijx_ix_j`, direct expansion gives

```math
H_A(x^S)
=Q-2\sum_(i\in S)\ell_i
  +4\sum_(\{i,j\}\subseteq S)s_ij.                           \tag{LA.5}
```

Since `Q-|t|<=Q-t` for every `t<=Q`, equations (LA.4)--(LA.5) imply

```math
Q-|H_A(x^S)|
\le2\sum_(i\in S)\ell_i+4\binom {|S|}2
\le {4kQ\over n}+2k(k-1).
```

For `k<=n-1`, distinct subsets of `I` remain distinct modulo global sign.
Finally, the product of an odd number of vectors `x^(S_j)` is
`x^(S_1 triangle ... triangle S_(2r+1))`, proving closure. `square`

For every bounded-cap sequence `Q(A_n)=O(n^(3/2))`, taking
`k=o(n^(3/4))` makes `Delta_k=o(n^(3/2))`.  In particular
`k=Theta(sqrt n)` gives an affine algebra with `2^(Theta(sqrt n))` near-top
spins and the much sharper absolute error `Delta_k=O(n)`.

## 2. A majority port frame whose whole orbit stays in the coset

Assume `k` is even, so `p=k+1` is odd.  Enumerate
`I={i_1,...,i_k}` and define the Boolean port frame

```math
W=(w_0,w_1,...,w_k),
\qquad w_0=x,
\qquad w_j=x^{\{i_j\}}.                                      \tag{LA.6}
```

Let `tau` be odd majority on `p` inputs.  For an endpoint
`epsilon in {+-1}^p`, put

```math
x_epsilon=\operatorname {sgn}(W\epsilon)
           =\tau(W\mathbin\odot\epsilon)                     \tag{LA.7}
```

coordinatewise.  No coordinate ties occur because `p` is odd.

### Lemma LA.2 (the nonlinear selector orbit is affine)

Every endpoint selector `x_epsilon` belongs projectively to `mathcal C_I`.
Every odd Fourier product of the ports in (LA.6) also belongs to
`mathcal C_I`.  Consequently all those selectors and products have
absolute energy at least `Q-Delta_k`.

#### Proof

Let

```math
t=\epsilon_0+\sum_(j=1)^k\epsilon_j.
```

Outside `I`, the row sum is `x_i t`.  At vertex `i_j`, it is
`x_(i_j)(t-2epsilon_j)`.  Relative to the common global sign `sgn(t)`, the
selector therefore differs from `x` only on a subset of `I`, proving the
first assertion.  (If one of the latter sums changes sign, record that
coordinate in the subset.)

For an odd set of port columns, the product contains an odd number of base
factors `x`; its remaining single-coordinate flips XOR to a subset of `I`.
This is exactly a member of (LA.2). `square`

## 3. Cap-relative response compression

For an integer field `g`, define the absolute trust response

```math
\mathcal B_A(g)
=\max_(y\in\{+-1\}^n,\ \sigma\in\{+-1\})
 \{\sigma H_A(y)+g\mathbin\cdot y\}.                          \tag{LA.8}
```

### Theorem LA.3 (uniform response formula on the affine interface)

For every integer `m>=0` and every endpoint `epsilon`,

```math
\boxed{
0\le
Q+m\|W\epsilon\|_1-\mathcal B_A(mW\epsilon)
\le\Delta_k.}                                               \tag{LA.9}
```

After switching by `diag(x)`, the projective row histogram of `W` has only
`k+1` types: one all-positive type of multiplicity `n-k`, and `k` types with
one negative coordinate, each of multiplicity one.  Thus the scalar
`||W epsilon||_1` and the complete uniform error certificate (LA.9) have an
`O(k log n)`-bit presentation rather than an `n by (k+1)` table.

#### Proof

The cap and Holder give the upper bound

```math
\mathcal B_A(mW\epsilon)\le Q+m\|W\epsilon\|_1.
```

The Boolean spin `x_epsilon=sgn(W epsilon)` realizes the entire field norm.
Lemma LA.2 gives `|H_A(x_epsilon)|>=Q-Delta_k`; choose `sigma` to orient its
quadratic energy positively.  This proves the lower bound.  Switching sends
`x` to the all-positive vector, after which (LA.6) has exactly the displayed
row types. `square`

With `k=Theta(sqrt n)` and `Q=O(n^(3/2))`, Theorem LA.3 has

```math
\text{port arity }p=Theta(\sqrt n),
\quad
\text{declared-frame state bits }O(\sqrt n\log n),
\quad
\text{uniform response error }O(n)=o(n^{3/2}).                \tag{LA.10}
```

This is a strict sub-landscape collective state for an already declared
projective frame at the target interface width and does not use conference,
Walsh, tensor, or near-minimizer structure.  Reconstructing the labelled
frame independently of `A` additionally requires the projective ground-state
gauge (up to `n-1` bits).  The theorem survives the sparse-flip coherence
obstruction because it certifies the nonlinear selectors themselves, not
only their Fourier-product marginals.

## 4. Quantitative shell packing

### Corollary LA.4 (a universal target-distance affine shell packing)

Suppose `Q(A)<=C n^(3/2)`.  Fix `c>0`, and let `k` be the largest even
integer at most `c sqrt n` (and at most `n-1`).  The affine shell contains
`exp(Omega_c(sqrt n))` spins with one common energy orientation, all at
absolute deficit `O_(C,c)(n)`, whose corresponding oriented signed cut
words have pairwise edge-Hamming distance `Theta_c(n^(3/2))`.

Indeed, color every one of the `2^k` affine masks by a sign orienting its
energy.  One color class has at least `2^(k-1)` members.  Greedily delete
mask-Hamming balls of radius `floor(k/4)-1`; the entropy bound

```math
\sum_(j<k/4)\binom kj\le2^{h_2(1/4)k}
```

leaves `exp(Omega(k))` masks at mutual distance `d>=floor(k/4)`.  Two such
spins differ projectively on `d<=k` vertices, so their same-orientation cut
words differ on exactly

```math
d(n-d)=Theta(n^(3/2))
```

edges.  The common deficit is at most
`(4cC+2c^2+o(1))n` by LA.1.

This is target-scale **word distance**, not a contextual response packing:
the entire family still has one `k`-generator affine presentation.  It
therefore supplies structured shell geometry, not a lower bound on the bits
needed by a composable state.

## 5. Scope audit

The theorem is a real response-compression result, but it does not yet move
the cross-order convergence arrow.

1. The interface is **designed from one ground state and its low-field
   vertices**.  It does not approximate an arbitrary dense bridge.
2. The full endpoint language contains
   `||W epsilon||_1=Theta(np)`.  At `p=Theta(sqrt n)` this is a fixed leading
   `Theta(n^(3/2))` contribution.  There is also a selected `t=3` endpoint
   of norm exactly `3n-4` whose selector is `x`; a microcanonical compiler
   can expose that one anchor.  It collapses competing endpoints to one
   scalar segment and retains the ordinary `Theta(n sqrt p)` random-bridge
   error, rather than physicalizing the whole affine response language.
3. Recomputing a new low-field coset after composition gives no congruence
   between successive interfaces.  The response error is depth-independent
   for one declared interface, not automatically reusable.
4. The switching gauge can be applied to the child signing at no cap cost,
   but a labelled external application would otherwise need to know the
   ground-state switch.  The compact state is therefore constructive/gauged,
   not a gauge-free invariant of `A`.

The selected subleading endpoint is therefore understood.  The next exact
question is whether its `Theta(n sqrt p)` balanced residual can cancel
jointly against the child quadratic energy, or whether a globally correlated
cross-level construction can avoid one-step error accumulation.  Without
such a mechanism, LA.1--LA.4 are a Level-5 necessary/response theorem and a
new benchmark, not the missing recurrence.
