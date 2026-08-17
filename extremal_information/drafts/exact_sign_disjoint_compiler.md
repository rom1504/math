# A compiler hierarchy and a character-preserving lock ceiling

**Status.**  Rigorous task-local draft.  This note freezes the weakest
compiler notion that transfers a positive response-information rate and then
closes one natural stronger architecture.  It complements
`algebraic_exact_sign_locking.md`: query-dependent coordinate pins solve the
weak metric problem, while the theorem below shows that a fixed one-layer
bridge cannot robustly lock a copy of the whole quadratic-character algebra.

## 1. The weakest rate-preserving notion

Let `H_z` be landscapes on `Omega_k`, and let

```math
q_t(z)=\max_{x\in\Omega_k}\{H_z(x)+K_t(x)\}          \tag{CD.1}
```

be a predeclared one-sided response language.  Its contextual metric is

```math
d_k(z,z')=\sup_t|q_t(z)-q_t(z')|.                   \tag{CD.2}
```

An **exact-sign modular compiler** at order `N` assigns a hollow complete
signing `P_(z,t)` to every child--query pair, together with a query-only
calibration `b_t`, subject to an edge-ownership rule: the sign of each parent
edge is a function of `z` alone, of `t` alone, or is public, but is never a
joint function of `(z,t)`.  Child and newly appended query vertices are
disjoint.  Its declared response is

```math
\widetilde q_t(z)=\max_s H_(P_(z,t))(s)-b_t.         \tag{CD.3}
```

It has distortion `eta` and gain `lambda>0` if

```math
|\widetilde q_t(z)-\lambda q_t(z)|\le\eta           \tag{CD.4}
```

uniformly in `(z,t)`.  A query-only calibration is legitimate here because
it cancels in every contextual distance.  It is not legitimate when the
quantity being minimized is the uncalibrated absolute cap of one parent.

### Lemma CD.1 (the minimal scale condition)

Suppose a subfamily `C_k` is pairwise separated by

```math
d_k(z,z')\ge\delta k^{3/2}.                          \tag{CD.5}
```

Then every compiler satisfying (CD.4) has compiled separation at least

```math
\lambda\delta k^{3/2}-2\eta.                        \tag{CD.6}
```

Consequently a gain bounded below, `eta=o(k^(3/2))`, and `N=O(k)` preserve
a positive response gap in units of `N^(3/2)` and preserve
`log_2|C_k|` information bits.  Conversely, if `N/k` tends to infinity while
the unscaled signal remains `O(k^(3/2))`, its normalized parent-scale gap
vanishes.

#### Proof

Choose a query witnessing (CD.5) up to an arbitrarily small error and apply
the triangle inequality twice in (CD.4).  This gives (CD.6).  The scale
claims are immediate. `square`

Thus a pointwise reconstruction of every overlay is stronger than necessary.
For response incompressibility it is enough to expose the metric by some
predeclared exact-sign queries.  The coordinate-pin compiler of EL.1 does
exactly this with `N=2k`.  By contrast, the sparse edge-variable identity

```math
\max_y\sum_{i<j}y_(ij)(x_i-T_(ij)x_j)
={k\choose2}-H_T(x)                                 \tag{CD.7}
```

uses `N=k+binom(k,2)` vertices, so its `k^(3/2)` signal has only
`N^(3/4)` size and fails the rate condition.

## 2. What it means to preserve the quadratic query algebra

A Boolean encoding `phi:{+-1}^k->{+-1}^k` is
**quadratic-character preserving** if, for every `a!=b`,

```math
phi_a(x)phi_b(x)=epsilon_(ab)x_(i_(ab))x_(j_(ab))   \tag{CD.8}
```

for a sign `epsilon_(ab)` and two distinct input coordinates.  This is the
minimal exact pullback property for all pair characters: it permits a common
state-dependent global gauge, which every quadratic query forgets.

### Lemma CD.2 (rigidity of exact quadratic-character pullback)

For `k>=5`, every quadratic-character-preserving encoding has the form

```math
phi(x)=g(x)DPx,                                     \tag{CD.9}
```

where `g:{+-1}^k->{+-1}` is arbitrary, `D` is a diagonal signing, and `P`
is a permutation matrix.

#### Proof

Ignore the signs in (CD.8) and denote the two-element support of the
character `phi_a phi_b` by `E_(ab)`.  Fix output coordinate zero and put
`F_a=E_(0a)`.  Character multiplication gives

```math
E_(ab)=F_a\mathbin\triangle F_b.                   \tag{CD.10}
```

The `F_a` are distinct two-subsets, and every two intersect in exactly one
point, because their symmetric difference is again a two-subset.  A family
of at least four pairwise-intersecting two-subsets is a star: after two sets
`{u,v}` and `{u,w}` are fixed, any set not containing `u` must be `{v,w}`,
and then no fourth distinct set exists.  Hence

```math
F_a={u,v_a}
```

with all `v_a` distinct.  There are `k-1` of them, so they exhaust the other
input coordinates.  The same cocycle identity for the signs gives
`epsilon_(ab)=s_as_b`.  Taking
`g(x)=phi_0(x)/(s_0x_u)` now yields (CD.9). `square`

The exceptional triangle at fewer than four nonbase coordinates explains
the harmless restriction `k>=5`.  Lemma CD.2 is a rigidity statement about
the full exact pair-character language, not about one structured subfamily
of quadratic queries.

## 3. A leading one-layer locking defect

Let `R in {+-1}^{k times k}` be a fixed complete bipartite bridge.  For an
old configuration `x`, unrestricted optimization of the new shore earns

```math
\max_y x^TRy=||R^Tx||_1.                            \tag{CD.11}
```

If the new shore is supposed to encode the same quadratic state through
`phi`, its global sign is irrelevant to every internal quadratic query.  The
best intended reward is therefore `|x^TRphi(x)|`.  Define the lock defect

```math
Delta_(R,phi)(x)=||R^Tx||_1-|x^TRphi(x)|\ge0.        \tag{CD.12}
```

### Theorem CD.3 (a fixed exact bridge cannot lock the full pair algebra)

For every `k>=5`, every complete sign bridge `R`, and every exact
quadratic-character-preserving encoding `phi`,

```math
\max_x Delta_(R,phi)(x)
\ge k\sqrt{k/2}-\sqrt3 k.                           \tag{CD.13}
```

In particular the worst locking defect is at least
`(1/sqrt2-o(1))k^(3/2)`.  It cannot be an
`o(k^(3/2))` stability error in a rate-preserving compiler.

#### Proof

By Lemma CD.2 and the absolute value in (CD.12), the gauge disappears and
the intended reward is `|x^TCx|` for another sign matrix
`C=RDP`.  For uniform Boolean `x`, every column of `R` has a Rademacher-sum
inner product, so the sharp `p=1` Khintchine lower bound gives

```math
E||R^Tx||_1
=k E|epsilon_1+...+epsilon_k|
\ge k\sqrt{k/2}.                                    \tag{CD.14}
```

On the other hand,

```math
x^TCx=tr C+\sum_{i<j}(C_(ij)+C_(ji))x_ix_j.
```

Orthogonality of cube characters and exact signs give

```math
E(x^TCx)^2
=(tr C)^2+\sum_{i<j}(C_(ij)+C_(ji))^2
\le k^2+4{k\choose2}<3k^2.                         \tag{CD.15}
```

Thus `E|x^TCx|<sqrt3 k`.  Average (CD.12) and use that its maximum is at
least its mean. `square`

This no-go allows the intended copy to choose its favorable antipode and so
is stronger than testing the literal duplicate `y=x`.  A scalar calibration
does not improve the selector gap (CD.12).  The theorem does **not** say
that every contextual response compiler must explicitly lock a copy; it
says that the standard proof plan “use one exact bridge as a dominant
equality penalty, then evaluate every quadratic future on the locked shore”
already loses the entire leading response scale.

### Theorem CD.4 (bounded coordinate replication does not help)

Let `m>=k`, let `pi:[m]->[k]` be onto with fibre sizes at most `L`, and let

```math
phi_a(x)=g(x)s_ax_(pi(a)).                           \tag{CD.16}
```

For every `R in {+-1}^{k times m}`,

```math
\max_x\{||R^Tx||_1-|x^TRphi(x)|\}
\ge m\sqrt{k/2}-\sqrt{m^2+2kLm}.                   \tag{CD.17}
```

Hence if `m=Theta(k)` and `L=O(1)`, the defect is still at least
`Omega(k^(3/2))`.

#### Proof

Again the gauge disappears.  Put

```math
C_(ij)=\sum_(a:pi(a)=j)R_(ia)s_a,
```

so the intended reward is `|x^TCx|`.  If `n_j=|pi^{-1}(j)|`, then

```math
||C||_F^2
\le k\sum_j n_j^2\le kLm,
\qquad |tr C|\le m.                                 \tag{CD.18}
```

The same Fourier calculation as (CD.15) yields

```math
E(x^TCx)^2
\le m^2+2||C||_F^2\le m^2+2kLm.                    \tag{CD.19}
```

Every one of the `m` bridge columns still contributes at least
`sqrt(k/2)` in expectation to the unrestricted roof.  Subtract the
Cauchy--Schwarz bound from that expectation. `square`

### Theorem CD.5 (universal robust one-layer pins are rank one)

Fix a target `u in {+-1}^k` and a complete sign bridge
`R in {+-1}^{k times k}`.  Put

```math
F_R(x)=\max_y x^TRy=||R^Tx||_1.                     \tag{CD.20}
```

Suppose `u` is a global maximizer of

```math
H_A(x)+F_R(x)                                       \tag{CD.21}
```

for **every** exact hollow sign child `A` satisfying
`||A||_(2->2)<=9sqrt(k)`.  Then

```math
R=u s^T                                             \tag{CD.22}
```

for a sign vector `s`.  Consequently

```math
||R||_(infinity->1)=k^2,                            \tag{CD.23}
```

and every complete quadratic parent containing this bridge has Boolean cap
at least `k^2`.

#### Proof

For every coordinate `i`, there is an exact hollow sign child `A_i` with

```math
u_i(A_iu)_i=-(k-1),
\qquad ||A_i||_(2->2)<=9sqrt(k).                    \tag{CD.24}
```

Prescribe the `i`th star by `A_(ij)=-u_i u_j`; this star matrix has norm
`sqrt(k-1)`.  A symmetric Rademacher completion of the remaining principal
block has norm at most `8sqrt(k)` for some realization, by the standard
`1/4`-net quadratic-form bound.  The triangle inequality proves (CD.24).

Let `u^(i)` flip coordinate `i`.  The child in (CD.24) gains exactly
`2(k-1)` under this flip.  Universal maximality of `u` therefore forces

```math
F_R(u)-F_R(u^(i))\ge2(k-1)                          \tag{CD.25}
```

for every `i`.  Write

```math
h_a=\sum_jR_(ja)u_j,
\qquad q_(ia)=R_(ia)u_i.
```

The contribution of column `a` to the left side of (CD.25) is

```math
|h_a|-|h_a-2q_(ia)|.                                \tag{CD.26}
```

It is `2` when `q_(ia)` has the sign of `h_a` and `|h_a|>=2`, is zero
when they agree and `|h_a|=1`, and is `-2` when they disagree or `h_a=0`.
Since the sum of `k` such terms is at least `2k-2`, no `-2` term is
possible, and each row has at most one zero term.  In particular every
`h_a` is nonzero and every `q_(ia)` has its sign, for every `(i,a)`.
Hence

```math
R_(ia)=u_i sign(h_a),
```

which is (CD.22) (and also forces `|h_a|=k`, so the apparent zero-term
exception cannot occur).  Evaluating the bridge at `(u,s)` proves (CD.23).

Finally, in an arbitrary quadratic parent, flipping the whole new shore
reverses the cross energy and preserves both internal energies.  One of the
two values has absolute magnitude at least the cross term.  Maximizing it
gives parent cap at least `||R||_(infinity->1)=k^2`. `square`

Thus the high `Theta(k^2)` calibration in the exact coordinate pin EL.1 is
not an accident within the bare one-layer architecture: universal robustness
against the natural spectrally flat child class forces precisely the
rank-one bridge and hence a quadratic cap.  CD.5 does not cover a correlated
auxiliary Hamiltonian; such an interaction can change the effective roof
away from (CD.20).

## 4. Consequences for the exact-sign closure target

There are now three sharply separated facts.

1. **Metric closure is possible.**  A query-dependent rank-one coordinate
   pin exposes the complete quadratic contextual metric with `2k` exact-sign
   vertices, albeit above a common `Theta(k^2)` calibration.
2. **Sparse pointwise closure is possible but loses rate.**  Edge selectors
   compile an oriented negative clone exactly with `Theta(k^2)` vertices.
3. **Universal one-layer equality locking fails at leading scale.**  A fixed
   complete bridge cannot both choose its natural optimizer and retain an
   exact copy of all quadratic pair characters with `o(k^(3/2))` error.
   If it must robustly pin one exposed state against every bounded-operator
   sign child, it is forced all the way to a rank-one bridge and quadratic
   parent cap.

Therefore the remaining pointwise/low-cap route must use at least one
ingredient absent from CD.3: a query-dependent pin, a nonlinear encoding
whose quadratic pullback is proved only for the narrower alternating-form
language, correlated interacting auxiliaries, or cancellation with the
child and query energies before the lock is separately optimized.  This is
exactly a joint same-switch obligation; replacing it by a stronger public
equality penalty is rigorously insufficient.
