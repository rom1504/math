# A public continuation cannot amplify an already-collapsed affine shell code

**Status.** Rigorous task-local theorem, independently audited.  This is a
no-go for the affine child encoding used in AO.2, not for every possible
near-minimizer response encoding.

The result separates two issues which can otherwise be conflated.

1. A common continuation is nonexpansive in the uniform Boolean Hamiltonian
   distance between its children, even when it contains arbitrarily many
   auxiliary spins and arbitrary joint interaction among them.
2. If two AO affine children come from same-orientation shell words at
   `o(n^2)` Hamming distance (in particular at `O(M_n)` distance), then their
   Boolean Hamiltonian distance is already `o(n^(3/2))`.  No common future can
   turn that pair into a leading response gap.

The present fractional-reservoir theorem proves a *lower* separation
`(1/4-o(1))M_n`; it does not prove that the resulting family has
`O(M_n)` diameter.  Thus the theorem below kills amplification of the
genuinely mesoscopic branch while leaving open the possibility that those
words contain a fixed-ambient-scale subpacking.

## 1. Conventions

For a hollow real symmetric coefficient block `D` on `n` old spins, write

```math
 H_D(x)=\sum_{i<j}D_{ij}x_ix_j,
 \qquad
 \|D\|_{\rm B}=\max_{x\in\{\pm1\}^n}|H_D(x)|.
                                                        \tag{SA.1}
```

The coefficient vector of an augmented cut is

```math
 z=\sigma c(x),\qquad c(x)_{ij}=x_ix_j,
 \qquad \sigma\in\{\pm1\}.
```

Consequently

```math
 \|D\|_{\rm B}=\max_z\langle D,z\rangle.             \tag{SA.2}
```

For a common continuation landscape `J(x,y)`, with an arbitrary auxiliary
state space for `y`, define

```math
 \mathcal Q(D;J)=\max_{x,y}|H_D(x)+J(x,y)|.           \tag{SA.3}
```

An exact-sign quadratic future on `m` new spins has

```math
 J(x,y)=x^TBy+K(y),
 \qquad B\in\{\pm1\}^{n\times m},
 \qquad K=H_C,                                      \tag{SA.4}
```

where `C` is hollow.  The formulation (SA.3) is deliberately more general.

## 2. Universal public-continuation nonexpansiveness

### Lemma SA.1

For any two old children `D,D'` and every *common* continuation `J`,

```math
 \boxed{
 |\mathcal Q(D;J)-\mathcal Q(D';J)|
 \le \|D-D'\|_{\rm B}.}                            \tag{SA.5}
```

The same estimate holds with the outer maximum absolute value replaced by
a one-sided maximum or minimum.

#### Proof

The two full landscapes differ pointwise by `H_(D-D')(x)`, whose uniform
absolute value is `\|D-D'\|_B`.  Maximum, minimum, and maximum absolute
value are one-Lipschitz in uniform norm. `square`

Thus arbitrary optimizer switching and arbitrary joint cancellation among
the auxiliary variables are already included.  The optimizing `(x,y)` may
be completely different for `D` and `D'`; no common optimizer is assumed.

For a bank of common queries `J_q`, (SA.5) holds separately for every `q`:

```math
 \sup_q|\mathcal Q(D;J_q)-\mathcal Q(D';J_q)|
 \le\|D-D'\|_{\rm B}.                              \tag{SA.6}
```

## 3. A complementary literal-shore bound

The nonexpansiveness theorem is decisive once the old-child metric is
known.  The following bound is useful earlier, when one asks whether a
shore can give two designated old witnesses a large relative score.

Let `B\in[-1,1]^{n\times m}` and let `K` be any even new-spin landscape.
For the outer absolute-value channel `\sigma\in\{\pm1\}`, put

```math
 G_\sigma(x)=\max_y\sigma\{x^TBy+K(y)\}.            \tag{SA.7}
```

### Lemma SA.2

For all old spins `x,x'`,

```math
 |G_\sigma(x)-G_\sigma(x')|
 \le2m\,d_{\rm proj}(x,x'),                        \tag{SA.8}
```

where

```math
 d_{\rm proj}(x,x')=\min\{d_H(x,x'),n-d_H(x,x')\}.
```

For opposite outer channels,

```math
 |G_+(x)-G_-(x')|
 \le2m\,d_{\rm proj}(x,x')+2\|K\|_\infty.          \tag{SA.9}
```

#### Proof

For a fixed channel, comparison of the two maxima costs at most

```math
 \max_y|(x-x')^TBy|\le2m d_H(x,x').
```

Evenness of `K` implies `G_\sigma(-x)=G_\sigma(x)`, so one may choose the
better representative of the projective spin.  Also

```math
 G_-(x')=\max_y\{x'^TBy-K(y)\}
```

after replacing `y` by `-y`.  Comparing this with `G_+(x)` gives (SA.9).
`square`

Suppose `z=\sigma c(x)` and `z'=\sigma c(x')` have the same augmented
orientation.  If `d=d_proj(x,x')` and `h=d_H(z,z')`, then exactly

```math
 h=d(n-d),\qquad d\le{2h\over n}.                  \tag{SA.10}
```

Hence

```math
 |G_\sigma(x)-G_\sigma(x')|\le{4mh\over n}.        \tag{SA.11}
```

If both augmented words lie in the positive deficit-`s` shell of a signing
`a`, their old oriented energies differ by at most `s`, and their completed
designated-witness scores satisfy

```math
 \left|
 \{\sigma H_a(x)+G_\sigma(x)\}
 -\{\sigma H_a(x')+G_\sigma(x')\}
 \right|
 \le s+{4mh\over n}.                               \tag{SA.12}
```

Thus, when `s=o(n^(3/2))`, `h=O(n^(3/2))`, and `m=o(n)`, a sublinear shore
cannot give these *designated* witnesses a leading relative score.  A
leading optimized gap could still use remote old optimizers; Lemma SA.1 is
what closes that loophole once the entire old-child distance is small.

For a desired designated-witness margin `Delta`, (SA.12) requires

```math
 m\ge {n(\Delta-s)\over4h}.                        \tag{SA.13}
```

At `h,Delta=Theta(n^(3/2))`, literal shore selection therefore needs
`m=Omega(n)`.  The evenness hypothesis is essential only for the projective
reduction in this section; Lemma SA.1 needs no such hypothesis.

## 4. Collision of the AO affine children

Let `a` be a base signing and let `z^u` be positive augmented-cut words.
Suppose exact-sign children `b^u` obey, uniformly over every augmented cut
`z`,

```math
 \left|\langle b^u,z\rangle-
 \{(1-p)\langle a,z\rangle+p\langle z^u,z\rangle\}
 \right|\le\rho.                                   \tag{SA.14}
```

This is precisely the deterministic event (AO.6); its probabilistic origin
is irrelevant below.

### Theorem SA.3 (affine mesoscopic children are publicly invisible)

For every pair `u,v`, with `h_uv=d_H(z^u,z^v)`,

```math
 \boxed{
 \|b^u-b^v\|_{\rm B}\le2p h_{uv}+2\rho.}           \tag{SA.15}
```

Consequently every common continuation, of arbitrary order and arbitrary
internal interaction, obeys

```math
 |\mathcal Q(b^u;J)-\mathcal Q(b^v;J)|
 \le2p h_{uv}+2\rho.                               \tag{SA.16}
```

The same bound holds coordinatewise for every member of an arbitrary common
query bank.

#### Proof

Subtract (SA.14) for `u` and `v`.  For every augmented cut `z`,

```math
 |\langle b^u-b^v,z\rangle|
 \le p|\langle z^u-z^v,z\rangle|+2\rho
 \le p\|z^u-z^v\|_1+2\rho
 =2ph_{uv}+2\rho.
```

Maximizing over augmented cuts and using (SA.2) proves (SA.15).  Lemma SA.1
then proves (SA.16). `square`

For AO.2,

```math
 p={\alpha\over\sqrt n},
 \qquad
 \rho=O(\sqrt\alpha\,n^{5/4}+n).                  \tag{SA.17}
```

For fixed `alpha`, if `h_uv=o(n^2)`, then

```math
 2ph_{uv}+2\rho=o(n^{3/2}).                        \tag{SA.18}
```

In particular, if `h_uv=O(M_n)=O(n^(3/2))`, the first term is only
`O(alpha n)` and the total is `O(n^(5/4))` for fixed `alpha`.

The orientation restriction matters when one starts only with projective
distance.  For same-orientation words `z=\sigma c(x)` and
`z'=\sigma c(x')`,

```math
 d_H(z,z')\le d_{\rm P}(z,z')+\lfloor n/2\rfloor. \tag{SA.19}
```

Indeed `d_H(z,z')=d(n-d)\le\lfloor n^2/4\rfloor`, so its excess above its
edge complement is at most `floor(n/2)`.  Therefore a same-orientation
family of projective diameter `o(n^2)` also has actual diameter `o(n^2)`.
Every growing augmented-cut family has a growing same-orientation subfamily
by the two-colour pigeonhole principle.

### Corollary SA.4 (the exact scope for the fractional-reservoir packing)

Let `U_n` be any same-orientation subfamily of the FR.5 shell and construct
its AO children with fixed `alpha`.  If

```math
 \max_{u,v\in U_n}d_{\rm P}(z^u,z^v)=o(n^2),        \tag{SA.20}
```

then its response diameter under every bank of common public futures is
`o(n^(3/2))`.

FR.5 itself proves only

```math
 \min_{u\ne v}d_{\rm P}(z^u,z^v)
 \ge(1/4-o(1))M_n.                                 \tag{SA.21}
```

It does not prove (SA.20), nor its negation.  If a growing fixed-ambient
subpacking is present, the audited AO.2 compiler already applies to that
subpacking.  SA.4 instead rules out the proposal to retain (AO.6) and ask a
more elaborate public future to amplify pairs which remain only
mesoscopically separated.

## 5. A deterministic bounded-cap affine ceiling

The Bernoulli implementation is not essential.  Suppose, for one common
`p\in[0,1]`, that exact-sign children satisfy

```math
 \left\|b^u-\{(1-p)a+pz^u\}\right\|_{\rm B}\le\rho,
 \qquad Q(b^u)=O(n^{3/2}),                         \tag{SA.22}
```

where each `z^u` is positive for `a`, and `rho=o(n^(3/2))`.  Evaluation at
`z^u` gives

```math
 pE\le Q(b^u)+\rho=O(n^{3/2}),
 \qquad p=O(n^{-1/2}).                             \tag{SA.23}
```

The same subtraction as in SA.3 gives

```math
 \|b^u-b^v\|_{\rm B}\le2p h_{uv}+2\rho.           \tag{SA.24}
```

Thus every bounded-cap affine interpolation encoder with a common mixing
coefficient has `o(n^(3/2))` pair distance whenever `h_uv=o(n^2)`.  This is
a deterministic architecture-level ceiling, not a concentration artefact.

## 6. What is ruled out, and what is not

### Ruled out

* Replacing only the rank-one AO shore by a more complicated *common*
  quadratic continuation while retaining (AO.6).
* Using more auxiliary spins, an interacting auxiliary block, an
  unbounded-cap auxiliary block, or adaptive optimizer switching to amplify
  an AO pair whose source Hamming distance is `o(n^2)`.
* Any deterministic bounded-cap common-`p` affine interpolation encoder of
  mesoscopically close same-orientation shell words, followed by a common
  public continuation.

The argument does not use scalar-channel decomposition, ordinary
polarization, same-map Krivine rounding, or a separately paid left/right
bound.  It optimizes the entire common continuation jointly and then uses
uniform nonexpansiveness.

### Not ruled out

* A fixed-ambient-scale projective subpacking already present in an exact
  minimizer shell.
* A new non-affine child encoding whose child landscapes are themselves
  `Omega(n^(3/2))` apart in uniform Boolean norm.
* Broadcasting the shell label into any state-dependent child-owned block
  (old--old, owned auxiliary, or owned interface coefficients) before the
  common future is attached.
* A child--query jointly owned coefficient block.  This is not a modular
  common continuation and is outside contextual nonexpansiveness.
* Orientation-sensitive encodings using opposite augmented orientations;
  internal-future cap then matters as in OV.1--OV.3.
* A construction whose optimized response gap is already present before the
  public future, as in the BCX switching family.

The necessary metric precursor for any modular public compiler is therefore
best stated for the *entire child landscape*, not only its old block:

> **Mesoscopic state broadcast (`L_broadcast`).**  From the positive
> `o(M_n)` shell of an exact minimizer, construct `K_n\to\infty` exact-sign
> bounded-cap child landscapes `D^1,\ldots,D^{K_n}` on a common domain, by a
> rule not invoking target-order optimization, such that
> `\|D^i-D^j\|_\infty\ge c n^(3/2)` for every pair, together with a finite
> common-query certificate exposing a fixed fraction of that distance.

Here `\|D^i-D^j\|_\infty` is uniform Hamiltonian distance; when the state is
stored only in a hollow old block it is exactly `\|A^i-A^j\|_B`.  A still
more checkable version supplies, for every pair, one explicit configuration
at which the difference has magnitude at least `c n^(3/2)`.  This precursor
is necessary by SA.1 and is strictly weaker than maximizing the target
parent, but the exposure certificate remains essential: a large child
metric alone need not furnish a low-cap public query that realizes it.

The clean trichotomy is therefore:

1. find a growing fixed-ambient shell subpacking, to which AO.2 applies;
2. abandon the collapsed affine child state and prove a non-affine
   `L_broadcast`-type theorem;
3. abandon modular public continuations by allowing child--query jointly
   owned coefficients.

This is a scalable no-go for the live affine implementation, not a no-go for
all possible near-minimizer contextual packings.
