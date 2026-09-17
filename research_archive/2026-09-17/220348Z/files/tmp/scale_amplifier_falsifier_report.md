# A public continuation cannot amplify the mesoscopic affine shell code

## Status

Task-local theorem/no-go report.  No canonical repository file is changed.
The result has two parts:

1. every quadratic continuation on `m` new spins is only `2m`-Lipschitz
   per changed old spin, even with arbitrary interaction among the new
   spins;
2. more decisively, the affine sparse-flip children used by AO.2 are already
   `o(n^(3/2))` apart in the Boolean quadratic norm when their source shell
   words are only `O(M_n)` apart.  Since a common continuation is
   one-Lipschitz in the child landscape, **no continuation whatsoever** can
   restore a leading gap for those children.  This conclusion does not use
   a bound on the number of new spins or on the continuation cap.

Thus the current energy-scale shell packing cannot be amplified by merely
replacing the rank-one AO shore while retaining the AO affine child
encoding.  An escape must change the state-dependent old block itself (or
obtain fixed-edge-scale shell directions), not just improve the public
future.

## 1. Notation

For a hollow real symmetric old block `D`, put

```math
 H_D(x)=\sum_{i<j}D_{ij}x_ix_j,
 \qquad
 \|D\|_{\mathrm B}=\max_{x\in\{\pm1\}^n}|H_D(x)|.
```

For an old child `A`, an old--new block
`B\in[-1,1]^{n\times m}`, and an arbitrary even new-spin landscape `K`,
write

```math
 \mathcal Q(A;B,K)
 =\max_{x,y}|H_A(x)+x^TBy+K(y)|.                 \tag{SA.1}
```

Every quadratic auxiliary block gives an even `K`.  Exact signs are a
special case of the displayed coefficient bound.

## 2. Two universal non-amplification inequalities

### Lemma SA.1 (public-continuation nonexpansiveness)

For every pair of old children `A,A'` and every common future `(B,K)`,

```math
 |\mathcal Q(A;B,K)-\mathcal Q(A';B,K)|
 \le \|A-A'\|_{\mathrm B}.                       \tag{SA.2}
```

The same inequality holds for a one-sided maximum or minimum.

#### Proof

The two full landscapes differ pointwise by `H_(A-A')(x)`.  Maximum,
minimum, and maximum absolute value are one-Lipschitz in uniform norm.
`square`

This permits arbitrary interaction and cancellation among all auxiliary
spins.  Their number and their internal cap do not enter.

### Lemma SA.2 (sublinear shores cannot select mesoscopically close cuts)

For an outer channel `sigma\in{\pm1}`, define the effective future

```math
 G_\sigma(x)=\max_y \sigma\{x^TBy+K(y)\}.          \tag{SA.3}
```

If `d_proj(x,x')=min(d_H(x,x'),n-d_H(x,x'))`, then

```math
 |G_\sigma(x)-G_\sigma(x')|
 \le 2m\,d_proj(x,x').                             \tag{SA.4}
```

For opposite outer channels one also has

```math
 |G_+(x)-G_-(x')|
 \le2m\,d_proj(x,x')+2\|K\|_\infty.               \tag{SA.5}
```

#### Proof

For the same channel, maxima of the two affine families differ by at most

```math
 \max_y |(x-x')^TBy|\le2m d_H(x,x').
```

Because `K` is even, `G_sigma` is even in `x`, which permits the projective
minimum.  After changing `y` to `-y`, the two functions at opposite outer
channels differ only by replacing `K` with `-K`; this costs at most
`2||K||_infinity`. `square`

Let `z=\sigma c(x)` and `z'=\sigma c(x')` have the same positive-energy
orientation, and put `h=d_H(z,z')`.  If `h<E/2`, then

```math
 h=d(n-d),\qquad d=d_proj(x,x')\le {2h\over n}.    \tag{SA.6}
```

Consequently

```math
 |G_\sigma(x)-G_\sigma(x')|\le {4mh\over n}.      \tag{SA.7}
```

Suppose both are in a positive shell of deficit at most `s`, so their old
oriented energies differ by at most `s`.  Their completed **anchor scores**
therefore obey

```math
 \left|\{\sigma H_A(x)+G_\sigma(x)\}
       -\{\sigma H_A(x')+G_\sigma(x')\}\right|
 \le s+{4mh\over n}.                              \tag{SA.8}
```

In the hard mesoscopic regime `s=o(n^(3/2))`, `h=O(M_n)=O(n^(3/2))`, and
`m=o(n)`, the right side is `o(n^(3/2))`.  A growing shell family contains a
growing same-orientation subfamily by pigeonhole.  Hence no sublinear shore,
regardless of its internal interaction or cap, can give leading relative
margin to the designated shell witnesses.  A leading optimized gap must
evacuate at least one designated witness and use genuinely different old
optimizers, or must already be present in state-dependent old-block
coefficients.

The width exponent is sharp for literal witness selection: (SA.8) requires

```math
 m\ge {n(Delta-s)\over4h}                          \tag{SA.9}
```

to create anchor-score margin `Delta`; when `h=Theta(n^(3/2))` and
`Delta=Theta(n^(3/2))`, this is `m=Omega(n)`.  The amplitude-`n` coordinate
pin in PR.1 uses exactly this linear-width scale.

## 3. The decisive collision with the AO affine child encoding

The preceding anchor result still allows a remote optimizer.  The next
result does not.

Let `a` be the base signing and let `z^u` be same-orientation augmented-cut
words.  Suppose exact-sign children `b^u` satisfy, uniformly over all
augmented cuts `z`,

```math
 |\langle b^u,z\rangle-
   \{(1-p)\langle a,z\rangle+p\langle z^u,z\rangle\}|
 \le\rho.                                          \tag{SA.10}
```

This is exactly the deterministic event (AO.6), independently of how it is
obtained.

### Theorem SA.3 (affine mesoscopic children are publicly invisible)

If `h_uv=d_H(z^u,z^v)<E/2`, then

```math
 \boxed{\|b^u-b^v\|_{\mathrm B}\le2p h_{uv}+2\rho.} \tag{SA.11}
```

Consequently, for **every** common exact-sign physical future, of arbitrary
order and with arbitrary auxiliary interaction,

```math
 |\mathcal Q(b^u;B,K)-\mathcal Q(b^v;B,K)|
 \le2p h_{uv}+2\rho.                               \tag{SA.12}
```

The bound holds separately for every member of an arbitrary common query
bank.

#### Proof

Subtract (SA.10) for `u` and `v`.  For every augmented cut `z`,

```math
 |\langle b^u-b^v,z\rangle|
 \le p|\langle z^u-z^v,z\rangle|+2\rho
 \le p\|z^u-z^v\|_1+2\rho
 =2ph_{uv}+2\rho.
```

Maximize over `z` to get (SA.11), then apply Lemma SA.1. `square`

For AO.2,

```math
 p={\alpha\over\sqrt n},
 \qquad
 \rho=O(\sqrt\alpha\,n^{5/4}+n).                  \tag{SA.13}
```

If the new exact-minimizer packing remains in the hard branch
`h_uv=O(M_n)=O(n^(3/2))`, then

```math
 2ph_{uv}+2\rho
 =O(\alpha n+\sqrt\alpha\,n^{5/4}+n)
 =o(n^(3/2)).                                      \tag{SA.14}
```

Thus changing the rank-one AO shore to a more sophisticated joint,
interacting, exact-sign continuation cannot possibly help while (AO.6) is
retained.  The child response metric has already collapsed before the
future is attached.

There is a useful architecture-level variant.  Suppose instead that
`b^u` uniformly approximates the affine interpolation
`(1-p)a+pz^u` to Boolean error `rho=o(n^(3/2))` and has
`Q(b^u)=O(n^(3/2))`.  Evaluating at `z^u` shows

```math
 pE\le O(n^(3/2))+\rho,
```

so `p=O(n^(-1/2))`.  Equation (SA.11) again gives `o(n^(3/2))` pair
distance whenever `h_uv=O(n^(3/2))`.  Therefore *every bounded-cap affine
interpolation encoder*, not only the independent Bernoulli implementation,
has the same mesoscopic ceiling.

## 4. Archive collision map

| Archived result | What it rules out | Relation to SA.1--SA.3 | Remaining escape |
|---|---|---|---|
| UP.1, universal pin barrier | one future robustly pins a prescribed old state against every complete-sign child; forces quadratic cap | SA.2 is narrower in child class but stronger in auxiliary interaction/cap for close-anchor margin: it needs only the cross-edge Lipschitz law | child-dependent optimizer switching on a restricted family |
| CD.3--CD.4, exact copy/character lock ceiling | fixed one-layer bridge preserving the full pair-character algebra loses `Omega(n^(3/2))` | compatible; SA.2 needs no character-preserving copy and applies to arbitrary interacting `K` | do not copy the old state; use an anti-pin/remote optimizer |
| SC.3--SC.5, disjoint-star and universal selector barriers | separately paid stars and uniform approximation of the all-positive cut shell | SA.3 is not a scalar-channel argument: arbitrary joint auxiliary cancellation is allowed and then removed by public nonexpansiveness | state-dependent old/cross coefficients outside the affine encoder |
| TC.1 / Axiom 61, public interaction cannot amplify absent state | public terms cancel when state enters only through bounded onsite blocks | SA.1 is the one-old-block specialization at the exact Boolean quadratic norm | broadcast the shell label into a leading old-block metric |
| BCX.2--BCX.3, bounded-cap anti-pin | a structured Hadamard switching family *can* be exposed with `sqrt(n)` new spins by optimizer switching | shows why no universal no-go is possible; BCX children are already `Theta(n^(3/2))` apart in the old Boolean norm and violate SA.3's mesoscopic affine premise | find analogous pre-existing old-block separation for near-minimizers |
| OV.1--OV.3, orientation visibility | opposite old signs can be exposed with an `n^(3/4)` high-internal-cap shore | irrelevant after taking a growing same-orientation shell subfamily; SA.2 then has no `Q(C)` term | orientation is not the hard-branch amplifier |
| CP.2 / projective response modulus | a particular spherical Gram response has a square-root modulus | SA.3 acts before relaxation and is exact for every common continuation | abandon the affine AO children or force fixed ambient separation |

Ordinary polarization and same-map Krivine are not used.  The estimate pays
no scalar channels separately; it allows the complete joint auxiliary
optimization and cancels it only because it is public to the two children.

## 5. Sharp remaining lemma

The hard branch cannot honestly be called a **future-amplification** problem
inside the current encoding.  A public future never increases old-child
uniform distance.  The first missing statement is instead a broadcast/state
construction:

> **Mesoscopic old-block broadcast (`L_broadcast`).**  From the positive
> `o(M_n)` shell of an exact minimizer, construct `K_n->infinity` exact-sign
> bounded-cap child blocks `A^1,...,A^K`, by a rule not using target-order
> optimization, such that
> `||A^i-A^j||_B >= c n^(3/2)` for all `i!=j`, while retaining a finite
> certificate that lets a common exact-sign query bank expose a constant
> fraction of this metric.

The norm-packing part is strictly weaker than parent maximization: it is a
same-support Boolean quadratic norm of pair differences, with no auxiliary
variables, no target-order optimum, and no requirement to compute `M_n`.
It can be weakened further to an explicitly certified version: for every
pair supply one displayed spin `x_(ij)` with
`|H_(A^i-A^j)(x_(ij))|>=c n^(3/2)`.  Verification is then a single quadratic
evaluation rather than a Boolean maximization.
It is also necessary for every public continuation by SA.1.  The additional
exposure certificate cannot be omitted—UP.1 shows that universal pinning is
too expensive, whereas BCX shows that a restricted anti-pin certificate can
work.

There are now only three honest exits from the hard branch:

1. prove a fixed-ambient shell direction (then AO.2 already works);
2. prove `L_broadcast`, necessarily abandoning the bounded-cap affine
   interpolation (AO.6);
3. allow child--query jointly owned coefficients, which is not a modular
   continuation and therefore does not count as a public compiler.

This is a scalable no-go for the live affine implementation, not a no-go for
all possible near-minimizer contextual packings.
