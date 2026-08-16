# Multichannel mixed holonomy has full macroscopic response rate

**Status.** The theorem below is proved analytically.  The finite identities
and a nontrivial `[16,5,8]` instance are checked by
[`verify_phase3_multichannel_holonomy_packing.py`](../experiments/verify_phase3_multichannel_holonomy_packing.py).

The one-channel packing in Theorem MC.6 of
[`phase3_mixed_circuit_hierarchy.md`](phase3_mixed_circuit_hierarchy.md)
shows that one new mixed circuit can carry linearly many response bits.  This
note answers the corresponding multichannel question.  If composition creates
`k` independent mixed circuits with `D`-bit holonomies, then, throughout the
regime `k <= D/32`, a constant fraction of all `Dk` raw holonomy bits remains
necessary even when:

* the only queries are **unlabelled kernel endpoints**;
* both input fragments are individually shear-trivial; and
* answers need only have a fixed macroscopic additive accuracy.

The proof is a packing of subspaces inside one asymptotically good binary
linear code.  It is not a dimension count on labelled channels: every pair is
separated by an actual shortest-path query `u in F_2^D`.

## 1. The two-fragment family

Let

```math
W=\mathbb F_2^D,\qquad Q=\mathbb F_2^k,
```

with coordinate bases `B={b_1,...,b_D}` of `W` and
`{q_1,...,q_k}` of `Q`.  For an ordered independent tuple

```math
V=(v_1,\ldots,v_k)\in W^k,
```

form two fragments

```math
P=B\cup\{(0,q_j):1\le j\le k\},\qquad
R_V=B\cup\{(v_j,q_j):1\le j\le k\}.          \tag{MP.1}
```

The quotient columns within either fragment are independent.  Therefore each
offset assignment is removable by a linear shear, and the two fragments have
the same separate kernel-fixing gauge classes for every `V`.  In their union,
however, the two lifts over each `q_j` form a mixed parallel-pair circuit with
holonomy `v_j`.

For `u in W`, quotient elimination gives the rooted kernel-endpoint profile

```math
F_V(u)
=\ell_{P\cup R_V}(u,0)
=\min_{J\subseteq[k]}
  \left(2|J|+\left|u+\sum_{j\in J}v_j\right|\right).          \tag{MP.2}
```

Let `C_V=span{v_1,...,v_k}` and write `d(u,C)` for Hamming distance
from `u` to a code `C`.  Since `V` is a basis, (MP.2) immediately implies

```math
d(u,C_V)\le F_V(u)\le d(u,C_V)+2k.            \tag{MP.3}
```

The upper bound is deliberately crude but uniform: represent a nearest
codeword using at most `k` basis vectors.

### Lemma MP.0 (response geometry is Grassmannian Hausdorff geometry)

For two `k`-dimensional subspaces `C,C' <= W`, equipped with arbitrary
ordered bases `V,V'`, let `d_H(C,C')` be their Hausdorff distance in the
ambient Hamming metric.  Then

```math
\left|\,\|F_V-F_{V'}\|_\infty-d_H(C,C')\,\right|\le2k.       \tag{MP.3a}
```

#### Proof

In any finite metric space, distance-to-set functions satisfy

```math
\|d(\mathord\cdot,C)-d(\mathord\cdot,C')\|_\infty=d_H(C,C').
```

The upper bound is the triangle inequality.  For the reverse bound, evaluate
at a point of either set attaining the directed Hausdorff distance.  By
(MP.3), write `F_V=d(.,C)+e` and `F_{V'}=d(.,C')+e'`, where both
`e,e'` take values in `[0,2k]`.  Therefore
`|(e-e')(u)|<=2k` pointwise, and the reverse triangle inequality for the
sup norm proves (MP.3a). `square`

Thus, when `k=o(D)`, the entire rooted response metric on this family is
asymptotically the Hamming--Hausdorff metric on the binary Grassmannian.  The
packing below is not an isolated witness construction; it is a packing of
that intrinsic response geometry.

## 2. A good ambient code

We record a self-contained weak Gilbert--Varshamov estimate with comfortable
constants.

### Lemma MP.1 (one good linear host)

For all sufficiently large `D`, there is a binary linear code
`C_0 <= F_2^D` such that

```math
\dim C_0=\lfloor D/4\rfloor,
\qquad d(C_0)>D/8.                             \tag{MP.4}
```

#### Proof

Let `r=floor(D/4)` and choose an `r`-dimensional subspace uniformly from the
Grassmannian.  A fixed nonzero vector belongs to it with probability

```math
\frac{2^r-1}{2^D-1}\le 2^{r-D+1}.             \tag{MP.5}
```

The number of nonzero vectors of weight at most `floor(D/8)` is at most

```math
\sum_{i\le D/8}{D\choose i}
\le 2^{H_2(1/8)D},                            \tag{MP.6}
```

where `H_2` is binary entropy.  Thus the expected number of such vectors in
the random subspace is at most

```math
2^{(H_2(1/8)-3/4)D+1}=o(1),                  \tag{MP.7}
```

because `H_2(1/8)=0.5435...<3/4`.  Some subspace contains none, which is
(MP.4). `square`

No effective construction of `C_0` is required for the information lower
bound.  Standard explicit asymptotically good codes could replace this
probabilistic host.

## 3. The full-rate response packing

### Theorem MP.2 (macroscopic `Dk` response packing)

For all sufficiently large `D` and every integer

```math
1\le k\le\lfloor D/32\rfloor,                 \tag{MP.8}
```

there is a family `mathcal V_{D,k}` of independent ordered `k`-tuples in
`F_2^D` such that

```math
|\mathcal V_{D,k}|\ge 2^{3Dk/16}              \tag{MP.9}
```

and, for distinct `V,V'` in the family,

```math
\boxed{\ \|F_V-F_{V'}\|_\infty>D/16.\ }       \tag{MP.10}
```

Every profile is realized by the union of the two individually
shear-trivial fragments in (MP.1).

#### Proof

Take `C_0` from Lemma MP.1 and put `r=dim C_0`.  For every `k`-dimensional
subspace `C <= C_0`, choose one ordered basis `V_C`.  The number of choices is
the Gaussian binomial coefficient, and

```math
{r\brack k}_2
=\prod_{i=0}^{k-1}\frac{2^r-2^i}{2^k-2^i}
\ge 2^{k(r-k)}.                               \tag{MP.11}
```

Indeed, every factor in the product is at least `2^{r-k}`.  When `D>=32`
and (MP.8) holds,

```math
r-k\ge D/4-1-D/32\ge 3D/16,                  \tag{MP.12}
```

which proves (MP.9).

Now take distinct subspaces `C,C' <= C_0` and choose
`c in C setminus C'`.  For every `c' in C'`, the vector `c+c'` is a
nonzero word of `C_0`; hence

```math
d(c,C')\ge d(C_0)>D/8.                        \tag{MP.13}
```

At the single unlabelled endpoint query `u=c`, (MP.2)--(MP.3) give

```math
F_{V_C}(c)\le2k,
\qquad
F_{V_{C'}}(c)\ge d(c,C')>D/8.                \tag{MP.14}
```

Since `2k<=D/16`, their difference is greater than `D/16`.  This proves
(MP.10).  Formula (MP.2) already established the two-fragment realization.
`square`

The same proof gives the flexible asymptotic statement: if a binary
`[D,r,d]` host has `r>=RD` and `d>=delta D`, then its `k`-subspaces give

```math
\log_2 |\mathcal V|\ge k(RD-k),
\qquad
\|F_V-F_{V'}\|_\infty\ge\delta D-2k.          \tag{MP.15}
```

Thus any asymptotically good code yields a constant-rate `Dk` packing for
all `k<=cD`, with `c` below both its rate and half its relative distance.

### Corollary MP.3 (deterministic and mutual-information charges)

Fix `epsilon<1/32`.  Any deterministic summary and decoder which uniformly
answer **every** kernel-endpoint query of (MP.2) to error at most `epsilon D`
must have, on this family,

```math
\log_2|\operatorname{range}S|\ge 3Dk/16       \tag{MP.16}
```

for all sufficiently large `D`.  Indeed, profiles assigned the same state
would be within `2epsilon D<D/16`, contrary to (MP.10).

More generally, let `V` be uniform on `mathcal V_{D,k}` and let a randomized
message `S` permit reconstruction of the entire profile to sup error
`epsilon D` with probability at least `1-eta`.  Nearest-profile decoding and
Fano's inequality give

```math
I(V;S)\ge(1-\eta)\log_2|\mathcal V_{D,k}|-H_2(\eta)
\ge (1-\eta)3Dk/16-H_2(\eta).                \tag{MP.17}
```

Here the decoder need not reveal labels for the `k` circuits.  The separated
objects are functions on the common query set `W`.

## 4. Interpretation and scope

The exact mixed-holonomy space in this example is `Hom(F_2^k,W)`, with
`Dk` binary degrees of freedom.  Theorem MP.2 proves that its macroscopic
response quotient still has `Theta(Dk)` bits in a nontrivial linear regime.
Consequently, the `D kappa` gluing count in MC.4 is not merely exact-state
overhead that disappears at coarse response resolution.  A constant fraction
can remain query-visible without labelled-channel access.

The mechanism is code separation.  The response profile is a bounded-cost
distance transform of the holonomy image code.  Distinct subspaces inside one
good host code contain a word macroscopically far from the other subspace, and
that word itself is the rooted witness.  This is more than ordinary convex
duality or a restatement of rate--distortion: the theorem identifies a
composition-created algebraic feature, constructs exponentially many
locally gauge-indistinguishable inputs, and exposes each difference through
the original shortest-path query family.

The theorem also has a direct systematic-code reading.  Appending two sets of
columns with the same independent quotient syndromes creates `k` parity
relations; the systematic syndromes of those relations form `C_V`.  Even
coarse prediction of all kernel-syndrome coset-leader weights requires a
constant fraction of the complete relation-syndrome matrix.

What the theorem does **not** say is that every high-nullity composition is
incompressible.  Structured holonomy images (low-rate families, synchronized
sections, or codes with special decoders) can have smaller response quotients.
It says that no general compression law may replace `Dk` by `o(Dk)` solely
because answers are approximate and channels are unlabelled.  Any such law
needs an additional structural hypothesis.
