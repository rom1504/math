# Marginal entropy cannot control the joint switching gain

Status: **proved scalable no-go for marginal entropy/rearrangement closures of
the exact switch-convolution certificate**.  Two pairs of strictly positive
Gibbs densities can have identical marginal value distributions (hence the
same Shannon, Rényi, Orlicz, and all scalar level-set data) while their best
translation gains differ by a linear amount.  The obstruction already occurs
on the rank-one switching group at every even group dimension.

This result does not obstruct a genuinely relational Fourier/subgroup
argument using the labelled alignment of the child and bridge kernels.  It
shows that such relational information is necessary: reverse
hypercontractivity, entropy transport, or rearrangement estimates fed only
the two marginal complexity profiles cannot yield the desired cross-order
cancellation.

## 1. Exact place in the cross-order certificate

Let `G` be the rank-one switching group and let `a,b:G->(0,infinity)` be
normalized densities,

```math
\mathbb E_Ga=\mathbb E_Gb=1.
```

In the exact bridge-switch identity from
`cross_order_exceptional_switch_convolution.md`,

```math
a={w_\epsilon\over\mathbb E w_\epsilon},
\qquad
b={k_B\over\mathbb E k_B},
```

and

```math
L_{\epsilon,g}(B)
=\log\mathbb Ew_\epsilon+\log\mathbb Ek_B
 +\log(a*b)(g),                                    \tag{1.1}
```

where convolution is normalized Haar convolution.  Consequently, for exact
own-scale child minimizers,

```math
\boxed{
E_{m,n}(\beta)
\le \log\mathbb Ew_\epsilon+\log\mathbb Ek_B
-P_m(\beta)-P_n(\beta)+\min_g\log(a*b)(g).}        \tag{1.2}
```

Thus a negative linear value of `min log(a*b)` is exactly the joint
cancellation which could remove a linear radial/bridge payment.  The theorem
below proves that no independently computed marginal entropy data determine
this term even to `o(|log |G||)` accuracy.

## 2. Complementary-subgroup construction

Write the switching group additively as

```math
G=\mathbb F_2^{2d}=H\oplus K,
\qquad \dim H=\dim K=d.                             \tag{2.1}
```

This is exactly `G_(m,n)` whenever `m+n-1=2d`; only the abstract group law is
used.  Let

```math
h_H=2^d\mathbf1_H,
\qquad h_K=2^d\mathbf1_K,                           \tag{2.2}
```

so both have Haar mean one.  For `0<eta<1`, form the strictly positive
two-level densities

```math
a_H=(1-\eta)h_H+\eta,
\qquad a_K=(1-\eta)h_K+\eta.                       \tag{2.3}
```

Each density has precisely the same value distribution:

```math
\eta+(1-\eta)2^d \quad\hbox{on a fraction }2^{-d},
\qquad
\eta \quad\hbox{on a fraction }1-2^{-d}.           \tag{2.4}
```

In particular `a_H` and `a_K` have identical:

- all positive and negative moments for which they are finite;
- Shannon and every Rényi divergence from Haar;
- decreasing rearrangements and all level-set cardinalities;
- unlabelled Fourier-coefficient multisets (an automorphism of `G` sends
  `H` to `K`).

The normalized subgroup densities satisfy the exact convolution identities

```math
h_H*h_H=h_H,
\qquad h_H*h_K=1.                                   \tag{2.5}
```

The second identity follows because every `g in G` has a unique
decomposition `g=h+k`.  Expanding (2.3) and using (2.5) gives

```math
\boxed{
a_H*a_K=1,}                                         \tag{2.6}
```

whereas

```math
\boxed{
a_H*a_H=(1-\eta)^2h_H+(2\eta-\eta^2).}             \tag{2.7}
```

Therefore

```math
\boxed{
\min_g\log(a_H*a_K)(g)=0,\qquad
\min_g\log(a_H*a_H)(g)=\log(2\eta-\eta^2).}      \tag{2.8}
```

Both inputs in the two pairs `(a_H,a_K)` and `(a_H,a_H)` have identical
marginal value data, but the joint translation gain is different.

Taking

```math
\eta_d=e^{-cd},\qquad c>0,                          \tag{2.9}
```

yields the scalable separation

```math
\boxed{
0-\log(2\eta_d-\eta_d^2)=cd-\log2+o(1).}           \tag{2.10}
```

Since `log|G|=2d log2`, this is a fixed positive fraction of the full group
dimension.

## 3. The reverse-KL/geometric-orbit version

Let `Pi_(a,b)` have density `a*b` relative to Haar.  The exact geometric
orbit identity uses

```math
D(U_G\Vert\Pi_{a,b})=-\mathbb E_G\log(a*b).         \tag{3.1}
```

For the complementary pair, (2.6) gives

```math
D(U_G\Vert\Pi_{a_H,a_K})=0.                         \tag{3.2}
```

For the aligned pair, (2.7) gives

```math
\begin{aligned}
D(U_G\Vert\Pi_{a_H,a_H})
=-&(1-2^{-d})\log(2\eta-\eta^2)\\
  &-2^{-d}\log\{(1-\eta)^22^d+2\eta-\eta^2\}.
                                                               \tag{3.3}
\end{aligned}
```

With (2.9),

```math
\boxed{
D(U_G\Vert\Pi_{a_H,a_H})=cd-\log2+o(1),}           \tag{3.4}
```

while every marginal divergence appearing in either pair is identical.
Thus the amount recovered by the geometric-orbit correction also cannot be
bounded, to sublinear accuracy, from any collection of separate scalar
entropies of the two kernels.

## 4. Formal method-class no-go

Call a proposed switching estimate **marginal-rearrangement based** if its
input is invariant under independent measure-preserving relabellings of
the values of `a` and `b`.  This class contains every estimate using only
their value distributions, all scalar `L^p`/Orlicz norms, all marginal
Rényi or Shannon entropies, effective support cardinalities, or any
combination of those data.

**Theorem 4.1.**  No marginal-rearrangement-based functional can estimate
either

```math
\min_g\log(a*b)(g)
\quad\hbox{or}\quad
D(U_G\Vert(a*b)U_G)                                 \tag{4.1}
```

with uniform additive error `o(log|G|)` over strictly positive densities on
finite rank-one switching groups.

*Proof.*  Feed the functional the two pairs in Section 2 with
`eta=e^{-cd}`.  Their marginal inputs agree exactly, whereas (2.10) and
(3.2)--(3.4) differ by `Theta(log|G|)`. `square`

Combining Theorem 4.1 with the direct arrow (1.2) gives the cross-order
conclusion.  A marginal entropy/rearrangement proof cannot certify the
linear negative correction needed to turn a linear bridge/radial bound
into `E_(m,n)=o(m+n)`: its uncertainty about that correction is itself
linear.  The missing datum is not another scalar entropy.  It is the
**labelled relative geometry** of the two densities (equivalently here,
the position of their Fourier/subgroup supports).

## 5. Scope

The theorem is deliberately a method-class obstruction, not a claim that
the two-level densities (2.3) are produced by actual pressure-minimizing
children and sign bridges.  A successful proof may exploit the special
quadratic origin of `w_epsilon` and `k_B`, but then it must control their
relative labelled geometry.  The following routes remain outside the
no-go:

1. a Fourier or representation-theoretic theorem tying the labelled
   spectra of the actual child and bridge kernels;
2. a joint transport/rearrangement inequality using their shared Boolean
   coordinates, rather than independent marginal summaries;
3. a law-specific signed interpolation controlling the actual Gibbs path.

What is ruled out is the broad proposed escape in which one replaces the
separate Hölder payments by separate entropy, effective-support, or
reverse-hypercontractive summaries and expects convolution alone to supply
a universal cancellation.  Even complete knowledge of both marginal
value distributions leaves a leading-order ambiguity.
