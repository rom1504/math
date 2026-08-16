# Mesoscopic induced sampling universalizes to iid signs

Date: 2026-08-16.

Status: **verified**, independently audited.  This note disproves the uniform
mesoscopic induced-submatrix recovery lemma proposed in
`ar_action_independent_proposals.md`.

## 1. Fixed samples of every operator-bounded signing are asymptotically iid

Let `A_N` be any deterministic symmetric hollow signing satisfying

```math
\|A_N\|_{op}\le C\sqrt N.                                \tag{MS.1}
```

Represent it by the signed step graphon `W_N` on `[0,1]^2`.  For measurable
sets `S,T`, bilinearity reduces the blockwise cut supremum to vertex subsets,
and Cauchy--Schwarz gives

```math
\left|\frac1{N^2}\mathbf1_S^{\mathsf T}A_N\mathbf1_T\right|
\le\frac{\|A_N\|_{op}\sqrt{|S||T|}}{N^2}
\le\frac C{\sqrt N}.                                    \tag{MS.2}
```

Thus `||W_N||_square->0`.  The bounded signed-graphon counting lemma implies,
for every fixed nonempty simple graph `F`,

```math
t(F,W_N)\longrightarrow0.                                \tag{MS.3}
```

The same is true for injective densities, since collisions have probability
`O_F(1/N)`.

For the cut/counting equivalence in this bounded dense setting, see
Borgs--Chayes--Lovasz--Sos--Vesztergombi, *Convergent sequences of dense
graphs I* ([arXiv:math/0702004](https://arxiv.org/abs/math/0702004)); the
signed version used here follows by the same telescoping argument, or by
writing `A=2G-(J-I)`.

Fix `m` and sample distinct uniform vertices `I_1,...,I_m`.  For any prescribed
edge-sign pattern `sigma in {+1,-1}^(E(K_m))`, expand

```math
\begin{aligned}
&\mathbb P\{a_{I_iI_j}=\sigma_{ij}\text{ for every }i<j\}\\
&\quad=2^{-\binom m2}
\sum_{F\subseteq E(K_m)}
 \left(\prod_{e\in F}\sigma_e\right)
 \mathbb E_{\rm inj}\prod_{e\in F}a_e.                  \tag{MS.4}
\end{aligned}
```

The empty term equals one and every other term tends to zero by (MS.3).
There are finitely many patterns at fixed `m`, so

> **Fixed-sample universality.**  For every fixed `m`, the induced signing
> `A_N[{I_1,...,I_m}]` converges in total variation, as `N->infinity`, to a
> symmetric hollow order-`m` matrix `G_m` with independent Rademacher upper
> edges.                                                         \(\tag{MS.5}\)

No action convergence of the parent sequence is used; the operator bound
alone forces (MS.5).

## 2. The iid cap is above the extremal upper scale

For `G_m`, expose vertices in order and choose spins greedily:

```math
x_1=1,
\qquad
x_j=\operatorname{sign}\left(\sum_{i<j}(G_m)_{ij}x_i\right).
```

The resulting energy is a sum of independent absolute simple random walks.
The standard expectation asymptotic and concentration give

```math
\frac{Q(G_m)}{m^{3/2}}
\ge\frac23\sqrt{\frac2\pi}-o_{\mathbb P}(1)
=0.5319230405\ldots-o_{\mathbb P}(1).                   \tag{MS.6}
```

This is strictly larger than the rigorous all-order upper scale `1/2`.
Equivalently, with probability tending to one,

```math
\Phi(T_{G_m})
\ge\frac43\sqrt{\frac2\pi}-o(1)
=1.0638460811\ldots-o(1).                                \tag{MS.7}
```

The proof of (MS.6), including independence of the sequential walk fields,
is also recorded in `exchangeable_recovery_obstruction.md`.

## 3. Failure of mesoscopic action recovery

Let a purified near-minimizing sequence satisfy (MS.1) and converge along a
subsequence to an action object `T`.  Since the all-order construction gives
`limsup M_n/n^(3/2)<=1/2`, its tolerance can be chosen so that

```math
\Phi(T)\le1+\eta                                         \tag{MS.8}
```

for arbitrarily small fixed `eta>0`.

Suppose uniform induced `m`-submatrices satisfied the mesoscopic recovery
claim: in the iterated limit `N->infinity` followed by `m->infinity`, their
directed one-profiles approached that of `T` and their normalized operator
norms were bounded by one finite `D`.  Directed one-profile continuity would
then give

```math
\Phi(T_{A_N[S_m]})\le1+\eta+o_{\mathbb P}(1).             \tag{MS.9}
```

But (MS.5)--(MS.7) imply that the same random variables exceed any fixed
number below `1.063846...` with probability tending to one.  Choose
`eta<0.063846...`; (MS.9) is a contradiction.

Therefore:

```math
\boxed{
\text{Uniform mesoscopic induced sampling cannot recover an extremal
bounded-operator action object.}}                         \tag{MS.10}
```

This rules out sampling plus correction whenever the correction is
`o(m^(3/2))` in Boolean/bilinear norm: such a correction cannot erase the
fixed iid gap in (MS.7).  It does not rule out a deliberately optimized
subset or a leading-scale correlated correction, but either requires a new
global selection mechanism rather than raw sampling.
