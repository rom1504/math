# PSD gluing fibres have bounded fixed-scale compatibility entropy

**Status.** Rigorous task-local parametrization, uniform upper cover, and
Boolean-realizable lower packing.  Given fixed marginal Gram--Rayleigh
sectors, the exact gluing fibre can have continuously many parameters and
growing rank.  Nevertheless, in the collective Boolean quadratic metric its
fixed-accuracy covering complexity is independent of the number of ports.

The proof identifies the exact compatibility object: a contraction between
the support Hilbert spaces of each pair of marginal PSD sectors.  Spectral
truncation shows that only finitely many macroscopic support modes survive at
fixed response accuracy.

## 1. Exact contraction parametrization

Let `K_L succeq 0` and `K_R succeq 0` have sizes `p` and `q`.  A cross block
`X in R^(p times q)` is compatible if

```math
\widehat K_X=
\begin{pmatrix}K_L&X\\X^T&K_R\end{pmatrix}\succeq0.       \tag{PF.1}
```

### Theorem PF.1 (the compatibility fibre is a contraction ball)

Equation PF.1 holds if and only if

```math
X=K_L^{1/2}CK_R^{1/2}                                \tag{PF.2}
```

for a contraction `C` from `supp(K_R)` to `supp(K_L)`.  Equivalently, the
compressed operator

```math
C_X=K_L^{dagger/2}XK_R^{dagger/2}                   \tag{PF.3}
```

has norm at most one and `X` has the corresponding row and column range
conditions.  The compressed `C_X` is unique; extensions away from the two
supports do not change `X`.

#### Proof

If PF.2 holds, then

```math
\widehat K_X=
\begin{pmatrix}K_L^{1/2}&0\\0&K_R^{1/2}\end{pmatrix}
\begin{pmatrix}I&C\\C^T&I\end{pmatrix}
\begin{pmatrix}K_L^{1/2}&0\\0&K_R^{1/2}\end{pmatrix}.    \tag{PF.4}
```

The middle matrix is PSD exactly when `||C||<=1`, by its Schur complement.

Conversely, factor PF.1 as a joint Gram matrix.  If `U,V` are the two
families of Gram vectors, their polar decompositions have the form

```math
U=Q_LK_L^{1/2},
\qquad V=Q_RK_R^{1/2},                               \tag{PF.5}
```

where `Q_L,Q_R` are isometries on the respective supports into the common
Gram space.  Hence

```math
X=U^TV=K_L^{1/2}(Q_L^TQ_R)K_R^{1/2},                \tag{PF.6}
```

and `Q_L^TQ_R` is a contraction.  Pseudoinversion gives PF.3 and the range
conditions. `square`

For a Gram--Rayleigh pair there are two independent sectors
`K^+-=(G+-R)/2`.  Fixing left and right marginals therefore leaves a product
of two contraction fibres, one for `K^+` and one for `K^-`.

## 2. The inherited collective response metric

Let `P=p+q`.  For two compatible choices `X=(X^+,X^-)` and
`Y=(Y^+,Y^-)`, the collective metric CP.1 on the joined pair reduces exactly
to

```math
d_f(X,Y)
={4\over P^2}\max_{epsilon\in\{+-1\}^p,
                         eta\in\{+-1\}^q}
 \max_{tau\in\{+,-\}}
 |epsilon^T(X^tau-Y^tau)eta|.                       \tag{PF.7}
```

Indeed the marginal blocks cancel, a full-sector quadratic form contributes
twice its cross bilinear form, and CP.2 contributes the second factor two.

The contraction parametrization also gives the immediate diameter bound

```math
d_f(X,Y)
\le {8\over P^2}\max_tau
 \sqrt{A_L^tau A_R^tau},                            \tag{PF.8}
```

where

```math
A_L^tau=\max_epsilon epsilon^TK_L^tauepsilon,
\qquad
A_R^tau=\max_eta eta^TK_R^taueta.                   \tag{PF.9}
```

This already shows the tradeoff: a high-rank isotropic marginal has a small
fibre diameter, while a coherent rank-one marginal can have constant
diameter but only scalar compatibility.

## 3. A uniform fixed-accuracy cover

Assume the marginal sectors are Gram--Rayleigh admissible, so each is PSD,
each diagonal lies in `[0,1]`, and the two sector diagonals sum to one.

### Theorem PF.2 (fixed-scale compatibility costs no port-rate)

There is an absolute `C` such that, for every `0<zeta<1/4`, every choice of
marginals, and every `p,q`,

```math
\log Cov_zeta(\mathfrak F,d_f)
\le C zeta^{-4}\log(C/zeta).                       \tag{PF.10}
```

In particular, at fixed collective-response accuracy the compatibility
fibre costs `O_zeta(1)` states as `p+q` tends to infinity, not
`O_zeta(p+q)` or worse.

#### Proof

Put

```math
theta=(zeta/16)^2.                                  \tag{PF.11}
```

For every marginal sector retain only eigenvectors with eigenvalue at least
`theta P`.  Since every sector has trace at most its number of ports, each
retained space has rank at most

```math
{tr K\over theta P}\le theta^{-1}.                 \tag{PF.12}
```

Write `K^(1/2)=L_top+L_tail` for the spectral split.  If
`X=K_L^(1/2)CK_R^(1/2)`, replace it by

```math
X_0=L_top C R_top.                                  \tag{PF.13}
```

This remains compatible because the compressed contraction
`P_top C P_top` has norm at most one.  For Boolean `epsilon,eta`,

```math
\begin{aligned}
|epsilon^T(X-X_0)eta|
&\le ||L_tail epsilon||\,||K_R^{1/2}eta||
   +||L_top epsilon||\,||R_tail eta||\\
&\le q\sqrt{theta Pp}+p\sqrt{theta Pq}
\le2\sqrt{theta}P^2.                               \tag{PF.14}
\end{aligned}
```

Here `epsilon^TKepsilon<=p^2` follows from PSD Cauchy--Schwarz and
`K_ii<=1`; the tail operator norm is at most `theta P`.  Equations
PF.7 and PF.11 show

```math
d_f(X,X_0)<=8\sqrt{theta}=zeta/2.                  \tag{PF.15}
```

It remains to cover the top-space contractions.  Both their row and column
dimensions are at most `theta^{-1}`.  If two contractions differ in
operator norm by at most `xi`, then

```math
\begin{aligned}
|epsilon^TL_top(C-C')R_top eta|
&\le xi\sqrt{epsilon^TK_Lepsilon}
          \sqrt{eta^TK_Reta}\\
&\le xi pq,
\end{aligned}                                      \tag{PF.16}
```

so their fibre distance is at most `4xi pq/P^2<=xi`.  Take
`xi=zeta/2`.

An `r_L by r_R` operator-norm unit ball lies in a Frobenius ball of radius
`sqrt(min(r_L,r_R))`.  A Euclidean volume net in Frobenius norm therefore
has size at most

```math
\left(1+{2\sqrt{min(r_L,r_R)}\over xi}\right)^{r_Lr_R}.
                                                               \tag{PF.17}
```

Frobenius accuracy implies operator accuracy.  Using PF.12 in PF.17, and
covering the two sectors independently, gives PF.10. `square`

This theorem is a statement about approximate compatibility only.  An exact
contraction matrix can require quadratically many real parameters.

## 4. Boolean-realizable lower packing

The independence from port count is not the same as a one-state theorem:
the dependence on accuracy is necessarily nontrivial.

Take `r` orthogonal Boolean top eigenvectors of a tensor regular-Hadamard
matrix.  Duplicate each vector `L` times on both shores, so each shore has
`p=rL` labelled ports.  The fixed marginal positive-sector Gram matrix is

```math
K_L^+=K_R^+=I_r\otimes J_L,
\qquad K_L^-=K_R^-=0.                               \tag{PF.18}
```

For a permutation `pi in S_r`, order the right-hand mode frames by `pi`.
The cross block is

```math
X_pi^+=P_pi\otimes J_L,
\qquad X_pi^-=0.                                   \tag{PF.19}
```

Every joined pair is the exact Gram--Rayleigh pair of Boolean top
eigenvectors, and the marginals do not depend on `pi`.

### Theorem PF.3 (signed-frame lower entropy)

For every sufficiently small `zeta`, there are fixed Boolean-realizable
marginals whose compatibility fibre contains a `zeta`-packing of size

```math
\exp\left(c zeta^{-1}\log(1/zeta)\right).          \tag{PF.20}
```

The construction remains valid for arbitrarily large port counts by
increasing the duplication factor `L`.

#### Proof

For two permutations, let `d_H(pi,sigma)` be the number of positions on
which they differ.  Put `rho=sigma^{-1}pi`.  For a right sign vector `y`,

```math
||(P_pi-P_sigma)y||_1
```

is twice the number of edges of the cycle permutation `rho` cut by `y`.
Every nontrivial even cycle can be cut on every edge; an odd cycle of length
`ell>=3` can be cut on `ell-1>=2ell/3` edges.  Hence

```math
||P_pi-P_sigma||_(infinity->1)
\ge {4\over3}d_H(pi,sigma).                        \tag{PF.21}
```

Duplication multiplies this norm by `L^2`.  Since the joined system has
`2rL` ports, PF.7 gives

```math
d_f(X_pi,X_sigma)
={||P_pi-P_sigma||_(infinity->1)\over r^2}
\ge {4d_H(pi,sigma)\over3r^2}.                    \tag{PF.22}
```

A greedy permutation code with relative Hamming distance `1/2` has size

```math
\exp\left({1\over2}r\log r-O(r)\right).           \tag{PF.23}
```

Indeed, the number of permutations differing from a fixed one on fewer
than `r/2` positions is at most
`sum_(k<r/2) binom(r,k)k! <= exp((r/2)log r+O(r))`, whereas
`r!=exp(rlog r-O(r))`.  Choose `r=floor(c_0/zeta)` with a sufficiently small
absolute `c_0`; PF.22 then gives pairwise distance at least `zeta`, and
PF.23 proves PF.20. `square`

A simpler diagonal-sign subfamily has

```math
X_s^+=D_s\otimes J_L,
\qquad
d_f(X_s,X_t)={2d_H(s,t)\over r^2},                 \tag{PF.24}
```

and already gives `log Pack_zeta=Omega(zeta^{-1})`.  These are scalar
alignment coordinates in a fixed labelled frame; permutations add
frame-matching information.  Section 5 distinguishes these frame signs from
the separate marginal-sector antipode in OC.2.

Together PF.2--PF.3 classify the port-count regime:

```math
\boxed{\text{fixed-accuracy PSD compatibility has }O_zeta(1)
       \text{ information cost as }p+q\to\infty.}  \tag{PF.25}
```

The current accuracy exponents are not claimed sharp.

## 5. Exact boundary with the regular-Hadamard orientation/cycle fibre

There are two different `Z_2` operations, and identifying them would be an
overclaim.  A sign change of a chosen Gram **factor** is gauge.  The OC.2
relative antipode instead reverses the sign of a marginal energy channel;
in Gram--Rayleigh coordinates it swaps the two spectral sectors:

```math
R\longmapsto-R,
\qquad
(K^+,K^-)\longleftrightarrow(K^-,K^+).             \tag{PF.26}
```

Thus, if an absolute closed marginal retains only the unordered sector pair
`{K^+,K^-}`, gluing `s` pieces requires choosing an ordering for every pair.
A simultaneous swap is the one global output antipode, leaving exactly
`2^(s-1)` relative orderings.  This is the precise PSD-sector realization of
the `s-1` relative marginal antipodes in OC.12.  PF.1 fixes ordered sectors,
so those bits must be adjoined before applying its contraction fibre.

The OC **cycle** variable has a different status.  To see the distinction
exactly, give every marginal a one-dimensional Gram support and suppose all
edge compatibilities are isometric.  Numerically one can write

```math
C_e\in O(1)=\{+-1\},
\qquad
C_(uv)\longmapsto d_uC_(uv)d_v                    \tag{PF.27}
```

under a local factor-frame change.  This is the same switching action as on
signed edges.  But a family of such edge blocks comes from **one global PSD
Gram matrix** only when there are global factor signs `q_v` with

```math
C_(uv)=q_uq_v.                                     \tag{PF.28}
```

Consequently every cycle must satisfy

```math
\prod_{e\in cycle}C_e=+1.                          \tag{PF.29}
```

Conversely, on a connected graph, the positive-cycle condition lets one
define the `q_v` along paths and gives a global rank-one Gram realization.
Therefore a negative OC cycle holonomy is **not** an element of the PSD
Gram compatibility fibre.  It is coefficient-side interaction information
which survives precisely because regular-Hadamard bridge signs need not be
mutual inner products of one family of marginal Gram factors.

The same distinction persists in higher rank.  If full-dimensional
isometric edge alignments arise from global frames `Q_v`, then

```math
C_(uv)=Q_u^TQ_v
```

and every ordered cycle product is the identity (up to the harmless choice
of base-frame coordinates).  A nontrivial orthogonal holonomy obstructs one
global Gram realization; it is not a free PF.1 compatibility coordinate.

The exact comparison is therefore:

1. **Relative OC antipodes:** unordered-to-ordered sector choices, PF.26,
   giving `s-1` bits modulo the global swap.
2. **PSD two-shore alignment:** the contraction `C` in PF.2, including the
   labelled sign and permutation families PF.19/PF.24.
3. **OC cycle holonomy:** a coefficient-side fibre absent from a globally
   PSD Gram gluing whenever its cycle product is nontrivial.

This is a structural boundary, not a metaphorical analogy.  PF.2 proves
bounded fixed-scale entropy for item 2 only.  It does not compress away the
OC cycle bits of item 3.  OC.3--OC.4 use fixed graphs where either a relative
sector sign or a coefficient cycle sign has a leading cap witness; PF.2
instead describes the approximate dilution of many PSD-realizable alignment
modes under one globally normalized collective metric.

## 6. Research consequence

The PSD compatibility fibre is not the source of a port-extensive
fixed-accuracy obstruction under `d_q`.  Marginal spectral mass can support
only `O_zeta(1)` macroscopic alignment modes, and the remaining contraction
coordinates are collectively invisible at scale `zeta`.

The unresolved cost in a collective Gram carrier must therefore lie in one
of three places outside this theorem: covering the marginal PSD response
images themselves, realizing their covers by exact Boolean ports, or
passing the spherical carrier through its Boolean integrality gap.
