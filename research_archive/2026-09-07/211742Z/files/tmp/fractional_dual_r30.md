# Wave 30 Route 1: mesoscopic common-coset collision

Status: the finite-game reduction and the signature-conditioned planted-hash
theorem below are proved.  The final random-label construction is an abstract
mechanism wall: its prescribed child labels need not be realizable as the
ground states of the principal submatrices of one exact minimizer.

Put `L=n^(3/4-c)` and use `k=Theta(L/log n)` balanced blocks.  Let
`mathscr A_C` be the eligible row-good block cosets in (10.873), let
`h_w=max_a w(H_a)`, and let

```math
J_r=\{(S_1,\ldots,S_r):\text{some }a\in\mathscr A_C
\text{ has }S_1,\ldots,S_r\in H_a\}.
```

## 1. Exact mesoscopic collision sandwich

If `M=|mathscr A_C|`, then, pointwise for every selector law `w`,

```math
h_w^r\le \Pr_{w^{\otimes r}}(J_r)\le M h_w^r.                 \tag{A}
```

The lower bound retains one coset attaining `h_w`; the upper bound is the
union bound over cosets.  Since `delta_*=inf_w h_w`, taking infima and
`r`-th roots gives

```math
M^{-1/r}\inf_w\Pr_{w^{\otimes r}}(J_r)^{1/r}
\le\delta_*
\le\inf_w\Pr_{w^{\otimes r}}(J_r)^{1/r}.                    \tag{B}
```

There are at most `k^n` labelled block maps and at most `2^n` diagonal
translates for each map, so (also after balance and row pruning)

```math
M\le (2k)^n.                                                   \tag{C}
```

Take

```math
r=\left\lceil\frac{n\log(2k)}L\right\rceil
 =\Theta(n^{1/4+c}\log n).
```

Then `M^(1/r)<=e^L`.  Consequently the following statement is equivalent,
up to changing the constant in the exponent, to (10.875):

> **Mesoscopic common-coset lemma.**  Uniformly for every selector law `w`,
> `Pr_{w^r}(J_r)>=exp(-K r L)`.

Indeed this lemma and (B) give `delta_*>=exp(-(K+1)L)`, hence
`tau_*<=exp((K+1)L)`.  Conversely (10.875) makes `h_w>=exp(-K L)`
for every `w`, and the left inequality in (A) gives the mesoscopic lemma.
A law with `Pr(J_r)=exp(-omega(rL))` therefore falsifies (10.875) for the
fixed candidate family, because the right inequality in (B) gives
`delta_*=exp(-omega(L))`.

This is stronger than self-incidence in a useful way.  Self-incidence asks
only that one sampled selector plant its own coset.  Here a single row-good
coset must simultaneously hit a mesoscopic sample of
`Theta(n^(1/4+c) log n)` selectors, with only `exp(O(rL))` loss.

## 2. A checkable planted-completion condition for the common event

The common event has a sufficient formulation that exposes exactly what a
joint completion theorem must prove.  Fix selectors `S_1,...,S_r` and full
oriented completions `(sigma_ell,x^(ell))`.  Use `x^(1)` as the base word,
put `D=diag(x^(1))`, and partition the coordinates into their
relative-signature cells

```math
C_j=\left\{i:(x_i^{(2)}x_i^{(1)},\ldots,
x_i^{(r)}x_i^{(1)})=\eta_j\right\},\qquad j=1,\ldots,J.
```

Write `N_j=|C_j|`, `mu_j=AD 1_(C_j)`, and define

```math
C_{\rm sig}=\max_{\epsilon\in\{\pm1\}^J}
 \left\|\sum_j\epsilon_j\mu_j\right\|_2^2.                 \tag{D}
```

This is exactly the maximum row square of the coarsest signature coset.  It
is also a necessary lower bound for every containing block refinement,
because every such refinement contains all coarsest-coset words.

Here is a sharp-up-to-slack planted-refinement theorem.  Fix a maximum block
size `b>=2` and put

```math
k_j=\begin{cases}1,&N_j\le b,\\
\lceil2N_j/b\rceil,&N_j>b,\end{cases}
\qquad K_0=\sum_jk_j\le J+2n/b.                              \tag{E}
```

If `b>=9 log(4K_0)`, there is a partition refining the signature cells,
with at most `K_0` nonempty blocks and maximum block size `b`, such that

```math
\boxed{
\max_z\|ADPz\|_2^2
\le\left[\sqrt{C_{\rm sig}}+\sqrt{K_0}
 \left\{\sqrt{2\nu u}+\frac23\sqrt{n-1}\,u\right\}\right]^2,}
                                                                    \tag{F}
```

where

```math
\nu=\max\{\|A\|_{op}^2,b(n-1)/2\},
\qquad u=\log(4(n+K_0)).
```

For the proof, hash every large cell independently into its `k_j` private
labels and leave a small cell as one block.  With `g(i)` the private label,
local all-ones vectors `1_j`, and `v_i=x_i^(1)Ae_i`, decompose

```math
ADP=M_0+Z,
\quad M_0[:,(j,a)]=\mu_j/k_j,
\quad Z=\sum_i v_i(e_{g(i)}-k_j^{-1}\mathbf1_j)^T
\quad(i\in C_j).
```

The centered summands satisfy

```math
\|Z_i\|_{op}\le\sqrt{n-1},\qquad
\sum_i\mathbb E Z_iZ_i^T\preceq A^2,
```

```math
\left\|\sum_i\mathbb E Z_i^TZ_i\right\|_{op}
\le (n-1)\max_{j:k_j\ge2}N_j/k_j\le b(n-1)/2.
```

Rectangular matrix Bernstein bounds `||Z||_op` by the braces in (F),
except on probability at most `1/4`.  For a large cell,
`N_j/k_j in (b/3,b/2]`; binomial Chernoff and a union bound put every load
at most `b`, except on probability at most
`K_0 exp(-b/9)<=1/4`.  Hence the two events intersect.

For every block sign `z`,

```math
M_0z=\sum_j a_j\mu_j,\qquad
a_j=k_j^{-1}\sum_{a=1}^{k_j}z_{j,a}\in[-1,1].
```

The norm is convex on this cube, so its maximum is attained at a vertex and
is at most `sqrt(C_sig)`; also `||Zz||<=sqrt(K_0)||Z||_op`.  This proves
(F).  In expanded form its right side is at most

```math
2C_{\rm sig}+8K_0\nu u+\frac{16}{9}K_0(n-1)u^2.             \tag{G}
```

At the target scales, take
`b=Theta(n^(1/4+c) log n)` and assume

```math
J=O(L/\log n),\qquad C_{\rm sig}=O(n^{9/4-c}).               \tag{H}
```

Then `K_0=O(L/log n)`.  Exact minimality gives
`||A||_op^2=O(n^(3/2))`, and all Bernstein terms in (G) are
`O(n^(9/4-c))`.  Thus (H) produces a containing row-good block coset at
the required entropy scale.

This improves the initially tempting allocation-field condition

```math
k\sum_j\frac{\|\mu_j\|_2^2}{k_j}=O(n^{9/4-c}).              \tag{I}
```

Indeed weighted Cauchy--Schwarz gives
`C_sig<=k sum_j ||mu_j||^2/k_j`, so (I) is stronger than the necessary
coarsest-coset condition.  Conversely the crude implication from `C_sig`
to the allocation expression can lose a factor `k`.  The vertex argument
above removes that loss, making `C_sig` the sharper hypothesis.

It remains to connect favorable completions to hits.  Put

```math
\delta_S(\sigma,x)=Q(A[S])-\sigma x_S^TA[S]x_S.
```

Direct expansion and `sigma x^TAx<=q_n` give

```math
\sigma x^TH_Sx\ge Y_A(S)+B_{n,m}-\delta_S(\sigma,x).         \tag{J}
```

Thus `delta_S<=B_(n,m)+t` makes every block coset containing `x` hit `S`
at tolerance `t`; exact child completions have `delta_S=0`.

It follows that the following is an exact structural sufficient lemma for
(10.875), and is the concrete endpoint of the planted-hash route:

> **Mesoscopic planted-completion lemma.**  For every selector law `w`, an
> `r`-sample has probability at least `exp(-K rL)` of admitting completions
> satisfying `delta_(S_ell)<=B_(n,m)+t` whose relative-signature family
> satisfies (H).

Unlike the uniform-mean identities (10.823)--(10.824), this is a nonlinear
high-order incidence statement.  It also does not use the low-information
channel extraction (10.863): no selector-to-output mutual information is
estimated.  Its missing content is genuine mesoscopic coherence of favorable
completions and their coarsest-coset row field.

## 3. Rigorous wall for local planting alone

The preceding coherence cannot follow from the bare fact that every selector
has a favorable child completion.  This can be made rigorous even if every
row cap is deleted.

Fix a density `p in (0,1)`, let `m=floor(pn)`, and independently prescribe
for every `m`-selector `S` a uniform projective label
`y_S in {+1,-1}^S/{+/-}`.  In this abstract incidence model, declare a full
spin favorable for `S` exactly when its restriction is `+/- y_S`; hence all
`2^(n-m)` outside completions are available, just as in (10.825).

For any fixed block coset with at most `k=o(n)` blocks, at most `2^k`
projective restrictions occur on `S`.  Therefore, over the random labels,
its probability of hitting a fixed `S` is at most

```math
p_0=2^{k-m+1}=\exp(-p(\log2)n+o(n)).                          \tag{K}
```

The hit indicators are independent over selectors.  There are at most
`(2k)^n=exp(O(n log n))` block cosets altogether.  Choose

```math
0<\alpha<\min\{p\log2,h(p)\},\qquad \delta=e^{-\alpha n},
```

where `h` is binary entropy.  For `N=binom(n,m)`, the binomial Chernoff bound
and a union bound give

```math
\Pr\{\text{some coset hits at least }\delta N\text{ selectors}\}
\le (2k)^n
\left(\frac{e p_0}{\delta}\right)^{\delta N}=o(1).           \tag{L}
```

Indeed `delta/p_0=exp((p log2-alpha)n+o(n))`, while
`delta N=exp((h(p)-alpha)n+o(n))`, so the negative logarithm from the last
factor is exponential in `n` and dominates `O(n log n)`.

Thus some deterministic system of local favorable labels has, for the
uniform selector law and its abstract completion-hit sets `\widetilde H_a`,

```math
\max_a w(\widetilde H_a)\le e^{-\alpha n}=e^{-\omega(L)}.     \tag{M}
```

For the mesoscopic `r` above, the analogue of (A) also gives
`Pr(\widetilde J_r)<=M e^{-alpha n r}=exp(-omega(rL))`.  This wall permits every
outside completion, every partition, every diagonal translate, and every
coset as row-good.  Hence neither jointly randomizing the locally planted
completion nor hashing it can by itself yield a stronger-than-self-incidence
bound.  A proof must import a minimizer-specific correlation theorem strong
enough to establish the mesoscopic planted-completion lemma (or make common
hits occur without containing the chosen ground completions).

The wall is deliberately scoped: independently prescribed labels need not be
the exact child grounds of one signing.  It therefore does not falsify
(10.875); it identifies the precise additional structure that an exact
minimizer theorem must supply.
