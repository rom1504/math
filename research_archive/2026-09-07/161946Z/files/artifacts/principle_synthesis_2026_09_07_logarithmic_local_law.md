# Low-cap switching laws are locally iid through logarithmic vertex order

2026-09-07. New synthesis-track theorem. The elementary derivation below is
complete; independent audit is requested. It constructs an actual law supported
on the switching/permutation orbit of ANY prescribed signing. In particular,
for an exact minimizer every output remains an exact minimizer. It does not
construct an optimizer at a new order and does not prove convergence.

Throughout, `H_A(x)=sum_(i<j) A_ij x_i x_j`, `Q(A)=max_x |H_A(x)|`, and
`M_n=min_A Q(A)`. Relative entropy and logarithms are in natural units; total
variation is one half of the `l^1` distance.

## 1. Finite theorem, with explicit error

Let A be any order-n hollow symmetric full signing and put B=A+I, so B is a
full symmetric sign matrix. Let L_(A,k) be the law on the edges of K_k obtained
by choosing k distinct vertices of A uniformly in order and switching them by
k independent fair signs. Let U_k be the iid fair-sign law on K_k. Set

```math
\tau_4(B)=\frac{\operatorname{tr}(B^4)}{n^4},\qquad
\Lambda_k=2^{k-1}-k\quad(k\ge1).
```

Then, for `1<=k<=n`,

```math
\boxed{\quad
d_{\rm TV}(L_{A,k},U_k)
\le \sqrt{\frac{\Lambda_k\tau_4(B)}2}
     +\frac{k(k-1)}{2n}.
\quad}                                                     \tag{1}
```

The expression `Lambda_k` is nonnegative for every k>=1. The square-root
term comes from an exact relative-entropy bound for sampling WITH replacement:

```math
D(\widetilde L_{B,k}\Vert U_k)
\le \Lambda_k\tau_4(B).                                  \tag{2}
```

There is no claim that the without-replacement law satisfies (2); its coupling
to the with-replacement law is paid separately in (1).

For any fixed C, the established bounded-entry fourth-moment theorem gives

```math
Q(A)\le Cn^{3/2}
\quad\Longrightarrow\quad
\operatorname{tr}(B^4)
\le K_G^2\bigl(4Cn^{3/2}+n\bigr)^2.                       \tag{3}
```

Indeed the full-matrix bilinear norm is at most `4Q(A)+n`, the maximum entry
of B is one, and `||B||_(S4)^2<=K_G beta(B)`. Therefore for every fixed
`delta>0`, uniformly over ALL such actual signings and
`k<=floor((1-delta)log_2 n)`,

```math
d_{\rm TV}(L_{A,k},U_k)
=O_C(n^{-\delta/2})+O((\log n)^2/n)=o(1).                 \tag{4}
```

The original minimum has a uniform upper cap bound, so (4) applies uniformly
to every exact minimizer and every additive near-minimizer. Selection of a
particular representative is not needed.

## 2. One-row Fourier calculation

For the with-replacement model, let V_1,...,V_k be independent uniform
vertices of [n] and let S_1,...,S_k be independent uniform signs, independent
of the vertices. Set

```math
Z_{ij}=S_iS_j B_{V_i,V_j},\qquad i<j.
```

Reveal the edge rows successively. Conditional on all previously exposed
latent variables `V_1,...,V_t,S_1,...,S_t`, the next row is

```math
R_j=S_{t+1}S_j B_{V_{t+1},V_j},\qquad 1\le j\le t.
```

For a subset F of [t], its Fourier coefficient is the conditional mean of
`prod_(j in F) R_j`. It vanishes whenever |F| is odd. For nonempty even F it is

```math
\widehat p(F)=\Bigl(\prod_{j\in F}S_j\Bigr)
        \frac1n\sum_{v=1}^n\prod_{j\in F}B_{v,V_j}.
```

Fourier Parseval for the density relative to the uniform t-bit law gives

```math
\chi^2(p\Vert U_t)=\sum_{F\ne\varnothing}\widehat p(F)^2.
```

Put `r_vw=(B B^T)_vw/n`; then |r_vw|<=1. Averaging over the independent old
vertex labels gives the exact identity

```math
\begin{split}
\mathbb E\chi^2(p\Vert U_t)
 &=\frac1{n^2}\sum_{v,w}
   \left[\frac{(1+r_{vw})^t+(1-r_{vw})^t}{2}-1\right]\\
 &=\frac1{n^2}\sum_{v,w}
      \sum_{\substack{2\le j\le t\\j\ {\rm even}}}
              {t\choose j}r_{vw}^j.
\end{split}                                               \tag{5}
```

Every surviving power is even. Hence for t>=1 the last line is at most

```math
\bigl(2^{t-1}-1\bigr)\frac1{n^2}\sum_{v,w}r_{vw}^2
=\bigl(2^{t-1}-1\bigr)\frac{\operatorname{tr}(B^4)}{n^4}.
                                                               \tag{6}
```

For t=1 both sides are zero. At t=0 the row is empty and its entropy cost is
zero; (6) is not invoked there. Symmetry of B is used only to identify
`||BB^T||_F^2` with `tr(B^4)`.

## 3. Conditioning correctly and summing entropy

Let E_t denote the already observed edge array on the first t vertices and
let G_t denote their full latent labels and switches. The array E_t is a
deterministic function of G_t. Therefore conditional entropy decreases when
the full latent past is revealed:

```math
H(R\mid E_t)\ge H(R\mid G_t).
```

Equivalently, convexity of relative entropy gives

```math
\mathbb E D(\mathcal L(R\mid E_t)\Vert U_t)
\le \mathbb E D(\mathcal L(R\mid G_t)\Vert U_t).
```

The latent labels remain iid when the right side is averaged unconditionally.
No assertion that they remain iid conditional on E_t is made or needed.
Using `D(p||U_t)<=log(1+chi^2(p||U_t))<=chi^2(p||U_t)`, the chain rule and
(6) give

```math
\begin{split}
D(\widetilde L_{B,k}\Vert U_k)
&\le\tau_4(B)\sum_{t=1}^{k-1}(2^{t-1}-1)\\
&=\bigl(2^{k-1}-k\bigr)\tau_4(B).
\end{split}
```

This proves (2). Pinsker's inequality gives its square-root total-variation
bound. Finally the with-replacement sample has a repeated vertex with
probability at most `binom(k,2)/n`; conditional on no repetition it is exactly
the ordered distinct-vertex model. On that event the diagonal B_ii is never
queried. Coupling proves (1).

## 4. Exact cap preservation, pairwise independence, and loop dependence

Choose a uniform permutation pi of all n vertices and independent vertex
switches S. The full output

```math
C_{ij}=S_iS_j A_{\pi(i),\pi(j)},\qquad i\ne j,
```

has Q(C)=Q(A) for EVERY output, not merely on average. Its induced k-vertex
marginal is L_(A,k). If A is exact minimizing, the full law is supported
entirely on exact minimizers.

Distinct edges have zero covariance: their product involves at least one
vertex switch to an odd power. Each individual edge is fair. More generally,
the edge-product expectation vanishes unless its support is Eulerian. For
Eulerian support the switches cancel and the expectation retains the signed
subgraph statistic of A averaged over distinct labels. Thus identity edge
covariance and asymptotically iid logarithmic vertex marginals coexist with
an exactly preserved non-iid global Boolean cap.

In particular, covariance-only or fixed-local-law information does not
justify replacing arbitrary low-cap gauge laws by independent disorder. This
does NOT contradict the Gaussian-sign universality theorem: that theorem
assumes a specified bounded-spectrum nondegenerate Gaussian-sign law, not an
arbitrary law with its covariance. Nor does this show that compressed local
information is required for an eventual convergence proof.

## 5. A matching logarithmic-order upper boundary for this sampler

For every fixed base signing A, the support of L_(A,k) has size at most

```math
(n)_k 2^k\le n^k2^k.
```

The iid law has `2^(k(k-1)/2)` equally likely arrays. If T is the support of
L_(A,k), then

```math
\boxed{\quad
d_{\rm TV}(L_{A,k},U_k)
\ge 1-2^{\,k(\log_2n+1)-k(k-1)/2}.
\quad}                                                     \tag{7}
```

For `k>=ceil((2+delta)log_2 n)` and fixed delta>0, the right side tends to one.
This support bound needs neither low cap nor spectral control and applies
equally to the orbit of an exact minimizer. It is a boundary for this fixed-
seed switching/permutation sampler, not for mixtures over many distinct seeds.

Equations (4) and (7) locate the indistinguishability threshold at logarithmic
vertex order, up to a factor of two. Nothing here resolves the interval between
`(1-delta)log_2n` and `(2+delta)log_2n`.

## 6. Consequence for actual principal restrictions

Let `k_n` tend to infinity with `k_n<=(1-delta)log_2 n`, and choose a uniform
k_n-vertex subset I of ANY bounded-cap A_n. Although (1) concerns the switched
law, cap is exactly switching-invariant. The established iid sequential-greedy
construction has

```math
\frac{Q(G_k)}{k^{3/2}}
\ge c_{\rm gr}-o_{\mathbb P}(1),\qquad
c_{\rm gr}=\frac23\sqrt{\frac2\pi}=0.5319230405\ldots .
```

For completeness, greedily assigning vertex t makes its contribution the
absolute value of a sum of t-1 fresh iid signs. Those contributions are
independent, have total mean `(c_gr+o(1))k^(3/2)`, and total variance at most
k(k-1)/2. Chebyshev's inequality proves the displayed lower bound. Thus (1)
gives, uniformly over the actual original input class,

```math
\Pr\!\left\{\frac{Q(A_n[I])}{k_n^{3/2}}
                   \ge c_{\rm gr}-\eta\right\}\longrightarrow1
\quad\text{for every fixed }\eta>0.                       \tag{8}
```

For exact minimizing inputs this has a fixed normalized loss relative to the
current all-order upper bound `U=0.49360809358874865...`; the difference is
`c_gr-U=0.0383149469...`. This improves the archived random-restriction scale
from order `sqrt(log n)` to order `log n`. It does not exclude exceptionally
good selected subsets, positive retention ratios, a different sampling law,
or global changes to internal edges.

## 7. Dependencies and status

The new steps are the row Fourier identity (5), latent-past chain-rule bound,
and their quantitative logarithmic-order combination. The fourth-moment input
is reconstructed in
`flatify_construct_2026_09_07_spectral_fourth_moment.md`; its maximum-entry
factor and full-vs-hollow normalization are retained in (3).

Earlier fixed-profile and induced-sampling artifacts were inspected, including
`phase2_profile_state_main_audit.md`, `ar_director_synthesis.md`, and
`transfer_adversary_random_restriction_2026_09_06.md`. They do not supply the
present growing-k total-variation estimate. No literature novelty is asserted.

The associated exact finite checker is
`computations/principle_synthesis_2026_09_07_logarithmic_local_check.py`.
It tests the finite identity and normalization only; the proof does not depend
on finite tests or on small-order optimum certificates.
