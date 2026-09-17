# Growing dependent blocks: an all-offset Gaussian replacement theorem

2026-09-17. **Proof reconstructed below; independent audits of the
independent-block statement, backward-conditional extension, and
linear-block counterexample all passed.** This extends the fixed-size block argument in
[the archived frame audit](principle_construct_2026_09_07_fixed_frame_gaussian_obstruction.md),
not a claim to have invented the Lindeberg method. External novelty has
not been established. Relevant primary background is Chatterjee,
[*A simple invariance theorem*](https://arxiv.org/abs/math/0508213) and
[*A generalization of the Lindeberg principle*](https://arxiv.org/abs/math/0508519).

## 1. General finite theorem

Let `T={t_1,...,t_K}` be a finite set of vectors with coordinates in
`[-1,1]`, partitioned into blocks of sizes `b_l`. Let `h_j` be arbitrary
deterministic offsets. Suppose the centered blocks X_l are independent
and are L-subGaussian in the Euclidean metric, with a common finite L.
Let G_l be independent centered Gaussian blocks with the SAME covariance
as their respective X_l. Set `b=max b_l`, `S=sum b_l^(3/2)`, and
`Psi(X)=E max_j(h_j+<t_j,X>)`. For every `tau>0`,

```math
|\Psi(X)-\Psi(G)|
\le {2\log K\over\tau}
 +C L^3\tau^2 S\exp(4L^2\tau^2b).                 \tag{1}
```

Here C is absolute. The result permits arbitrary dependence WITHIN a
block; no coordinatewise independence or product Gibbs law is used.
The offsets are retained throughout, without their own norm bound.

**Proof of the decisive tilted-moment step.** Freeze every other block.
Write p_j for its normalized softmax weights and
`z_j=<t_(j,l),X_l>`, `bar z=sum p_j z_j`. Along interpolation parameter
`s in [0,1]`, the new weights are proportional to
`p_j exp(tau s z_j)`. Jensen gives their denominator at least
`exp(tau s bar z)`. Consequently

```math
\mathbb E\sum_j p_j(s)|z_j|^3
\le\sum_j p_j\mathbb E[|z_j|^3e^{\tau s(z_j-\bar z)}]
\le C L^3 b_l^{3/2}\exp(4L^2\tau^2b_l).             \tag{2}
```

For the last inequality, apply Cauchy--Schwarz. SubGaussian tails give
`E|z_j|^6 <= C L^6 b_l^3`; the coefficient vector of `z_j-bar z` has
norm at most `2sqrt(b_l)`. The square root of its exponential moment at
`2tau s` is at most `exp(4L^2 tau^2 b_l)`.

The third derivative, with respect to s, of the softmax divided by tau
is `tau^2` times the third Gibbs cumulant of z. Its absolute value is
at most `8tau^2 sum p_j(s)|z_j|^3`. Taylor expansion at zero has matching
first two EXPECTED terms for X_l and G_l. Integrating the third derivative
and using (2) bounds their replacement error by
`C L^3 tau^2 b_l^(3/2) exp(4L^2 tau^2 b_l)`.
The covariance of X_l is bounded by `L^2 I`, so G_l has the same
subGaussian bound. Sum the replacements and pay at most `logK/tau`
for each endpoint's maximum/softmax comparison. This proves (1).

Replacing a block by its worst possible l1 magnitude would cost `b_l^3`,
not `b_l^(3/2)`. The tilted-moment argument is the actual improvement;
it does not maximize over all witnesses before taking the block moment.

## 2. Dense signing and bridge normalization

For `d<=C_0 n^2` coordinates and at most `exp(C_1 n)` Boolean witnesses,
one has `S<=d sqrt(b)`. If `1<=b<=n`, choose
`tau=n^(-1/3)b^(-1/6)`. Then (1) gives

```math
|\Psi(X)-\Psi(G)|\le C_{L,C_0,C_1}\,n^{4/3}b^{1/6}.
                                                               \tag{3}
```

In particular `b=o(n)` gives `o(n^(3/2))` loss. If
`b<=n^(1-delta)`, the loss is `O(n^(3/2-delta/6))`.
The bound covers the **scalar absolute quadratic maximum**, by using
the two polarities of the usual Boolean edge-feature vectors. It also
covers a dense bridge with arbitrary fixed child-energy offsets.

If the X_l are centered full signs with covariance I, the Gaussian
reference is the ordinary independent-edge Gaussian model. A second
independent-sign comparison yields

```math
|\mathbb E Q(X)-\mathbb E Q(\varepsilon)|
\le C_L n^{4/3}b^{1/6}.                            \tag{4}
```

Thus bounded-size dependent constraints, and even uniformly subGaussian
independent blocks of size `o(n)`, cannot improve the random-sign leading
cap if they leave covariance unchanged. Conversely (3) is a positive
transfer theorem: a tractable covariance-matched Gaussian construction
with such a realizable full-sign law transfers its all-witness bound
at the displayed power-saving loss. The existence of that sign law and
the favorable Gaussian cap are separate obligations, not assumed solved.

At block size proportional to n the displayed error is leading order.
This is a LIMIT OF THE PROOF, not a counterexample to all larger-block
universality. Pairwise independence alone does not imply the uniform
subGaussian hypothesis.

## 3. Backward-conditional variant

Independence can be replaced by the following precise conditions. In a
fixed ordering, for every l and almost every future
`(X_(l+1),...,X_L)`, require the conditional law of X_l to have:

- mean zero;
- a deterministic covariance matrix Sigma_l, not depending on the future;
- the same Euclidean L-subGaussian MGF bound.

Compare with independent Gaussian blocks of covariances Sigma_l. In the
l-th telescoping term, earlier blocks are already independent Gaussians,
and later blocks retain their original marginal. Conditional on these,
the current block has exactly the stated three properties. Thus the
same Taylor cancellation and (2) apply, proving (1)--(4) unchanged.
This does NOT assume independence by another name: a uniform global
parity law, divided into blocks of at least five coordinates, supplies
a non-product example with deterministic conditional covariance I.
Its last block is unconstrained after future conditioning; earlier
conditional laws are mixtures of parity laws and remain uniformly
subGaussian. Precise parity examples are audited separately.

This conditional form does not hold automatically for a Gram--Schmidt
walk: its future-conditioned moments can differ from its forward
martingale moments. No GS Gaussian replacement is claimed on that basis.

## 4. Relationship to the other paper tracks

The bounded-subGaussian convex-order theorem retains arbitrary offsets
but costs a universal constant. Here offsets are also retained, while
many individually small blocks improve the comparison factor to one
with a quantitative additive error. The two statements have different
regimes; neither subsumes the other.

The information/retuning/dependence theorem shows that a successful
recovery law needs a macroscopic resource. Formula (4) shows that the
amount of dependence alone is not enough: many local parity constraints
can have extensive entropy cost while leaving the leading cap unchanged.
Their location and conditional structure matter.

No original bound or convergence result follows without a favorable
covariance/sign-realization construction. The theorem isolates a regime
where higher-than-pairwise block information can be discarded at the
required asymptotic accuracy, and gives its exact error exponent.

## 5. Sharpness at linear block size: actual sign rows

**Constructed and proved below; independent verification passed.**
There are centered sign vectors `X^(n) in {+-1}^n` such that

```math
\operatorname{Cov}(X^{(n)})=I_n,\qquad
\mathbb E e^{\langle t,X^{(n)}\rangle}
 \le e^{C\|t\|_2^2/2}
```

with an absolute C, but

```math
\liminf_n\left[
 \mathbb E\left|n^{-1/2}\sum_i\varepsilon_i\right|
 -\mathbb E\left|n^{-1/2}\sum_iX_i^{(n)}\right|\right]>0.
                                                               \tag{5}
```

Here n runs through all sufficiently large integers; a subsequence is
not required.

**Exact construction.** Under uniform signs put `Z=n^(-1/2)sum epsilon_i`.
Use three symmetric cells in |Z|:

```math
C_0=[0,1/10],\quad C_1=[9/10,11/10],\quad
C_2=[19/10,21/10].
```

Let `p_i=P(|Z| in C_i)`, `a_i=E[Z^2||Z| in C_i]`, and
`b_i=E[|Z|||Z| in C_i]`. For all large n each p_i is positive. Define
f to be zero off the three cells and to take values

```math
f_0=(a_2-a_1)/p_0,\quad f_1=-(a_2-a_0)/p_1,
\quad f_2=(a_1-a_0)/p_2.
```

Set `delta=p_1/[2(a_2-a_0)]` and change density by `w=1+delta f`.
Then `Ew=1`, `Ew Z^2=1`, and w is exactly `1/2` on the middle cell
and at least one elsewhere. The elementary central limit theorem and
bounded-cell moment convergence show that all p_i and a_i converge to
their corresponding standard Gaussian quantities. Thus delta tends to
a strictly positive number and `1/2<=w<=C` with an absolute finite C
for all sufficiently large n.

The new law is exchangeable and invariant under global sign reversal.
Its mean is zero. Since `Ew Z^2=1`, expanding the square of the spin sum
shows that each off-diagonal covariance is zero; diagonal entries are
one. This proves the exact covariance assertion, not merely a limit.

To prove the uniform subGaussian assertion, let W be any globally
symmetric law with density at most C relative to uniform signs. For any
t and integer k>=1,

```math
\mathbb E_W\langle t,X\rangle^{2k}
\le C\mathbb E\langle t,\varepsilon\rangle^{2k}
\le C(2k-1)!!\|t\|_2^{2k}.
```

The last inequality follows by expanding into even coordinate powers
and comparing Rademacher and Gaussian moments term by term. Odd moments
vanish. Consequently

```math
\mathbb E_W e^{\langle t,X\rangle}
\le1+C(e^{\|t\|_2^2/2}-1)\le e^{C\|t\|_2^2/2}.
```

Finally the change in E|Z| is delta times

```math
D=(a_2-a_1)b_0-(a_2-a_0)b_1+(a_1-a_0)b_2.
```

The cell endpoints alone give

```math
D\le-\tfrac45a_2+2a_1-\tfrac65a_0
\le-\tfrac45(19/10)^2+2(11/10)^2=-117/250<0.
```

This proves (5) with no numerical integration premise.

**Leading-scale Boolean consequence.** Form an `n by n` full sign
matrix B from independent rows with this law. For the pinned right word
`y=1`, the scalar maximum is exactly

```math
\max_{x\in\{\pm1\}^n}|x^TB\mathbf1|
=\sum_{i=1}^n\left|\sum_jB_{ij}\right|.
```

Its expected value differs from its covariance-matched independent
Gaussian counterpart by at least `c n^(3/2)` for some c>0 and all large
n. Indeed Gaussian row sums give `sqrt(2/pi)n^(3/2)`, while the ordinary
sign row sums converge to the same normalized value and (5) supplies
a fixed additional deficit.

Thus the linear-block boundary in Section 2 cannot in general be removed
using only covariance, full sign entries, and a uniform subGaussian
constant. The query set is a genuine pinned bipartite Boolean maximum,
with `2^n` left words. This is **not** a counterexample for the unrestricted
maximum over both words, nor a statement about actual minimizing children.
The distinction is essential.

A subsequent independently audited construction embeds this row-law
gap into the UNRESTRICTED scalar absolute cap of actual full sign
matrices with Q=O(N^(3/2)), using a weakly planted deterministic child:
[bounded-cap embedding](paper_localization_bounded_cap_embedding_2026_09_17.md).
It is not a near-minimizer example. The sharper all-offset comparison
for scalar sign-frame drivers is separate:
[common-Gibbs frame universality](paper_symmetric_frame_universality_2026_09_17.md).
