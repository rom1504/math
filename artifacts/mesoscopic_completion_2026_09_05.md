# Mesoscopic completion and separately paid bridges

Status: **Verified**, with independent proof audits. This improves a
near-order completion coefficient, not the rigorous limiting interval.

## Exact reversal identity

For hollow symmetric sign matrices A,D and an n by r sign bridge B, put
H_A(x)=sum_{i<j} a_ij x_i x_j and Q(A)=max_x |H_A(x)|. Then

```math
Q\!\begin{pmatrix}A&B\\B^{\mathsf T}&D\end{pmatrix}
=\max_{x,y}\left(|H_A(x)+H_D(y)|+|x^{\mathsf T}By|\right). \tag{MC.1}
```

Indeed y -> -y preserves both child energies and reverses the bridge;
max(|u+v|,|u-v|)=|u|+|v|. Any gain over separately bounding the two
nonnegative terms must control their joint large-value locations.

## Theorem 1: uniform rectangular estimate

Define

```math
R_{n,r}=\min_{B\in\{\pm1\}^{n\times r}}\|B\|_{\infty\to1},
\quad
\mu_r=\mathbb E\left|\sum_{j=1}^r\varepsilon_j\right|
=\frac r{2^{r-1}}\binom{r-1}{\lfloor(r-1)/2\rfloor}.
```

For all n,r >= 1,

```math
\max(n\mu_r,r\mu_n)\le R_{n,r}\le
\min\left\{
n\mu_r+\sqrt{2nr(r-1)\log2},
r\mu_n+\sqrt{2nr(n-1)\log2}
\right\}.                                                   \tag{MC.2}
```

**Proof.** For every fixed B, averaging over uniform y gives
E_y sum_i |(By)_i|=n mu_r, proving the first lower bound; transpose for
the other. For the upper bound take all nr entries independently uniform.
For fixed y, F_y=sum_i |(By)_i| has mean n mu_r. Changing one entry changes
F_y by at most 2. The bounded-difference exponential-moment inequality,
from the Doob martingale and Hoeffding's lemma, gives

```math
\mathbb E\exp(t(F_y-\mathbb EF_y))\le \exp(nrt^2/2).
```

There are 2^(r-1) distinct tests, since F_y=F_{-y}. Jensen and log-sum-exp
therefore give, for t>0,

```math
\mathbb E\max_y F_y\le
n\mu_r+\frac{(r-1)\log2}{t}+\frac{nrt}{2}.
```

Optimize t. When r=1, F_y=n identically, giving zero remainder. Some
deterministic bridge has norm at most the expectation. Transposition gives
the other alternative, possibly using a different bridge. QED.

## Corollary 2: exact optimizing-child completion

Choosing exact minimizers at orders n,r and applying (MC.1) yields

```math
M_{n+r}\le M_n+M_r+
\min\left\{
n\mu_r+\sqrt{2nr(r-1)\log2},
r\mu_n+\sqrt{2nr(n-1)\log2}
\right\}.                                                   \tag{MC.3}
```

Since mu_r >= sqrt(r/2), uniformly for r=o(n),

```math
R_{n,r}=n\mu_r\left(1+O(\sqrt{r/n})\right).                    \tag{MC.4}
```

If also r -> infinity, Stirling's formula gives
mu_r=sqrt(2r/pi)(1+O(1/r)). The established bound M_r=O(r^(3/2)) then gives

```math
M_{n+r}\le M_n+\left(\sqrt{2/\pi}+o(1)\right)n\sqrt r.         \tag{MC.5}
```

For every principal restriction, Q(A[S]) <= Q(A): extend each prescribed
spin assignment by independent uniform spins and take conditional
expectations. Thus M_n <= M_{n+r}. For c_n=M_n/n^(3/2),

```math
\left(\frac n{n+r}\right)^{3/2}c_n
\le c_{n+r}
\le \left(\frac n{n+r}\right)^{3/2}c_n+
\left(\sqrt{2/\pi}+o(1)\right)\sqrt{r/n}.                     \tag{MC.6}
```

For bounded r use (MC.3) instead of its central-limit coefficient.

## Comparison with the archive and limitation

The archived iid completion has leading cost sqrt(2 log 2) n sqrt(r).
Equation (MC.5) improves that coefficient to sqrt(2/pi). The archived
deterministic rectangular construction has error O(r 2^r), which can be
better for tiny or specially divisible r. Qualitative continuity across
o(n) changes of order was already known; it is not claimed as new.

At equal splits every bridge costs at least
(sqrt(2/pi)+o(1)) n^(3/2). A putative recurrence
b_{2n} <= 2b_n+o(n), b_n=M_n^(2/3), allows cap at most
2^(3/2) M_n+o(n^(3/2)). After separately paying 2M_n, the bridge allowance
is at most (sqrt(2)-1+o(1)) n^(3/2), using limsup c_n <= 1/2.

Equivalently the best separately paid certificate has b-scale defect at
least

```math
\left[(1+\sqrt{2/\pi})^{2/3}-2^{1/3}+o(1)\right]n
=(0.218646\ldots+o(1))n.                                    \tag{MC.7}
```

This method-class obstruction is already in
[cross_order_method_linear_floors.md](cross_order_method_linear_floors.md).
It is not a lower bound on the true composition defect.

A sufficient unresolved statement is: for some delta>0 and K<infinity,
all sufficiently large comparable m,n have exact minimizing children and
a sign bridge satisfying

```math
Q\!\begin{pmatrix}A&B\\B^{\mathsf T}&D\end{pmatrix}^{2/3}
\le M_m^{2/3}+M_n^{2/3}+K(m+n)^{1-\delta}.                   \tag{MC.8}
```

Take comparable to mean 1/2 <= m/n <= 2. Balanced merge trees of q equal
seed blocks of order k accumulate O(q k^(1-delta)) error: at merge scale
2^j k there are O(q/2^j) nodes, so the sum is a convergent geometric
series in 2^(-j delta). Dividing by qk gives O(k^(-delta)).
The leftover fewer than k vertices cost o(1) after normalization by
(MC.3). Choosing seed orders along a liminf subsequence then proves
convergence of b_n/n and hence c_n. Generic o(N) error without a summable
rate does not justify that argument.

No estimate of the form (MC.8) was obtained. Finite exhaustive checks in
the accompanying verifier test normalizations, not asymptotic truth.
