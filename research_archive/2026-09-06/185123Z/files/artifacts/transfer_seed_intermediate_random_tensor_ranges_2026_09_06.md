# Two rigorous intermediate random-tensor restriction ranges

Date: 2026-09-06. Seed-transfer track. All cap statements below concern
actual hollow sign matrices and actual Boolean witnesses. No spin-glass
or random-matrix universality theorem is assumed.

Fix a symmetric full sign seed `B` of order `d>=2`. Put

```math
N=d^r,\qquad K_r=B^{\otimes r},\qquad
\theta=\beta(B)/d^2,\qquad
\beta(B)=\max_{x,y\in\{-1,1\}^d}|x^TBy|.
```

Take a uniform principal subset of size `n` from `K_r`, and hollow it.
The following two actual lower mechanisms cover distinct retention
ranges. They do not cover every intermediate range.

## 1. Exponentially growing thin restrictions have cap above one half

If `B` is not rank one, then `theta<1`. If `n` tends to infinity and

```math
n^{14}\theta^r\longrightarrow0,                            (1)
```

then the normalized cap of a uniform `n`-principal restriction satisfies

```math
\frac{Q(\operatorname{hollow}((K_r)_T))}{n^{3/2}}
\ge\frac{377}{750}-o_P(1)
=0.5026666666\ldots-o_P(1).                                 (2)
```

In fact the proof has limiting coefficient greater than
`0.502762704756`, but (2) is the simple exact rational statement.
For the scale `n=floor(d^(alpha r))`, a sufficient range is

```math
0<\alpha<\alpha_{\rm bulk}(B)
:=\frac{\log_d(1/\theta)}{14}.                             (3)
```

The term 'bulk' is only a label for this explicit greedy-moment proof;
it is not a claim of a spectral or SK universality limit.

### 1.1 A cap-dependent absolute correlation estimate

For independently uniform columns `z_1,...,z_q` of the base seed, define

```math
\eta_q=\mathbb E\left|\frac1d\sum_a
                         \prod_{j=1}^qB_{a,z_j}\right|.
```

For every `q>=1`, `eta_q<=theta`. Indeed, condition on `z_1,...,z_(q-1)`
and put `f_a=prod_(j<q)B_(a,z_j)`, which is a sign vector. Averaging
the remaining absolute sum over `z_q` gives
`||B^T f||_1/d^2<=beta(B)/d^2`.

For `q` independently sampled tensor words `W_j` indexed by distinct
indices (their values may coincide), this implies

```math
\mathbb E\left|\mathbb E_W\prod_{j=1}^q(K_r)_{W,W_j}\right|
=\eta_q^r\le\theta^r.                                     (4)
```

The word distribution is a product over the `r` tensor coordinates,
so both the inner correlation and its expected absolute value factor
exactly. No asymptotic independence assertion is involved.

If `A=hollow(B)`, actual seed cap gives the useful upper estimate

```math
\theta\le\min\left(1,\frac{4Q(A)+d}{d^2}\right).             (5)
```

This follows from polarization `beta(A)<=4Q(A)` and the diagonal cost
`beta(B)<=beta(A)+d`. Thus (1) can be checked using only an actual cap
upper bound when the right side of (5) is less than one.

### 1.2 Uniform moment control even for adaptive greedy spins

First sample words independently, with replacement. After seeing the
first `i` words, let `x_1,...,x_i` be ANY signs depending on their entire
history. For an independent next word `W`, set

```math
F_i(W)=i^{-1/2}\sum_{j=1}^i x_j(K_r)_{W,W_j},\qquad
m_k(i)=\mathbb E(S_i/\sqrt i)^k,
```

where `S_i` is a sum of `i` independent unbiased signs. For even `k`,
expand the `k`th power into ordered tuples of indices. If all indices
have even multiplicity, the summand equals one, exactly as for `S_i`.
Otherwise, let `J` be the nonempty set of indices having odd
multiplicity. The summand is a coefficient of absolute value one times
`E_W prod_(j in J)(K_r)_(W,W_j)`.

Let `R_(i,k)` be `i^(-k/2)` times the sum of the absolute values of
these non-even-multiplicity correlations. Pointwise in the sampled
history and UNIFORMLY over every adaptive sign choice,

```math
|\mathbb E_WF_i(W)^k-m_k(i)|\le R_{i,k},\qquad
\mathbb E R_{i,k}\le i^{k/2}\theta^r.                       (6)
```

The expectation bound uses (4) for each index tuple; there are at most
`i^k` tuples. Adaptation is harmless because coefficient signs were
removed before expectation. This is the key distinction from an
unjustified conditional Gaussian approximation to the greedy process.

### 1.3 An exact global polynomial minorant of the absolute value

Let the seven positive rational numbers `s_j` be

```text
7991/10000, 16067/10000, 24324/10000, 32891/10000,
41962/10000, 51901/10000, 63639/10000.
```

There is a unique polynomial `q(t)=sum_(k=1)^14 c_k t^k` satisfying

```math
q(0)=0,\qquad q(s_j^2)=s_j,\qquad q'(s_j^2)=1/(2s_j).
```

Its coefficients are rational. It obeys `q(t)<=sqrt(t)` for every
`t>=0`, not merely on a sampled grid. To see this, Hermite-interpolate
`sqrt(t+epsilon)` at the same nodes, with a single node at zero and
double nodes at the seven positive squares. The order-15 interpolation
remainder is its positive 15th derivative at an intermediate point
times
`t prod_j(t-s_j^2)^2/15!`, which is nonnegative. Let `epsilon` decrease
to zero. The coefficients converge to those of `q`, proving the
global inequality. Consequently `P(x)=q(x^2)<=|x|` globally.

The accompanying exact rational verifier constructs `q`, checks all
interpolation equalities, and proves

```math
\sum_{k=1}^{14}|c_k|<9,\qquad
g_P:=\mathbb E P(G)=\sum_{k=1}^{14}c_k(2k-1)!!
>\frac{150828811427}{200000000000}>\frac{377}{500}.           (7)
```

Here `G` is standard normal only as a compact description of its exact
even moments; the displayed computation is rational. The downward
endpoint in (7) is `0.754144057135`.

### 1.4 Drift and martingale estimates for actual sequential greedy energy

Choose `x_1=1` and sequentially choose `x_(i+1)` to make its newly
added edge sum nonnegative. Write its increment as
`Y_i=sqrt(i)|F_i(W_(i+1))|` and its conditional mean as `a_i`.
Let

```math
D_n=\sum_{i=1}^{n-1}\sqrt i\sum_{k=1}^{14}c_k m_{2k}(i).
```

The pointwise minorant and (6) imply

```math
\sum_i a_i\ge D_n-H_n,\quad H_n\ge0,\qquad
\mathbb E H_n\le9 n^{3/2}n^{14}\theta^r.                    (8)
```

By the even-moment expansion of an iid sign sum,
`m_(2k)(i)->(2k-1)!!`. Therefore

```math
D_n/n^{3/2}\longrightarrow(2/3)g_P>377/750.                 (9)
```

The centered increments form a martingale difference sequence. From
the case `k=2` of (6),

```math
\mathbb E\left(\sum_i(Y_i-a_i)\right)^2
\le\sum_{i=1}^{n-1}(i+i^2\theta^r)
\le n^2/2+n^3\theta^r/3.                                  (10)
```

In particular, for every `h>0`, Markov and Chebyshev give the finite
bound

```math
\Pr\{\text{greedy energy}<D_n-2h n^{3/2}\}
\le\frac{9n^{14}\theta^r}{h}
  +\frac{1/(2n)+\theta^r/3}{h^2}.                          (11)
```

Finally condition on distinct sampled words. The added probability
error is at most `n(n-1)/(2d^r)`, and this conditional sample is a
uniform principal restriction. This error also vanishes under (1):
`beta(B)>=d`, so `theta>=1/d` and `n^14/d^r<=n^14 theta^r`.
Equations (7)--(11) prove (2).

## 2. A coherent inherited witness covers a complementary exponent range

Suppose now that `B` is non-Hadamard and non-rank-one. The explicit
vector witness from `transfer_seed_tensor_power_instability_2026_09_06.md`
is

```math
V(B)=d^{-1/2}\sum_j\|(B^TB)_{:,j}\|,
\qquad d^{3/2}<V(B)<d^2.
```

Write `v=log_d V(B)` and `ell=log_d ||B||op`, so `3/2<v<2` and
`ell<1`. For every even `r`, the full hollow tensor has an actual spin
witness with signed energy at least
`V(B)^r/(2pi)-d^r/4`. Choose a sign `sigma_r` so this inherited energy
is positive; it need not be the positive extremum before that choice.

Set

```math
\alpha_{\rm coh}(B)=4-2v\in(0,1).
```

For even `r` and `n=floor(d^(alpha r))`, with
`alpha_coh(B)<alpha<1`, the normalized cap of a uniform `n`-principal
restriction tends to INFINITY in probability.

Here is the retention-sensitive variance argument. For a fixed spin
witness of a hollow matrix `L` on `N` vertices, put `p=n/N`,
`pi_j=(n)_j/(N)_j`, and `P_L(x)=x^T Lx/2`. Its restricted signed energy `Z` has mean
`pi_2` times its original signed energy and variance

```math
\operatorname{Var}Z
=(\pi_2-2\pi_3+\pi_4)\binom N2
 +(\pi_3-\pi_4)\|Lx\|^2
 +(\pi_4-\pi_2^2)P_L(x)^2
\le p^2N^2/2+p^3N\|L\|_{op}^2.                            (12)
```

The last coefficient is nonpositive. The first two are nonnegative
and bounded by `pi_2,pi_3`, respectively; `pi_j<=p^j`.
Use `L=hollow(B^(tensor r))`, so
`||L||op<=||B||op^r+1`. Relative to the square of the inherited mean,
the two errors in (12) are bounded by constant multiples of

```math
d^{r(4-2\alpha-2v)},\qquad
d^{r(2-\alpha+2\ell-2v)}.                                  (13)
```

At `alpha=alpha_coh`, their exponents are `2v-4<0` and `2ell-2<0`.
They remain negative above that value. The inherited normalized mean
is at least

```math
(1/(2\pi)-o(1))\,d^{r(\alpha/2+v-2)},                     (14)
```

which diverges precisely above `alpha_coh`. Chebyshev proves the claim.

At the critical exponent there is also a quantitative statement: if
`n=floor(C d^(alpha_coh r))`, for fixed `C>0`, then the same proof gives
normalized cap at least `sqrt(C)/(2pi)-o_P(1)`. In particular `C>pi^2`
already forces a lower bound strictly above one half. The stronger
direct same-spin tensor witness can improve this prefactor, but is not
needed for either range theorem here.

## 3. What is and is not covered

For every non-rank-one seed, low positive exponents (3) have actual cap
above one half with high probability. For every non-Hadamard such seed,
high exponents above `alpha_coh` have diverging cap along even powers.
Rank-one seeds have cap `n(n-1)/2` under every restriction, separately.

There is generally a wide remaining window between the two proven
ranges. For a good large seed with `Q(hollow(B))<=c d^(3/2)`, (5)
gives a low-exponent range approaching `alpha<1/28` as `d` grows.
For the actual strict-upper weave seeds, the uniform bound
`V(B)/d^(3/2)>1.0154` gives
`alpha_coh<1-2 log(1.0154)/log d` for sufficiently large seed orders.
Neither statement closes the intervening range, rules out exceptional
deterministic selectors, or resolves convergence of the original minima.

The source
`computations/transfer_seed_intermediate_tensor_moments_2026_09_06.py`
contains the complete rational-node specification and exact replay.

### Bounded literature scope check

A primary-source abstract check found no imported theorem covering the
remaining window. [Lu--Yau, arXiv:2205.06308](https://arxiv.org/abs/2205.06308)
establishes spectral-distribution equivalence for spherical/Gaussian
inner-product kernels in polynomial sample/dimension scaling.
[Do--Vu, arXiv:1206.3763](https://arxiv.org/abs/1206.3763) addresses
spectral universality for inner-product and distance kernels. Those
stated conclusions do not give a Boolean-cap comparison for the present
fixed-alphabet product kernel. This was a bounded scope screen, not a
claim that the entire relevant literature has been excluded.
