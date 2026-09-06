# Uniform Gaussian row profiles and a lower bound on the actual annealed weave

Date: 2026-09-06. The finite smooth-test theorem and its consequences below
are proved directly. They concern arbitrary retained Hadamard rows, not only
the randomized recursive ensemble. The annealed conclusion is NOT a
high-probability lower bound on a constructed signing's cap.

## 1. A finite, basis-uniform smooth-test theorem

Let `H` be any real Hadamard matrix of order `m`, retain any `k` of its
rows to form `H_T`, and let `x` be uniform in `{+-1}^k`. Put

```math
v_j=k^{-1/2}\sum_{a\in T}H_{aj}x_a,\qquad
\nu_x=\frac1{2m}\sum_{j=1}^m(\delta_{v_j}+\delta_{-v_j}),\qquad
R=k^{-1}H_T^TH_T.
```

For an even function `f in C^4(R)` with
`M=max_(0<=r<=4)||f^(r)||_infinity<infinity`, let
`Y_f=integral f d nu_x`. If `G` is standard Gaussian, then

```math
\left|\mathbb E Y_f-\mathbb Ef(G)\right|\le\frac{M}{6k},
\qquad \operatorname{Var}(Y_f)\le\frac{4M^2}{k}.                 (1)
```

Both bounds are uniform in `m`, the Hadamard, and the retained rows.

Replace the `k` Rademacher inputs one at a time by independent standard
Gaussians. Their first three moments agree; Taylor expansion through degree
three with fourth-derivative remainder gives a replacement error at most
`(1+3)M/(24k^2)` per coordinate for a single test. This proves the mean
bound. For a product test `f(v_i)f(v_j)`, every directional fourth derivative
is bounded by `M^2(|u_1|+|u_2|)^4`. Each input direction has coordinates
of magnitude `1/sqrt(k)`, so the total product expectation error is at most
`8M^2/(3k)`. Comparing the two products of means costs at most `M^2/(3k)`.
Thus its Rademacher covariance is at most its Gaussian covariance plus
`3M^2/k`.

The Gaussian pair has correlation `R_ij`. The orthonormal Hermite expansion
of an even `f` contains only even degrees, so

```math
0\le\operatorname{Cov}(f(G_i),f(G_j))
=\sum_{r\ge1} f_{2r}^2 R_{ij}^{2r}
\le R_{ij}^2\operatorname{Var}(f(G))\le M^2 R_{ij}^2.
```

This identity also holds at correlations `+-1` by continuity. Exact row
orthogonality gives

```math
\sum_{ij}R_{ij}^2=\operatorname{tr}R^2=m^2/k.
```

Average the covariance bound over the `m^2` pairs. This proves (1), including
the diagonal pairs. No small bound on each individual correlation is assumed.

## 2. Typical empirical spectra converge in W2

As `k->infinity`, uniformly over the choices above,

```math
W_2(\nu_x,N(0,1))\longrightarrow0\quad\text{in probability}.     (2)
```

For any finite set of smooth even tests, (1) and Chebyshev give convergence
to their Gaussian values in probability. Smooth tests determine weak
convergence of symmetric measures. Tightness is automatic because the
second moment is EXACTLY one for every spin:

```math
\int z^2\nu_x(dz)=\|H_T^Tx\|^2/(km)=1.
```

To see the uniform W2 assertion without assuming tail uniformity, fix a large
continuity radius `R` and use a smooth bounded approximation from below to
`z^2 1_(|z|<=R)`. Its empirical integral is close to the Gaussian integral
with uniformly high probability by (1). Subtract from the exact moment one
to control the empirical tail energy. Finite smooth approximations of interval
indicators then control the compact part. First choose the truncation and
tests, then let `k` grow. This is the usual weak-plus-second-moment criterion
with every uniformity supplied by (1).

In particular for any fixed finite magnitude quantizer `q_R,delta` which
clips at `R` and rounds on a mesh of size at most `delta`, typical spins have
the corresponding Gaussian bin frequencies up to `o(1)`. Their mean squared
quantization error is, with uniformly high probability, at most the Gaussian
tail payment plus `O(delta^2)` and an arbitrarily small further error. Boundaries
are chosen at Gaussian continuity points. The proportion of such spins is
`1-o(1)`, so their count is `(1-o(1))2^k`, not just one specially chosen spin.

## 3. Actual joint annealed pressure, before Finner factorization

Use arbitrary deterministic order-m Hadamards `H_i`, arbitrary selectors
`T_i` of common size `k`, independent uniform column permutations, and
independent fair off-diagonal outer signs `S_ij=S_ji`. Diagonal outer signs
may be fixed. Construct the full symmetric weave `W` and let

```math
D_+(x)=\sum_{i,j}(h_i(j)-S_{ij}h_j(i))^2,
\qquad h_i=H_i[T_i,:]^Tx_i,\qquad N=mk.
```

Column permutations are incorporated in each `h_i`. Row orthogonality gives
`D_+(x)=2(m^2k-x^TWx)`. For fixed `t>0`, define the ACTUAL annealed partition

```math
\mathcal Z_{m,k}(t)=
\mathbb E_{S,\mathrm{perm}}\sum_{x\in\{\pm1\}^{mk}}
                      e^{-tD_+(x)/(2k)}.
```

If `k/m->p in (0,1]`, then, uniformly in the deterministic bases/selectors,

```math
\boxed{\quad
\liminf_{m\to\infty}\frac1{m^2}\log\mathcal Z_{m,k}(t)
\ge p\log2+g_t(1).
\quad}                                                       (3)
```

Here is the complete reduction to the finite endpoint-count theorem in
`transfer_reconstruction_homogeneous_edge_count_2026_09_06.md`.
First fix a large truncation and a fine finite quantizer. By (2), each fibre
has `(1-o(1))2^k` spins whose row magnitudes have nearly the same quantized
Gaussian histogram and mean squared quantization error at most `epsilon`.
All errors are uniform in the fibre. Choose a tolerance tending to zero
slowly enough that the histogram statement holds for this fraction of spins.

Condition each column permutation's diagonal coordinate to have magnitude
at most a fixed large `R_0`. Since each spectrum has squared norm `m`, this
event has probability at least `1-R_0^(-2)` in each row. Removing that bounded
coordinate changes the row histogram by at most one entry. The product
conditioning probability is `exp(-O(m))`, and its diagonal defect factor is
at least `exp(-2t R_0^2 m)`. Both payments are `o(m^2)` in the logarithm.
Conditional on the diagonal choices, the remaining row assignments are
independent uniform permutations of their remaining multisets.

Averaging an off-diagonal edge sign gives the positive folded kernel

```math
K_t(a,b)=\tfrac12(e^{-t(a-b)^2}+e^{-t(a+b)^2}).
```

Quantization changes the TOTAL logarithmic kernel weight by at most
`O(t sqrt(epsilon))m^2`, uniformly in every endpoint permutation. Indeed
`|partial_a log K_t(a,b)|<=2t(|a|+|b|)` and similarly for `b`.
Integrate along the segment between the original and quantized endpoints,
sum over edges, and apply Cauchy--Schwarz to directed endpoints. The total
squared original energy is at most `m^2`, total squared endpoint error is
at most `epsilon m^2`, and quantized energy is at most
`(1+sqrt(epsilon))^2 m^2`. No cancellation of signed edge errors is needed.

The finite endpoint-count theorem therefore gives a lower exponent
`Phi_t(nu_quantized)-O(t sqrt(epsilon))` for the actual edge expectation.
Its proof allows the small row-to-row histogram discrepancies present here.
The magnitude-kernel transport is exactly the full symmetric-source
transport after relative-sign folding. As truncation grows and mesh shrinks,
`nu_quantized -> N(0,1)` in W2, hence `Phi_t(nu_quantized)->g_t(1)`.

The product number of chosen spin blocks has logarithm
`mk log2 + m log(1-o(1))=p m^2 log2+o(m^2)`. Combining these bounds proves
(3). The limit order is fixed truncation/mesh, then `m->infinity`, then
truncation/mesh refinement. One can equivalently use `mk log2/m^2->p log2`.

## 4. What (3) excludes

For target full cap `c N^(3/2)`, the usual exponential union bound uses
`exp[t(1-2c sqrt(p))m^2] mathcal Z_(m,k)(t)`, with an additional factor two
if both energy signs are counted. Equation (3) gives its lower exponential
rate

```math
p\log2+g_t(1)+t(1-2c\sqrt p).
```

Its optimized threshold is the same `sqrt(15)/8` derived in the exact finite
permanent-floor theorem. Thus at fixed positive retention and temperature,
even EXACT evaluation of the joint annealed partition cannot beat a strictly
smaller constant. This is stronger in mechanism, but narrower in parameter
uniformity, than the finite pointwise row-certificate floor.

It does not prove that typical signings have cap at least this value.
The expected exponential sum can be dominated by rare signings, or by very
large clusters of correlated spins. A cap theorem would require a further
probabilistic argument. The original convergence question and optimality of
the weave ensemble remain open.
