# Very thin tensor restrictions: universality versus randomization

Date: 2026-09-06. Seed-transfer track. These are actual finite signing
statements, not bounds on an auxiliary Gaussian certificate. They concern
an extremely thin range of principal restrictions and leave the important
intermediate retention range open.

Let `B` be a fixed symmetric full sign matrix of order `d>=2`, and let
`A=hollow(B)`. Write `Q(A)=max_x |x^T A x|/2` and
`beta(B)=max_(x,y signs)|x^T B y|`.

## 1. Every non-rank-one seed is universal after sufficiently thin selection

If `B` has rank one, symmetry implies `B=+-ss^T`. Every principal
restriction of every tensor power remains of this form. Its hollow cap
at order `n` is exactly `n(n-1)/2`.

If `B` does not have rank one, then EVERY hollow signing `C` of order
`n` is, up to a switching and global negation, a principal submatrix of
`hollow(B^(tensor n))`. Consequently

```math
\min_{|T|=n}Q(\operatorname{hollow}((B^{\otimes n})_T))=M_n.       (1)
```

Here the ambient order is `d^n`, so the retention is `n/d^n`. This does
not construct an optimizing target from the seed: it shows that searching
all such exceptionally thin selectors already contains the original
optimization problem exactly.

Proof and an explicit encoding follow. Fix a symbol `a` and switch/negate
the seed to the anchored matrix

```math
H_{uv}=B_{uv}B_{ua}B_{av}B_{aa},\qquad H_{au}=1.
```

Because `B` is not rank one, some `b,c` satisfy `H_bc=-1`; both differ
from `a`, though `b=c` is allowed. Encode vertex `i=1,...,n` by a word
`w_i` with its `i`th symbol `b`, all later symbols `a`, and earlier
symbols chosen in `{a,c}`. Construct the words in increasing `i`.
When setting the `j`th symbol of `w_i`, with `j<i`, the interaction
between `w_j,w_i` before coordinate `j` is already determined.
Coordinate `j` multiplies it by either `H_ba=1` or `H_bc=-1`, so it can
be made equal to the desired `C_ji`. Later coordinates contribute one,
because `w_j` has symbol `a` there. Thus all off-diagonal signs match.

The words are distinct: for `i<j`, coordinate `j` is `a` in `w_i` and
`b` in `w_j`. Finally

```math
(B^{\otimes n})_{w_i,w_j}
=B_{aa}^{\,n}D_iD_j(H^{\otimes n})_{w_i,w_j},
\qquad D_i=\prod_{\ell=1}^n B_{(w_i)_\ell,a}.
```

Switching and global negation preserve the absolute hollow cap. This
proves the construction and (1). Padding every word with the symbol `a`
also proves universality in `B^(tensor r)` whenever `r>=n`.

## 2. Quantitative mixing for a uniform random very thin restriction

Set `theta=beta(B)/d^2`. Non-rank-one is equivalent to `theta<1`.
Choose `n` independent uniform words in `[d]^r`, initially allowing
collisions. The vector of their `e=n(n-1)/2` off-diagonal tensor signs
has total-variation distance from `e` independent unbiased signs at most

```math
\Delta_{n,r}(B)=\frac12\sqrt{2^{n(n-1)/2}-1}\,\theta^r.          (2)
```

To prove this, index the Fourier characters of the sign-edge cube by
simple graphs `F` on `[n]`. The corresponding Fourier coefficient is
`t_F(B)^r`, where

```math
t_F(B)=d^{-n}\sum_{z_1,...,z_n\in[d]}
                    \prod_{ij\in E(F)}B_{z_i,z_j}.
```

For any nonempty `F`, choose one edge `ij` and condition on all other
colors. The remaining integrand is a fixed sign times
`B_uv f(u)g(v)`, for sign-valued functions `f,g`. Therefore
`|t_F(B)|<=beta(B)/d^2=theta`. Parseval and Cauchy--Schwarz on the
finite edge cube give (2). No asymptotic graph universality theorem is
being imported.

The probability of any repeated word is at most
`n(n-1)/(2d^r)`. Conditioning on distinct words makes their unordered
set a UNIFORM principal subset of size `n`. Thus the total-variation
distance for that actual principal restriction is at most

```math
\Delta_{n,r}(B)+\frac{n(n-1)}{2d^r}.                         (3)
```

For a fixed non-rank-one seed, this vanishes whenever, for example,
`r/n^2 -> infinity`. The sharper sufficient condition is

```math
r(-\log\theta)-\frac{n(n-1)}4\log2\longrightarrow+\infty,
\qquad n^2/d^r\longrightarrow0.                            (4)
```

Actual seed cap information can be used in this rate:

```math
\theta\le\min\left(1,\frac{4Q(A)+d}{d^2}\right).            (5)
```

Indeed `beta(B)<=beta(A)+d<=4Q(A)+d`. For a good large seed with
`Q(A)<=c d^(3/2)`, the last bound is `4c/sqrt(d)+1/d`.
This is a Boolean-cap-dependent loss-of-information estimate, not a
seed-preserving upper construction.

## 3. A cap obstruction above one half, with an elementary finite bound

For an iid unbiased hollow signing `R` of order `n`, take `x_1=1` and
sequentially choose `x_i` to make
`x_i sum_(j<i)R_ij x_j` nonnegative. The resulting energy is a sum of
independent increments with laws `|S_1|,...,|S_(n-1)|`, where `S_m` is
the sum of `m` independent unbiased signs. Independence follows because
at step `i` the newly exposed row edges are independent of the complete
past; multiplying them by the previously selected spins preserves their
uniform product law.

Put

```math
\mu_n=\sum_{m=1}^{n-1}\mathbb E|S_m|,
\qquad \mathbb E|S_m|
=\frac{m}{2^{m-1}}\binom{m-1}{\lfloor(m-1)/2\rfloor}.
```

The variance of the energy is at most `sum_m m=n(n-1)/2`. Hence, for
every `h>0`,

```math
\Pr\{Q(R)<\mu_n-h\}\le\frac{n(n-1)}{2h^2}.                (6)
```

Combining (3) and (6), a uniform `n`-subset `T` of `[d]^r` satisfies

```math
\Pr\{Q(\operatorname{hollow}((B^{\otimes r})_T))<\mu_n-h\}
\le\frac{n(n-1)}{2h^2}
   +\Delta_{n,r}(B)+\frac{n(n-1)}{2d^r}.                    (7)
```

The exact binomial formula and Stirling's formula give

```math
\frac{\mu_n}{n^{3/2}}\longrightarrow
\frac23\sqrt{\frac2\pi}=0.5319230405...>\frac12.
```

Thus if `n` tends to infinity and (4) holds, the random restrictions
have normalized cap at least `0.5319230405-o_P(1)`. In particular,
very thin UNIFORM RANDOM sampling cannot preserve a strict-subhalf seed
coefficient. This conclusion is about actual spin witnesses.

## 4. Scope relative to the earlier tensor obstruction

The even self-tensor classification proves divergence for a fixed
non-Hadamard seed, including after fixed-positive random retention.
The present calculation addresses a separate extreme: retained order
only `O(sqrt(r))` or smaller, compared with ambient order `d^r`.
It gives an actual lower obstruction there as well.

However, carefully selected restrictions can reproduce every signing
already at retained order `n=r`; they are not ruled out. Nor do these
arguments treat the large intermediate range, such as `n=exp(alpha r)`
with `0<alpha<log d`. No new upper bound or convergence conclusion is
asserted. Exact selector encodings, Fourier coefficients, and small
greedy laws are checked by
`computations/transfer_seed_thin_tensor_verify_2026_09_06.py`.
