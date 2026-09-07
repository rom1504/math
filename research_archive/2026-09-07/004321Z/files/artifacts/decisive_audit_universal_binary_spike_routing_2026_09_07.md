# Universal binary-channel spike routing for arbitrary Hadamard bases

Date: 2026-09-07. Status: finite-channel proof below, independent audit
requested. This extends a restricted part of the Walsh routing theorem
to ARBITRARY Hadamard bases. It is a lower bound for a specified actual
annealed signing partition, not a lower bound for the original minimum.

## 1. Exact scope and conclusion

Fix `k>=2`, a symmetric full sign seed `A` of order `k` with `A_ii=1`,
`p in (0,1]`, and `t>0`. Let `m=kn` be Hadamard orders with `n` even
and tending to infinity, and set `ell=floor(pm)`. At each group `alpha`
use ANY real Hadamard `H_alpha` of order `m`, followed by an independent
uniform permutation of ALL its columns. Bases can differ by group or
be randomized; the estimates below are uniform in their choices.

Use the actual grouped weave and common selector in
`decisive_bridge_vector_seed_weave_unitary_obstruction_2026_09_06.md`.
Selectors are uniform `ell`-subsets, all retained `k`-spin tuples are
summed, and each inter-group edge has one independent shared fair sign.
Let `Z_full` be the resulting FULL positive defect partition, including
the within-group factors. Then, with binary entropy `h_2`,

```math
 \liminf_{n\to\infty}\frac{\log Z_{\rm full}}{n^2}
 \ge k^2\max_{0\le b\le1}
       \left\{p\,h_2\left(\frac{1+b}{2}\right)
                                      +g_t(1-pb^2)\right\}.           (1)
```

Here `g_t(v)` is the scalar Gaussian self-transport value. Equivalently,
define

```math
 E_t^{\rm mb}(\nu_p)=\max_{0\le b\le1}
 \left\{g_t(1-pb^2)
           -p\left[\log2-h_2\left(\frac{1+b}{2}\right)\right]\right\};
```

then the right side of (1) is `k^2[p log2+E_t^mb(nu_p)]`.
Section6 proves that this is EXACTLY the mean-variance envelope restricted
to channels revealing no information about the hole mask. No assertion
`E_t^mb=E_t` is made for a general ternary source.

## 2. A mean channel requiring only actual columns, not column products

Fix an interior `b`; endpoints follow by continuity. In each base choose
`2k` distinct columns. Write the first `k` as `e_i(g)=H(g,c_i)` and the
next `k` as `z_j(g)=H(g,d_j)`. Use their joint `2k`-bit pattern to partition
the physical rows into at most `2^(2k)` cells. Their joint empirical law
need not be uniform, independent, or even convergent.

At row `g`, draw a common mask `M` with probability `p_m=ell/m`. Given
`M=1`, draw the `k` signs independently with means `b e_i(g)`, and set
`Y_i=M sigma_i`. Conditional on the cell,

```math
 f_i(g):=\mathbb E Y_i(g)=p_m b e_i(g),
 H(Y(g)\mid\text{cell})=h_2(p_m)
                              +p_m k h_2((1+b)/2).                     (2)
```

The entropy in (2) is INDEPENDENT of the cell pattern. This is what
replaces the Walsh uniform-label requirement. The means are literal
single Hadamard columns, so no closure of products of columns is used.

In each cell impose an integer conditional type approximating this law,
with the total number of nonzero physical vectors exactly `ell`. A fixed
finite number of count corrections suffices after ordinary rounding.
The entire exact-type set has logarithmic cardinality

```math
 m[h_2(p_m)+p_m k h_2((1+b)/2)]-O_k(\log m).                          (3)
```

The constant may depend on fixed `p,b,k`, but not on the base or on the
cell sizes. Empty cells are ignored; small cells contribute only the
usual fixed-alphabet type-count logarithmic errors. Every resulting
configuration is an actual COMMON selector and actual retained spins.
Dividing by `binom(m,ell)` changes (3) to
`ell k h_2((1+b)/2)-O_k(log m)`.

## 3. Matched diagonal spikes and a controlled self block

Normalize spectra as `h_i=H^T Y_i/sqrt(ell)`. Exact conditional types
fix all sums tested against every chosen label column. The rounding
error in any such unnormalized sum is bounded by a constant depending
only on the finite cell alphabet. Hadamard orthogonality therefore gives

```math
 h_i(c_j)=\frac{p_m b m}{\sqrt\ell}\,1_{i=j}+O_k(m^{-1/2}),
 \qquad h_i(d_j)=O_k(m^{-1/2}).                                      (4)
```

Choose one perfect matching of the `n` group vertices. At each group
condition its full column permutation to place `c_i` in column `i` of
the partner block, and `d_j` in column `j` of the self block. This fixes
only `2k` columns, with probability at least `m^(-2k)` per group.
The total logarithmic routing cost is `O_k(n log m)=o(n^2)`.

The partner block is the SAME scalar diagonal matrix at both ends,
up to `O_k(m^(-1/2))`. Since `R_A` acts trivially on diagonal matrices,
each matched folded kernel is at least
`(1/2)exp[-O_(k,t)(1/m)]`. These `n/2` edges cost only `O(n)` in log
partition, although they carry an order-`n^2` total energetic mass.

The self block has squared norm `O_k(1/m)`. For any symmetric full
within-group signing its signed-transpose map is orthogonal. Its full
defect factor is at least `exp[-2t ||self block||^2]`. Thus the sum of
ALL within-group log losses is `O_(k,t)(n/m)=O_(k,t)(1)`. No omitted
within-group factor or uncontrolled tail deletion is used in (1).

## 4. Uniform Gaussian residual profile without Walsh algebra

Before exact-type repair, rows are independent bounded vectors. Let
`V_ij(g)=Cov(Y_i(g),Y_j(g)|cell)`. Then

```math
 V_{ii}=p_m-p_m^2b^2,\qquad
 V_{ij}(g)=p_m(1-p_m)b^2 e_i(g)e_j(g)\quad(i\ne j).                   (5)
```

The average offdiagonal covariance vanishes by orthogonality of the
active columns. The covariance between residual Fourier coordinates
is the corresponding entry of

```math
 C^{ij}=\ell^{-1}H^T\operatorname{diag}(V_{ij})H,
 \qquad \|C^{ij}\|_F^2
             =\frac{m^2}{\ell^2}\sum_gV_{ij}(g)^2=O_{p,k}(m).       (6)
```

Every same-frequency cross-row covariance is zero. For distinct random
frequencies, (6) makes all cross-covariances tend to zero in probability,
uniformly over the base. The diagonal variances tend to
`v=1-pb^2`. Random formation of residual ordered `k`-column blocks
therefore makes a generic block, or generic pair of blocks, have covariance
approaching `v I_(k^2)` or its two-block version.

Each such vector is a sum of `m` independent bounded vectors scaled by
`ell^(-1/2)`. The fixed-dimensional characteristic-function expansion has
uniform error `O(m^(-1/2))`. Together with (6) this proves the one-block
Gaussian limit and the two-block independent Gaussian limit. The empirical
profile converges weakly in probability, uniformly in the Hadamard bases.

Repair the preliminary physical sample to the prescribed exact cell
types by editing whole physical vectors. With probability tending to one,
only `O_k(sqrt(m log m))` edits are needed. Parseval shows that the mean
squared spectral perturbation is `o(1)`. The repair has `exp(o(m))`
preimages. Intersecting the original information-typical set with the
good-profile set and applying this repair therefore leaves weighted
exact-type count at least the exponential rate in (3), minus `o(m)`.
This count statement is uniform in the bases and includes averaging
over random residual block formation.

After repair, Parseval fixes total spectral energy to `km=k^2n`.
Equation (4) puts `k^2n pb^2+O_k(1)` of it in the matched spike block
and only `O_k(1/m)` in the self block. Thus the residual empirical second
moment tends to `k^2(1-pb^2)`. Weak convergence upgrades to W2 convergence.

Decompose the remaining column permutation into random formation of
ordered blocks and then a uniform permutation of their destinations.
The good-profile event concerns only the formed multiset, so the second
permutation remains uniform after conditioning. This is essential for
the next pressure step.

## 5. Residual pressure and normalization

For the positive folded kernel
`K_R(x,y)=exp[-t(||x||^2+||y||^2)]cosh(2t x^T R y)`,
the exact global estimate is

```math
 |\log w(X)-\log w(Z)|
       \le2t\|X-Z\|(\|X\|+\|Z\|).                                  (7)
```

Couple empirical row multisets by optimal matchings. Uniform W2 proximity
to a common Gaussian and the deterministic Parseval energy bound make
the normalized residual pressure uniformly close to the homogeneous
Gaussian pressure. To justify this without a boundedness assumption,
first approximate the Gaussian in W2 by a fixed bounded finite alphabet,
apply (7), and use the exact positive finite-profile theorem. The removed
one matching and self slots cause only `O(n)` error at this bounded stage:
insert dummy colors, delete the prescribed slots, and use permutation-
equivariant minimal row recoloring. This preserves uniform target row
arrangements and changes only `O(n)` incidences. Then remove the W2
approximation by (7).

The homogeneous folded Gaussian value is `k^2 g_t(v)`, independent of
the orthogonal seed reflection. One direct calculation unfolds the fair
sign by the sign/endpoint-reflection entropy averaging, and maximizes
the Gaussian correlation `rho` in

```math
 g_t(v)=\max_{0\le\rho<1}
       \{-tv(1-\rho)+\tfrac14\log(1-\rho^2)\}.                      (8)
```

For `v=0` use the limiting value zero. Multiplying the good source counts
from the `n` groups gives normalized source exponent
`(n ell k/n^2)h_2((1+b)/2)->pk^2h_2((1+b)/2)`. Add the residual value
and the subleading routing, matched-edge, and self-edge costs. This proves
(1) for fixed `b`; optimize and take endpoint limits. All approximation
parameters are fixed before the order limit.

## 6. Binary channels optimize exactly among mask-blind channels

For the scalar source `X=M sigma/sqrt(p)`, call a channel mask-blind if
`M` and its output `L` are independent. Put
`a(L)=E[sigma|M=1,L]`. Then

```math
 I(X;L)=p\,\mathbb E[\log2-h_2((1+a(L))/2)],\qquad
 \mathbb E\operatorname{Var}(X|L)=1-p\,\mathbb E a(L)^2.             (9)
```

The function `u->h_2((1+sqrt(u))/2)` is concave on `[0,1]`, as follows
from its convergent series
`log2-sum_(j>=1) u^j/[2j(2j-1)]` and continuity at the endpoints.
Thus at fixed `E a^2=b^2`, the information in (9) is minimized by
constant `|a|=b`. Equal reflected labels with posterior means `+b,-b`
attain that value and preserve the fair sign marginal. This is exactly
the binary symmetric channel used in Sections2--5. Hence (1) realizes
the complete mask-blind envelope, not merely a selected binary example.

A strict gain of the full `E_t(nu_p)` over this envelope MUST use information
about the mask. That statement does not prove such a gain exists at any
particular parameter, or that all mask-informative channels require Walsh
bases. The full finite-latent Walsh theorem remains strictly stronger in
its channel class, while (1) is strictly broader in its admissible bases.
