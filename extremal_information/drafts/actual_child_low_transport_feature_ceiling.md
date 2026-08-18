# Low-transport row features cannot expose extensive actual-child dependence

Status: **rigorous actual-law method ceiling**.  This note applies to the
negative-disorder bridge escort induced by the contracted-temperature
optimizing children.  It proves that the finite matching-bit certificates
in `actual_child_row_product_global_certificate.md` cannot amplify to a
linear directed row-product gap merely by repeating the same kind of
low-sensitivity feature at larger orders.

The result is not a surrogate counterexample.  It is a uniform theorem for
the exact actual-child escort.  It bounds what a rowwise data-processing
certificate can see; it does not upper-bound the full directed projection.

## 1. Setup

Let the bridge cube be

```math
\Omega=\prod_{i=1}^m\{-1,1\}^{n_i},
\qquad U=\bigotimes_iU_i,
```

and let `L` be the actual parent pressure in one fixed orientation at raw
temperature `t`.  Flipping one bridge coordinate changes `L` by at most
`2t`.  For `lambda>0`, write

```math
{dq\over dU}(B)={e^{-\lambda L(B)}\over E_Ue^{-\lambda L}}. \tag{LT.1}
```

Choose finite row maps

```math
\phi_i:\{-1,1\}^{n_i}\longrightarrow Y_i,
\qquad
Q=(\phi_1,\ldots,\phi_m)_\#q,
\qquad
\nu=(\phi_1,\ldots,\phi_m)_\#U=\bigotimes_i\nu_i.       \tag{LT.2}
```

Their exposed reverse-product information is

```math
J(Q)=\inf_{P=\otimes_iP_i}D(P\Vert Q).                  \tag{LT.3}
```

The image density is governed by the exact conditional negative moment

```math
{dQ\over d\nu}(y)
={W_\phi(y)\over E_\nu W_\phi},
\qquad
W_\phi(y)=E_U[e^{-\lambda L(B)}\mid\phi(B)=y].        \tag{LT.3a}
```

At the physical channel amplitude of RR.1--RR.5, `e^(-lambda L)` can be
replaced by the inverse actual-child likelihood `p_t^(-lambda)`; the peeled
sector constant cancels.  Thus (LT.3a), not a surrogate row law, is the
exact quantity whose coarse product projection is being bounded.

For two probability measures on a row cube, let `W_infinity` denote the
least essential supremum of Hamming distance over all couplings.  Define the
uniform-fibre transport radius

```math
\delta_i=\max_{y,y'\in\operatorname{supp}\nu_i}
 W_\infty\!\left(
 U_i(\,\cdot\mid\phi_i=y),
 U_i(\,\cdot\mid\phi_i=y')
 \right).                                               \tag{LT.4}
```

This quantity measures the geometric sensitivity of a declared row
feature, not the cardinality of its output alphabet.

## 2. Reverse-information ceiling

**Theorem LT.1 (fibre-transport ceiling).**  In the setup above,

```math
\boxed{
J(Q)\le D(\nu\Vert Q)
\le {\lambda^2t^2\over2}\sum_{i=1}^m\delta_i^2.}       \tag{LT.5}
```

In particular, at the contracted scale `t=beta/sqrt(N)`, a comparable
split with

```math
\sum_i\delta_i^2=o(N^2)                                \tag{LT.6}
```

has `J(Q)=o(N)`.  Hence a coarse data-processing proof of
`I_lambda^leftarrow=Omega(N)` requires
`sum_i delta_i^2=Omega(N^2)`.  For row-homogeneous feature architectures
this is root-mean-square fibre radius `Omega(sqrt(N))`; without homogeneity
the transport budget may instead be concentrated on fewer rows.

*Proof.*  Put

```math
g(y)=\log{dQ\over d\nu}(y)
=\log E_U[e^{-\lambda L(B)}\mid\phi(B)=y]
 -\log E_Ue^{-\lambda L}.                              \tag{LT.7}
```

Fix `y_(-i)` and two values `y_i,y_i'`.  Couple the two conditional uniform
row fibres with Hamming displacement at most `delta_i`, and use the identity
coupling on every other row.  Along every coupled pair, the coordinate-flip
bound and a Hamming path give

```math
e^{-2\lambda t\delta_i}e^{-\lambda L(B)}
\le e^{-\lambda L(B')}
\le e^{2\lambda t\delta_i}e^{-\lambda L(B)}.          \tag{LT.8}
```

Taking coupled expectations proves that the `i`th coordinate oscillation of
`g` is at most `2 lambda t delta_i`.

The bounded-difference exponential inequality on the product law `nu`
therefore gives

```math
\log E_\nu e^{g-E_\nu g}
\le{1\over8}\sum_i(2\lambda t\delta_i)^2.             \tag{LT.9}
```

Because `E_nu e^g=1`, the left side equals `-E_nu g=D(nu||Q)`.
Finally, `nu` is itself row-product, so it is an admissible competitor in
(LT.3).  This proves (LT.5). `square`

## 3. Coordinate-budget and linear-feature corollaries

There is a sharper bound when the row features inspect only declared bridge
coordinates.

**Corollary LT.2 (coordinate budget).**  Suppose `phi_i` depends only on a
set `S_i` of bridge coordinates in row `i`, and put
`K=sum_i|S_i|`.  Then

```math
\boxed{J(Q)\le{\lambda^2t^2K\over2}.}                 \tag{LT.10}
```

Thus, at contracted temperature, every feature family reading
`K=o(N^2)` bridge entries in total exposes only `o(N)` directed dependence.

*Proof.*  Let `q_S` be the marginal of `q` on the union of the selected
coordinates.  Marginalizing the positive density in (LT.1) preserves the
per-bit log-oscillation bound `2 lambda t`.  Apply (LT.9) with the selected
bits as the independent coordinates to obtain

```math
D(U_S\Vert q_S)\le\lambda^2t^2K/2.                    \tag{LT.11}
```

The image of `U_S` under the row maps is the row-product law `nu`.
KL data processing gives `D(nu||Q)<=D(U_S||q_S)`, and (LT.10) follows.
`square`

**Corollary LT.3 (bounded-rank Walsh features).**  If each `phi_i` records
at most `k_i` independent Walsh parities of its row, then one may choose
`delta_i<=k_i`.  Consequently

```math
J(Q)\le{\lambda^2t^2\over2}\sum_i k_i^2.              \tag{LT.12}
```

In particular, a fixed number of row parities on each of `Theta(N)` rows
has `J(Q)=O(1)`, even if those parities collectively inspect every bridge
entry.

*Proof.*  Row parities are a linear map over `F_2`.  Choose pivot
coordinates for a basis of its image.  Any difference of two attainable
syndromes can be corrected by flipping at most `k_i` pivots; translation by
that correction is a measure-preserving bijection of the corresponding
uniform affine fibres.  Hence (LT.4) has `delta_i<=k_i`, and LT.1 applies.
`square`

## 4. Consequence for the finite actual-child certificate

The order-eight certificate in GC.17 uses one Walsh parity in each bridge
row (in fact one selected bridge bit per row).  LT.2 gives the all-order
ceiling

```math
J(Q_N)\le {\lambda^2\beta^2\over2N}m=O(1)             \tag{LT.13}
```

for any repetition of that architecture at comparable splits.  Its
certified finite gap `1.075` is real, but it cannot become the required
`Omega(N)` gap by a block-amplification argument that retains only boundedly
many Walsh parities per row.  Such an argument would contradict (LT.12).

For reference, the rigorous coarse gaps from the actual-child finite
certificate and the corresponding LT.10 ceilings are

| actual optimized-child law | selected bits | certified `J(Q)` | LT.10 ceiling |
|---|---:|---:|---:|
| `N=8`, `beta=4`, `lambda=5.382104195764755` | 4 | `[1.075,1.075620]` | `115.8681822962743` |
| `N=9`, `beta=2`, `lambda=1` | 4 | `[.006637668616855,.006637668616858]` | `8/9` |
| `N=9`, `beta=4`, `lambda=1` | 4 | `[.027345865444539,.027345865444542]` | `32/9` |

The order-eight bound is loose because its raw amplitude is `sqrt(2)`, not
weak.  LT.10 is nevertheless decisive at fixed `beta,lambda` as `N` grows:
any `Theta(N)` matching-coordinate family has only `O(1)` exposed reverse
information.

The surviving scalable alternative is narrower than CA.2 suggested:

> **High-transport actual-child feature lemma.**  Either find row features
> with `Omega(sqrt(N))` uniform-fibre transport radius on a positive density
> of rows whose actual escort image has a linear reverse-product gap, or
> prove directly that the full actual escort admits an `o(N)` reverse
> product certificate without passing through such features.

High transport is necessary for this certificate architecture, not
sufficient.  The feature must still have a compressed law and a global
reverse-projection proof; merely choosing a sensitive function would retain
no less information than the original row response.

## 5. The transport threshold is a real structural boundary

The preceding theorem is about the actual escort.  The following auxiliary
example is not asserted to be child-induced; it shows why its
`Omega(sqrt(N))` transport alternative cannot be removed using only the
actual law's already proved local regularity properties.

Let `r` be odd, take `r` rows of `r` signs, pair all but one row (leaving the
last row fair and independent), and set

```math
F_r(R)=\tanh\left(r^{-1/2}\sum_{j=1}^rR_j\right).      \tag{LT.14}
```

Pair the rows and, on every pair, use the density

```math
{dq_{\gamma,r}^{(2)}\over d(U_r\otimes U_r)}(R,S)
={\exp\{\gamma F_r(R)F_r(S)\}\over Z_{\gamma,r}},
\qquad 0<\gamma\ll1.                                  \tag{LT.15}
```

Let `q_(gamma,r)` be the product of these pair laws.

**Proposition LT.4 (high-transport sharpness).**  The family (LT.15) has:

1. microscopic log-density flip oscillation at most `2 gamma/sqrt(r)`;
2. uniformly bounded conditional and marginalized row `D_2`;
3. unbiased individual signs and global sign symmetry; but
4. an extensive reverse row-product projection,

   ```math
   \inf_{P=\otimes_iP_i}D(P\Vert q_{\gamma,r})
   \ge c_\gamma r.                                    \tag{LT.16}
   ```

The dependence is exposed by the one-bit row feature
`sign(sum_j R_j)`.  Its two uniform fibres have transport radius at least
`(r+1)/2`, so it lies decisively outside LT.1's low-transport regime.

*Proof.*  The first claim follows from the `1`-Lipschitz property of `tanh`.
Every pair potential lies in `[-gamma,gamma]`; after conditioning outside a
row and marginalizing arbitrary coordinates within it, the remaining
density relative to its fair law is bounded between `e^(-2gamma)` and
`e^(2gamma)`.  Its Renyi-two divergence is therefore at most `4 gamma`.
Oddness of `F_r` proves the symmetry claims.

Push one row pair through the signs of its magnetizations.  Its law is

```math
Q_r(a,b)={1+\rho_rab\over4},
\qquad
\rho_r={E\sinh(\gamma U_rV_r)\over
             E\cosh(\gamma U_rV_r)},                 \tag{LT.17}
```

where `U_r,V_r` are independent copies of `|F_r(R)|` under the fair row
law.  The central limit theorem and bounded convergence give
`rho_r -> rho_gamma>0`.  Choose `gamma` small enough that
`rho_r<tanh(1)` for all sufficiently large `r`.  The binary rectangle
certificate CA.3 then makes the uniform product the global reverse-product
minimizer and gives

```math
J(Q_r)=-{1\over2}\log(1-\rho_r^2)\ge c_\gamma>0.      \tag{LT.18}
```

KL data processing and exact additivity across the independent row pairs
prove (LT.16).  Finally, the all-plus word in the positive-majority fibre is
at Hamming distance at least `(r+1)/2` from every negative-majority word.
It has positive conditional mass, so every essential-supremum transport
coupling has at least that radius. `square`

LT.4 proves a sharp qualitative point.  Microscopic flip bounds, tight row
Renyi complexity, and symmetry do not select between tight and extensive
directed dependence.  The missing actual-child information must control
the presence or absence of genuinely collective, high-transport row
features.

More precisely, the smallest remaining overlap/cumulant statement is a
**uniform conditional-negative-moment lemma**.  It must provide high-
transport, fixed-alphabet row maps and control the potential

```math
g_N(y)=\log E_U[p_t(B)^{-\lambda}\mid\phi(B)=y]       \tag{LT.19}
```

well enough to prove either bounded-block domination with `o(N)` total
remainder or a globally contractive coarse product problem with value
`cN-o(N)`.  The control must be uniform over whole feature fibres at
`t=beta/sqrt(N)`.  A finite Taylor expansion in child overlaps is
insufficient: `Theta(N^2)` weak bridge coordinates can leave an order-`N`
collective remainder, while specifying all such conditional moments simply
reconstructs the bridge cube.  LT.1--LT.3 therefore stop the matching-bit
route before that reconstruction occurs and isolate the genuinely
aggregate lemma that would be needed to escape.
