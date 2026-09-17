# Wave 22 memo: Johnson entropy and neighbor overlap

## Status

The channel inequalities and Gibbs neighbor identities below are **proved**.
The Johnson modified log-Sobolev and log-Sobolev constants are imported from
primary sources with their normalizations translated explicitly.  The `A_9`
table is **Numerical**, from exhaustive enumeration by
`tmp/johnson_overlap_r22_verify.py`.  Nothing here proves the required
rate--distortion bound; the fixed-temperature child-Gibbs specialization is a
scoped obstruction.

## 1. Johnson normalization and channel entropy

Let `Omega_{n,m}` be the uniform `m`-slice.  From `S`, choose `u in S` and
`v notin S` independently and uniformly, and put `S'=S-u+v`.  Thus `S'` is a
uniform Johnson neighbor and the degree is

```math
d_J=m(n-m).
```

For functions on the slice, Salez's normalization is

```math
\mathcal E_J(f,g)
=\frac1{2n}\sum_{1\le i<j\le n}
\mathbb E[(\nabla_{ij}f)(\nabla_{ij}g)]
=\frac{d_J}{2n}\mathbb E_{S,S'}
[(f(S')-f(S))(g(S')-g(S))].
```

In this normalization, the Poincare constant is exactly `tau_rel=1`, or,
equivalently, the discrete neighbor chain has spectral gap `n/d_J`.  The
optimal modified log-Sobolev constant obeys

```math
\frac12\le\tau_{mls}(n,m)\le1,
```

and

```math
\operatorname{Ent}_{U_m}(f)
\le\tau_{mls}\mathcal E_J(f,\log f).
```

These statements, including the normalization, are Lemma 1.3 and equations
(1.2)--(1.4) of Justin Salez,
[A sharp log-Sobolev inequality for the multislice](https://doi.org/10.5802/ahl.99),
Annales Henri Lebesgue 4 (2021), 1143--1161.  Its proof derives the upper
bound by coarsening the modified LSI of Caputo--Dai Pra--Posta.  Thus the
factor below is not an unnormalized appeal to Johnson expansion.

Let `K_S` be arbitrary probability laws on one common finite output space,
let

```math
\overline K=\mathbb E_SK_S,
\qquad f_y(S)=\frac{K_S(y)}{\overline K(y)}
```

on outputs with positive mixture mass.  Then `E_S f_y=1` and

```math
I(S;D)=\sum_y\overline K(y)\operatorname{Ent}_{U_m}(f_y).
```

Apply the modified LSI to every `f_y` and sum.  Since

```math
\sum_y\overline K(y)(f_y(S')-f_y(S))
(\log f_y(S')-\log f_y(S))
=J(K_S,K_{S'}),
```

where

```math
J(P,Q)=D_{KL}(P\Vert Q)+D_{KL}(Q\Vert P),
```

one obtains the exact channel inequality

```math
\boxed{
I(S;D)
\le\tau_{mls}(n,m)\frac{m(n-m)}{2n}
\mathbb E_{S,S'}J(K_S,K_{S'})
\le\frac{m(n-m)}{2n}\mathbb E_{S,S'}J(K_S,K_{S'}).}
```

The inequality remains true with value infinity if neighboring laws are
singular.  It applies in particular to any exact- or near-ground channel and
therefore to the common-prior objective after its channel is constructed.

There is a Hellinger version which remains finite under singularity.  Define

```math
H^2(P,Q)=1-\sum_y\sqrt{P(y)Q(y)}
=\frac12\sum_y(\sqrt{P(y)}-\sqrt{Q(y)})^2.
```

The ordinary log-Sobolev inequality gives

```math
\boxed{
I(S;D)
\le\tau_{ls}(n,m)\frac{m(n-m)}n
\mathbb E_{S,S'}H^2(K_S,K_{S'}).}
```

Lee and Yau's two-color Bernoulli--Laplace theorem, in the same Salez
normalization, supplies the explicit upper bound

```math
\boxed{
\tau_{ls}(n,m)
\le\frac2{\log2}\log\frac{n^2}{m(n-m)}.}
```

See Theorem 5 of Lee--Yau,
[Logarithmic Sobolev inequality for some models of random walks](https://doi.org/10.1214/aop/1022855885),
Annals of Probability 26 (1998), 1855--1873.  At ratios bounded away from
zero and one, `tau_ls=O(1)` and both channel inequalities have a prefactor of
order `n`.

Consequently, the Wave 21 information target
`I(S;D)=O(n^{1/2-2c})` follows from either

```math
\mathbb E_{adj}J(K_S,K_{S'})=O(n^{-1/2-2c})
```

or

```math
\mathbb E_{adj}H^2(K_S,K_{S'})=O(n^{-1/2-2c})
```

at a fixed ratio.  This is the exact local overlap scale to seek.

## 2. Exact neighboring child-Gibbs formulas

Use the common output space of full oriented projective cuts.  For an
`m`-set `S`, define the child Gibbs law, including uniform outside extension,
by

```math
K_S(\omega)
=\frac{e^{\beta c_A(S,\omega)}}{Z_S},
\qquad
Z_S=\sum_{\omega\in\mathcal D_n}e^{\beta c_A(S,\omega)}.
```

The energy ignores outside relative spins, so this definition is exactly the
local child Gibbs law followed by uniform outside extension.

For adjacent selectors

```math
S=C\cup\{u\},\qquad T=C\cup\{v\},
```

put

```math
h_u^C(x)=\sum_{i\in C}a_{ui}x_i,
\qquad
g(\omega)=c_A(S,\omega)-c_A(T,\omega)
=2\sigma[x_uh_u^C(x)-x_vh_v^C(x)].
```

Directly from the exponential densities,

```math
\boxed{
\begin{aligned}
D_{KL}(K_S\Vert K_T)
&=\beta\mathbb E_{K_S}g-\log\frac{Z_S}{Z_T},\\
J(K_S,K_T)
&=\beta(\mathbb E_{K_S}g-\mathbb E_{K_T}g).
\end{aligned}}
```

If

```math
K_t(\omega)\propto
\exp\{\beta[c_A(T,\omega)+tg(\omega)]\},
```

then differentiation of the log partition function gives

```math
\boxed{
J(K_S,K_T)=\beta^2\int_0^1
\operatorname{Var}_{K_t}(g)\,dt.}
```

The Bhattacharyya affinity is the midpoint partition ratio

```math
\boxed{
\operatorname{BC}(K_S,K_T)
=\frac{Z_{(c_S+c_T)/2}}{\sqrt{Z_SZ_T}}.}
```

It also has an exact common-core form.  Summing first over `x_u,x_v` and
cancelling all outside multiplicities gives

```math
\boxed{
\operatorname{BC}(K_S,K_T)
=\frac{
\sum_{\sigma,x_C}e^{\beta c_A(C,\sigma,x_C)}
\cosh(\beta h_u^C)\cosh(\beta h_v^C)}
{\sqrt{
\left[\sum_{\sigma,x_C}e^{\beta c_A(C,\sigma,x_C)}
\cosh(2\beta h_u^C)\right]
\left[\sum_{\sigma,x_C}e^{\beta c_A(C,\sigma,x_C)}
\cosh(2\beta h_v^C)\right]}}.}
```

This is the required overlap quantity.  Scalar child pressures give only the
two denominator factors, not the correlated numerator.

## 3. Average neighbor divergence is a scalar wrong-way identity

Let

```math
\overline E_\beta
=\mathbb E_{S\sim U_m}\mathbb E_{K_S}c_A(S,D).
```

There is an exact identity over directed Johnson edges:

```math
\boxed{
\mathbb E_{S,S'}D_{KL}(K_S\Vert K_{S'})
=\frac{2\beta}{m}\overline E_\beta,
\qquad
\mathbb E_{S,S'}J(K_S,K_{S'})
=\frac{4\beta}{m}\overline E_\beta.}
```

Proof: the directed edge law is reversible, so the average log-partition
difference vanishes and the two directed KL averages agree.  Given `S`, write
`S'=S-u+v`.  Under `K_S`, the outside sign `x_v` is independent uniform, so
the newly introduced star has mean zero.  Averaging the removed vertex `u`
removes exactly `2/m` of the ordered internal energy, because every internal
edge is incident to two possible removed vertices.  Hence

```math
\mathbb E_{S,S',D\sim K_S}[c_A(S,D)-c_A(S',D)]
=\frac2m\overline E_\beta.
```

Substitution proves the identity.

The modified-LSI certificate consequently collapses to

```math
\boxed{
I(S;D)
\le2\tau_{mls}(n,m)\frac{\beta(n-m)}n\overline E_\beta
\le2\frac{\beta(n-m)}n\overline E_\beta.}
```

This is the wrong scale at useful fixed temperature.  At fixed ratios and
competitive energy `overline E_beta=O(n^{3/2})`, it is `O(beta n^{3/2})`.
At high temperature, orientation symmetry and Parseval give the fixed-`m`
Taylor expansion

```math
\overline E_\beta
=2m(m-1)\beta+O(\beta^3)
```

as `beta -> 0`.  Thus neighbor Jeffreys divergence has infinitesimal leading
term `8(m-1)beta^2`.  This Taylor statement alone is not uniform in `m`, so
it must not be used at an `n`-dependent temperature without a remainder
bound.

There is nevertheless a direct, nonperturbative scale obstruction if the
Gibbs law is near-ground.  If `overline E_beta=Theta(n^{3/2})`, the exact
identity above and the required local `O(n^{-1/2-2c})` scale force

```math
\beta=O(n^{-1-2c}).
```

That forced temperature rigorously contradicts the near-ground premise.
Indeed, for every symmetric zero-diagonal `m` by `m` matrix `B`, the
Rademacher quadratic form `Q=x^T Bx` obeys

```math
\log\mathbb E_x e^{tQ}\le16t^2\lVert B\rVert_F^2
\quad\text{when}\quad
|t|\le\frac1{\sqrt{32}\lVert B\rVert_{op}}.
```

For completeness, randomly bipartition the coordinates and use Jensen to
decouple `Q` into four times one cross-block bilinear form.  Conditional
Rademacher subgaussianity bounds its mgf by
`E_x exp(8t^2 ||B_{I,I^c}x||^2)`.  Gaussian linearization and
`-log(1-z)<=2z` for `z<=1/2` then give the displayed bound.  Convexity of
the log-mgf consequently gives, with a harmless factor-two shrinkage of the
range,

```math
\mathbb E_{K_S}c_A(S,D)
\le64\beta\lVert A_S\rVert_F^2
\le64\beta m(m-1).
```

Thus `beta=O(n^{-1-2c})` implies
`overline E_beta=O(n^{1-2c})`, not `Theta(n^{3/2})`.  This excludes the
canonical child-Gibbs law as a Johnson-overlap certificate at the required
scale.  The fixed-`m` expansion is consistent with the rigorous obstruction.
Conversely, fixed-temperature Gibbs laws have subleading entropy distortion
but their neighbor overlap is not small.  This scopes out the direct
combination of child Gibbs, uniform outside extension, and Johnson MLSI.  It
does not rule out a different minimizer-specific near-ground channel whose
adjacent laws are coupled much more smoothly.

The pointwise inequality

```math
(a-b)(\log a-\log b)\ge4(\sqrt a-\sqrt b)^2
```

also gives `H^2(P,Q)<=J(P,Q)/8`; hence replacing modified LSI by the standard
Hellinger LSI does not repair the generic Gibbs scaling without additional
control of the common-core numerator.

## 4. `A_9`, `m=8` finite audit

Exhausting all 512 full oriented projective cuts and all 36 undirected
Johnson edges gives:

| `beta` | mean distortion | `I(S;D)` | mean `H^2` | mean Jeffreys |
|---:|---:|---:|---:|---:|
| `0.10` | `13.845448` | `0.109850` | `0.061672` | `0.507728` |
| `0.50` | `1.723287` | `0.867729` | `0.492830` | `5.569178` |
| `1.00` | `0.220752` | `1.253406` | `0.716891` | `11.889624` |
| `2.00` | `0.004026` | `1.372158` | `0.814450` | `23.995974` |

The verifier checks the exact average-KL identity numerically at every row.
As the Gibbs distortion becomes small, adjacent Hellinger distance becomes
order one and Jeffreys divergence grows.  This is a finite diagnostic, not an
asymptotic counterexample.

## Updated route conclusion

Johnson functional inequalities do convert local overlap into exactly the
sub-square-root information required by (10.694).  The necessary neighbor
scale is explicit.  However, the canonical child Gibbs channel has an exact
average divergence identity with the wrong sign and scale.  A positive
continuation must construct a different near-ground channel/common prior
whose outputs vary by only `O(n^{-1/2-2c})` divergence across a typical vertex
exchange, or prove a correlated common-core affinity estimate not visible to
scalar pressure.
