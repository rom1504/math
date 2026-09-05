# Finite ancestor anchors: a universal Gaussian fixed-point extension

Date: 2026-09-05. The fixed-point argument, finite formula, and explicit
numerical instance passed an independent variational-agent audit. This uses
the already proved universal tree-energy theorem, not AMP universality.

## 1. A finite anchored subspace

Let `U h_T=G_T` be the Gaussian isometry from the hierarchical tree theorem.
Choose a finite ancestor-closed family `A` of trees, including the edge.
Its coordinates `G=(G_T:T in A)` are independent standard Gaussians. Each
`h_T` is a normalized even Hermite monomial in ancestor coordinates, hence
is a function of `G`. The functions `h_T`, `T in A`, are orthonormal.

Choose real coefficients `rho_T` with `R²=sum rho_T²<1`, and write
`s=sqrt(1-R²)`. Let `h(G,Z)` be a jointly even Gaussian function satisfying

```math
E h^2=1,\qquad E[h h_T]=0\quad(T\in A),\qquad
D_Z:=E|\partial_Zh|^2<1,
```

where `Z` is a standard Gaussian independent of `G`. These are finite
conditions on a polynomial in the explicit construction below.

On the unit sphere of first Gaussian chaos orthogonal to all `G_T`, set

```math
\mathcal R(Z)=\mathcal U h(G,Z).
```

This preserves the sphere: isometry gives its norm and the displayed
orthogonality gives its orthogonality to each anchor. For two innovations
of correlation `q`, expand in the innovation coordinate only to obtain

```math
E h(G,Z)h(G,Z')=K(q)=\sum_{l\ge0}w_lq^l,\quad
w_l\ge0,\quad\sum_lw_l=1,\quad\sum_l l w_l=D_Z.
```

The anchor is independent of the joint pair `(Z,Z')`. Thus the map has
Lipschitz constant at most `sqrt(D_Z)<1` on this complete sphere. It has
a unique fixed point `Z_*`. Consequently the same-space unit Gaussian

```math
V=\sum_{T\in A}\rho_TG_T+sZ_*
  =\mathcal U g,\qquad
g=\sum_{T\in A}\rho_T h_T+s h(G,Z_*)
```

exists. In particular `Eg²=1`. The ancestor condition is important: it
makes every removed feature a function of the fixed anchor, independent
of whichever innovation is being iterated. Arbitrary Gaussian anchors
with inverse images depending on the innovation would not justify this
argument.

## 2. Exact finite resolvent, independent of ambient anchor count

Put `chi=1{|sum rho_T G_T+sZ|<=alpha}` and `p=2Phi(alpha)-1`.
In the normalized multivariate Hermite basis its coefficient at total
even degree `d` and multi-index `k` is

```math
c_k=\beta_d\sqrt{d!/\prod_j k_j!}\prod_j v_j^{k_j},\quad
\beta_0=p,\quad
\beta_d=-2\phi(\alpha)He_{d-1}(\alpha)/\sqrt{d!},
```

where `v=(rho,s)`. Fix an even cutoff `D` at least as large as every
anchor feature's degree. Delete exactly the Hermite monomials `h_T`,
`T in A`, and retain the other coefficients of degree at most `D`.
Divide each retained coefficient by `a+l`, where `l` is its innovation
degree and `a>0`, and normalize the resulting finite polynomial to norm
one. This defines `h` and ensures all required anchor orthogonalities.

Write `c_T=E[chi h_T]`. If the top of `T` has child multiplicities
`m_j` and `d_T=sum m_j`, then

```math
c_T=\beta_{d_T}\sqrt{d_T!/\prod_jm_j!}
                 \prod_j\rho_j^{m_j}.
```

For `T=edge`, the empty product gives `c_T=p`. Define grouped squared
coefficients

```math
B_l=\sum_{\substack{d\le D\text{ even}\\d\ge l}}
       \beta_d^2\binom dl(R^2)^{d-l}s^{2l}
       -1_{\{l=0\}}\sum_{T\in A}c_T^2.
```

They are nonnegative: they are sums of the remaining Hermite coefficient
squares, not a formal signed relaxation. Therefore, exactly,

```math
D_Z=\frac{\sum_l l B_l/(a+l)^2}{\sum_l B_l/(a+l)^2},\qquad
w:=E[V\mathcal U\chi]
 =\sum_T\rho_Tc_T+
 s\frac{\sum_l B_l/(a+l)}{\sqrt{\sum_l B_l/(a+l)^2}}.
```

The only exponential-looking coefficient sum has collapsed to `D+1`
scalar groups and the finite anchor list. This is an identity from the
multinomial theorem; it does not approximate the function by its bulk
moments. The actual function being defined is a finite polynomial.

## 3. Original signing consequence

Let `W=U chi`. The pair `(V,W)` is jointly Gaussian with variance `(1,p)`
and covariance `w`. Provided `0<w<sqrt(p)`, set `sigma=sqrt(p-w²)`.
The universal tree theorem gives

```math
\liminf_n\frac{M_n}{n^{3/2}}
\ge E|W|1_{\{|V|>\alpha\}}
=2w\phi(\alpha)[2\Phi(w\alpha/\sigma)-1]
 +4\sqrt p\,\phi(0)\overline\Phi(\alpha\sqrt p/\sigma).
```

The matrix-limit order is unchanged: first construct the Gaussian fixed
point, then approximate its bounded final mask by finitely many smooth
tree coordinates, then send the matrix order to infinity. No intermediate
unbounded polynomial is used as a spin mean. The contraction alone is
not being imported as a theorem about a growing-depth matrix algorithm.

## 4. Exact numerical instance

The exploratory programs `fresh_two_anchor_fixed_point_diagnostic.py` and
`fresh_finite_anchor_fixed_point_diagnostic.py` find finite candidates
above 0.4304 with nine anchors, and above 0.4307 with 21 anchors. These
floating-point searches are not certificates. The separate program
`computations/fresh_finite_anchor_fixed_point_certificate.py` fixes 21
explicit anchors, degree 200, `alpha=361/500`, `a=17/5`, and the rational
anchor coefficients displayed in its saved JSON. It uses exact Fraction
arithmetic and outward Gaussian-function intervals, with no quadrature
or omitted Hermite tail. It proves

```math
D_Z<0.990872910662429,\qquad
\liminf_n M_n/n^{3/2}\ge0.4306581794055286.
```

The full enclosed certificate interval is

```text
[0.430658179405528602724053804634711026327173238336325190455581,
 0.430658179405528602724053804634711026327173238336325190455840].
```

The independent audit reconstructed the anchor/innovation independence,
all removed Hermite coefficients, grouped norm and derivative sums, and
the original-limit normalization. An independent floating computation
agreed with the exact interval. The director executed the complete exact
certificate. Convergence remains open, and this construction remains
inside the previously proved nonnegative-mask method ceiling below 0.4495.
