# Anchored Gaussian tree fixed points and a finite certificate

Date: 2026-09-05. Status: the theorem, grouped coefficient formulas, and
exact-arithmetic certificate have passed an independent algebra-agent audit.
This is a Gaussian function-space construction, not an assumed AMP limit on
the original signing.

## 1. Common probability space and Gaussianization

Use the countable independent Gaussian family `Z_T` indexed by the rooted
odd-degree trees of the hierarchical tree-energy theorem. The normalized
child Hermite monomials `h_T` form the complete jointly-even Gaussian `L²`
basis. Let

\[
\mathcal U F=\sum_T E[Fh_T]Z_T
\]

be the isometry to first Gaussian chaos, and put `G_0=U1=Z_edge`.

Fix `0<rho<1` and `s=sqrt(1-rho²)`. Let `h(g,z)` be jointly even,
`h(-g,-z)=h(g,z)`, and suppose, for independent standard normals `G,Z`,

\[
Eh(G,Z)=0,\quad Eh(G,Z)^2=1,\quad
D_Z:=E|\partial_z h(G,Z)|^2<1.                        \tag{1}
\]

A Gaussian weak derivative suffices. No boundedness or pointwise smoothness
is required of this intermediate function, which is used only in Gaussian
space. The finite-polynomial example below satisfies all integrability
requirements automatically.

## 2. Exact anchored fixed point by contraction

Consider the unit sphere of first Gaussian chaos orthogonal to `G_0`. This
is a complete metric space in `L²`. For an innovation `Z` in that sphere,
define

\[
\mathcal R(Z)=\mathcal U[h(G_0,Z)].
\]

The input is jointly even under the global Gaussian sign flip. Its norm is
one because `G_0,Z` are independent standard normals. Its mean is zero, so
`R(Z)` is orthogonal to `U1=G_0`. Thus the map preserves the domain exactly.

For two such innovations `Z,Z'` with correlation `q`, the anchor `G_0` is
independent of their joint Gaussian pair. Expand `h` in the second variable:

\[
h(g,z)=\sum_{l\ge0} a_l(g)H_l(z)/\sqrt{l!}.
\]

Then

\[
K(q)=E h(G_0,Z)h(G_0,Z')=\sum_l\|a_l\|_2^2q^l,
\quad K(1)=1,\quad |K'(q)|\le D_Z<1
\]

on the whole interval `[-1,1]`. It follows that

\[
\|R(Z)-R(Z')\|_2^2
 =2[1-K(q)]\le 2D_Z(1-q)=D_Z\|Z-Z'\|_2^2.
\]

Banach's theorem gives a unique actual innovation fixed point `Z_*`, and
therefore a unique anchored unit Gaussian

\[
\boxed{\quad
V=\rho G_0+sZ_*,\qquad
V=\mathcal U\!\left[\rho+s h\!\left(G_0,
                 \frac{V-\rho G_0}{s}\right)\right].
\quad}                                               \tag{2}
\]

The equality holds in `L²` on the same countable Gaussian space. Every
iterate lies in first chaos, so all iterates, the anchor, and the limit are
jointly Gaussian; this is stronger than a distributional fixed-point claim.

## 3. Final Boolean certificate

For a fixed `alpha>0`, let `H=1_{|V|≤alpha}` and `W=UH`. Then `(V,W)` is
jointly Gaussian, with

\[
p:=\operatorname{Var}W=2\Phi(\alpha)-1,
\]

and, by the fixed-point identity and the isometry,

\[
\boxed{\quad
w:=E VW=\rho p+s E[h(G,Z)1_{|\rho G+sZ|\le\alpha}].
\quad}                                               \tag{3}
\]

The hierarchical tree-energy theorem and finite approximation of the
bounded mask give

\[
\liminf_n\frac{M_n}{n^{3/2}}
 \ge E|W|1_{|V|>\alpha}.
\]

Explicitly, for `w≥0` and `sigma=sqrt(p-w²)`, this is

\[
\boxed{\quad
L(p,w,\alpha)=
2w\phi(\alpha)[2\Phi(w\alpha/\sigma)-1]
 +2\sqrt{2p/\pi}\,
      \overline\Phi(\alpha\sqrt p/\sigma).
\quad}                                               \tag{4}
\]

All infinite operations precede finite approximation in Gaussian space.
For each desired accuracy, approximate the bounded final mask by a
finite-coordinate smooth one, choose a finite tree response, and only then
take the original matrix dimension limit. No intermediate polynomial `h`
is evaluated on the matrix, and no dimension-dependent recursion depth is
assumed.

## 4. Finite bivariate Hermite resolvent

Let `chi(g,z)=1_{|rho g+sz|≤alpha}`. With normalized Hermites, its coefficient
at `(k,l)`, `r=k+l`, is

\[
c_{k,l}=\beta_r\sqrt{\binom r l}\,\rho^k s^l,
\quad
\beta_0=p,\qquad
\beta_r=-\frac{2\phi(\alpha)H_{r-1}(\alpha)}{\sqrt{r!}}
\quad(r\ge2\text{ even}),
\]

and all odd-degree coefficients vanish. Fix an even cutoff `N` and `a>0`.
Remove the constant coefficient and set

\[
h_{k,l}=C\frac{c_{k,l}}{a+l},\qquad
2\le k+l\le N\text{ even},
\]

with `C` chosen for squared norm one. This gives exact mean zero and joint
evenness. Orthogonality gives

\[
Eh_Z^2=\frac{\sum_{k,l}l c_{k,l}^2/(a+l)^2}
                 {\sum_{k,l}c_{k,l}^2/(a+l)^2},\qquad
E h\chi=\frac{\sum_{k,l}c_{k,l}^2/(a+l)}
                  {\sqrt{\sum_{k,l}c_{k,l}^2/(a+l)^2}}.
                                                               \tag{5}
\]

The denominator `a+l`, with penalty on the innovation degree only, is
essential. The anchor degree `k` is not penalized by (1). Lagrange
multipliers for maximizing `E h chi` under the norm and conditional
Dirichlet constraints give precisely this resolvent form; the explicit
certificate uses (5) directly and does not need optimality.

For exact computation, group by `l` and define the nonnegative rational
numbers

\[
C_l=\sum_{\substack{2\le r\le N\\r\text{ even},\ r\ge l}}
 \frac{H_{r-1}(\alpha)^2}{r!}\binom r l
       \rho^{2(r-l)}(1-\rho^2)^l.                    \tag{6}
\]

Then `sum_k c_kl²=4phi(alpha)² C_l`. The common factor cancels from the
derivative ratio in (5), making that ratio an exact rational number when
`alpha,rho,a` are rational. The covariance becomes

\[
E h\chi=2\phi(\alpha)
 \frac{\sum_l C_l/(a+l)}{\sqrt{\sum_l C_l/(a+l)^2}}.
\]

There is no infinite Hermite tail: these formulas define the actual finite
bivariate polynomial used in the theorem.

## 5. Completed exact candidate and independent arithmetic audit

The candidate in `computations/fresh_anchored_fixed_point_certificate.py` is

\[
\alpha=73/100,\quad\rho=83/100,\quad a=9/2,\quad N=200.
\]

The exact script uses (6), rational outward intervals, and the already
independently audited extended Gaussian-CDF routine. Its completed output is

\[
\begin{split}
0.992417196066180527645430338491519097685719920657817216509478
\le D_Z\\
\le
0.992417196066180527645430338491519097685719920657817216509479<1,
\end{split}
\]

and

\[
\begin{split}
0.428376417656187756438676313934781469275119835457995857473025
\le L\\
\le
0.428376417656187756438676313934781469275119835457995857473259.
\end{split}
\]

Thus the certified bound is strictly greater than `4283/10000=0.4283`.
The conditional zero-degree mass is approximately 0.438816, and the
conditional derivative stability margin is approximately 0.0075828.
No infinite Hermite tail or numerical quadrature enters this certificate.

The algebra agent independently read the whole script and the whole
theorem. It verified the bivariate Hermite addition formula, exclusion of
the constant term, common-factor cancellation in the conditional derivative
ratio, covariance normalization, common-space contraction, and decorated
critical extension. The interval primitives and extended Gaussian-CDF tail
had already passed separate line-by-line audits in the preceding certificate.
No mathematical or arithmetic gap was found.

The director independently reconstructed the common-space contraction and
finite Hermite coefficient formulas, read the complete certificate, and
reran it at 21:26 UTC with the identical enclosing interval. This remains
an asymptotic analytic lower bound, not a convergence theorem.

## 6. Critical extension and its necessary exception

The strict criterion (1) is already sufficient for the candidate. For
completeness, an exact decorated-tree argument extends it to most critical
cases. Write

\[
h(G,Z)=\sum_{k,l}c_{k,l}
   \frac{H_k(G)H_l(Z)}{\sqrt{k!l!}},\quad
c_{00}=0,\quad k+l\text{ even},\quad\sum c_{k,l}^2=1.
\]

An innovation coefficient on a tree `T` is zero for the single-edge tree.
For every other tree, let `k` be its number of single-edge child branches
(fixed anchor leaves), and let the remaining child types have multiplicities
`m_tau`, with `l=sum m_tau`. The forced coefficient is

\[
u_T=c_{k,l}\sqrt{\frac{l!}{\prod_\tau m_\tau!}}
                     \prod_\tau u_\tau^{m_\tau}.
\]

The squared coefficients are decorated Galton--Watson tree probabilities:
a node has `l` innovation children and a mark of `k` anchor leaves, with
probability `c_kl²`. The innovation offspring law is
`Pr(N=l)=sum_k c_kl²`, of mean `D_Z`. The total finite-tree mass is its
extinction probability. Thus a unit innovation fixed point exists if
`D_Z≤1` and `Pr(N=0)>0`; the same coefficient recursion proves uniqueness.
For `D_Z>1`, or the critical deterministic-one-child case, the finite-tree
mass is less than one and no unit first-chaos solution exists.

The exception matters: `h(G,Z)=GZ` is jointly even, centered, normalized,
and has `D_Z=1`, but its offspring law is deterministic one-child and its
Gaussianization is a strict tree-index shift with no unit fixed point.
The resolvent candidates have positive zero-innovation-degree mass, so
their critical boundary is not this exceptional case.

## 7. Independent audit of the finite-ancestor extension

Later in the same session, the parent extended the single-anchor theorem
to 21 explicit ancestor-closed tree coordinates; see
`fresh_finite_anchor_fixed_point_2026_09_05.md` and its corresponding exact
certificate script. The following points were independently reconstructed.

For a finite ancestor-closed anchor set `A`, every inverse image `h_T` of
an anchor coordinate is a fixed normalized Hermite monomial of the anchor
vector. If `h(G,Z)` is orthogonal to these finitely many monomials, then
`U h(G,Z)` is orthogonal to every anchor coordinate for *every* unit
first-chaos innovation `Z` orthogonal to the anchors. Such an innovation
is jointly independent of all anchors, as is any correlated pair of
innovations. Thus the same conditional-Hermite contraction proof applies
on the complete orthogonal unit sphere, with Lipschitz constant
`sqrt(E h_Z²)`. There is no innovation-dependent constraint hidden in the
choice of removed monomials.

For `V=rho dot G+sZ`, the total degree-`d`, innovation-degree-`l` indicator
coefficient mass is exactly

\[
 \beta_d^2\binom dl\|\rho\|^{2(d-l)}s^{2l}.
\]

Removing the selected anchor inverse images affects only the `l=0`
group. For a top-child multiplicity vector `m`, its coefficient is

\[
 -2\phi(\alpha)\,
 \frac{He_{|m|-1}(\alpha)\prod_j\rho_j^{m_j}}
      {\sqrt{\prod_jm_j!}}.
\]

Hence the script's exact subtraction `raw²/prod(m_j!)`, after factoring
out `4phi(alpha)²`, is precisely correct. The edge constant is removed
separately by beginning the common-factor sum at total degree two. The
normalization, conditional derivative quotient, and covariance formula
all cancel this common factor correctly.

The actual 21-child list contains unique even child multisets with all
child indices strictly preceding the parent, so its finite closure and
orthonormal-monomial claims hold literally. The fixed degree 200 exceeds
every removed feature degree. An independent floating-point diagnostic,
using normalized Hermite recurrence rather than the exact script's raw
Hermite recurrence, gives

\[
 \|\rho\|^2=0.85835081,\quad
 D_Z=0.9908729106624282,\quad
 w=0.7004154308828989,
\]
\[
 J=0.43065817940552875.
\]

These decimals are a cross-check, not the arithmetic certificate. The
exact script targets `J>0.4306`; its interval primitives were previously
audited independently. No proof, dependence, coefficient-normalization,
or constant mismatch was found in this extension.
