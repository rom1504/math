# Mesoscopic row-correlation reduction

## Status

This note proves a quantitative bounded-degree reduction after deleting
an arbitrarily small fixed fraction of vertices.  This is the correct
two-limit statement for the asymptotic problem:
\[
n\to\infty,\qquad\text{then }\varepsilon,\eta\downarrow0.
\]
It is stronger and cleaner than trying to find one \(n-o(n)\) subset
whose degree bound is uniform in \(n\), which is not what the limiting
argument needs.

The remaining gap is no longer row-correlation concentration.  It is
the finite/local-type Boolean response theorem described at the end.

## 1. Normalization

Let \(A=A^\top\) be a zero-diagonal signing of order \(n\), and write
\[
Q(A)=\max_{x\in\{\pm1\}^n}|x^\top Ax|.
\]
Assume
\[
Q(A)\le(c_2+o(1))n^{3/2},
\qquad c_2=0.783387533648\ldots
\tag{1}
\]
in the doubled normalization.

For a principal signing \(A_R\) of order \(r\), put
\[
B_R=\frac{A_R}{\sqrt{r-1}},\qquad
C_R=B_R^2,\qquad
q_{ij}=(C_R)_{ij}.
\tag{2}
\]
Then \(C_R\succeq0\) and
\[
(C_R)_{ii}=1.
\tag{3}
\]

## 2. Small-fraction spectral removal

The Grothendieck--Pietsch factorization theorem gives the following
principal-submatrix estimate.  For every \(0<\varepsilon<1\), there
is \(R\subset[n]\), \(|R|=r\ge(1-\varepsilon)n\), such that
\[
\|A_R\|_{\rm op}
\le \frac{4K_GQ(A)}{\varepsilon n}.
\tag{4}
\]
Consequently, under (1),
\[
\boxed{
\|B_R\|_{\rm op}
\le
\frac{4K_G(c_2+o(1))}
{\varepsilon\sqrt{1-\varepsilon}},
}
\tag{5}
\]
and hence
\[
\boxed{
\|C_R\|_{\rm op}\le
L_\varepsilon+o(1),\qquad
L_\varepsilon=
\frac{(4K_Gc_2)^2}{\varepsilon^2(1-\varepsilon)}.
}
\tag{6}
\]
Using \(K_G<1.782214\),
\[
L_\varepsilon<
\frac{31.1884}{\varepsilon^2(1-\varepsilon)}.
\tag{7}
\]

For a fixed signing, passing to a principal submatrix cannot increase
the obstruction to a lower bound:
\[
Q(A)\ge Q(A_R).
\tag{8}
\]
Indeed, fix spins on \(R\) and average the full quadratic form over
independent spins on \(R^c\); its conditional mean is the quadratic
form on \(R\), so some extension has at least that absolute value.

Thus a lower bound proved on \(R\) loses only the normalization factor
\[
\left(\frac rn\right)^{3/2}\ge(1-\varepsilon)^{3/2},
\tag{9}
\]
which disappears when \(\varepsilon\downarrow0\).

## 3. Exact bounded-degree truncation

Since \(C_R\succeq0\), (3) and (6) imply, for every vertex \(i\),
\[
\sum_jq_{ij}^2
=(C_R^2)_{ii}
\le\|C_R\|_{\rm op}(C_R)_{ii}
\le L_\varepsilon+o(1).
\tag{10}
\]
Fix a fourth-moment error tolerance \(\eta>0\), and define
\[
\tau_{\varepsilon,\eta}
=\sqrt{\frac{\eta}{L_\varepsilon}}.
\tag{11}
\]
Let \(G_{\varepsilon,\eta}\) be the weighted graph on \(R\) containing
the pairs with
\[
|q_{ij}|\ge\tau_{\varepsilon,\eta}.
\]
Equation (10) gives the uniform degree bound
\[
\boxed{
\deg_{G_{\varepsilon,\eta}}(i)
\le
\frac{L_\varepsilon}{\tau_{\varepsilon,\eta}^2}
=\frac{L_\varepsilon^2}{\eta}
=:\Delta_{\varepsilon,\eta}.
}
\tag{12}
\]
Numerically,
\[
\Delta_{\varepsilon,\eta}
<
\frac{972.72}
{\varepsilon^4(1-\varepsilon)^2\eta}.
\tag{13}
\]

At the same time, all discarded correlations carry at most \(\eta\)
normalized fourth-moment mass:
\[
\boxed{
\frac1r
\sum_{\substack{i\ne j\\|q_{ij}|<\tau_{\varepsilon,\eta}}}
q_{ij}^4
\le
\tau_{\varepsilon,\eta}^2
\frac1r\sum_{i,j}q_{ij}^2
\le\eta+o(1).
}
\tag{14}
\]

Equations (12)--(14) are the desired quantitative mesoscopic
dichotomy.  There is no uncontrolled high-degree branch after the
small-fraction Pietsch removal: either a correlation survives in a
bounded-degree weighted graph, or its total relevance to the
rank-three/depth-two diagrams is at most \(\eta\).

## 4. Finite labels and local types

For fixed \(\varepsilon,\eta\), both
\(\Delta_{\varepsilon,\eta}\) and
\(\tau_{\varepsilon,\eta}\) are constants independent of \(n\).
Round each surviving \(q_{ij}\in[-1,1]\) to a mesh of width
\[
\rho<\frac{\zeta}{8\Delta_{\varepsilon,\eta}}.
\tag{15}
\]
Since \(|u^4-v^4|\le4|u-v|\) on \([-1,1]\), this changes the normalized
fourth-moment mass by at most \(\zeta/2\).  The resulting graph has:

- uniformly bounded degree;
- finitely many edge labels;
- compact rooted-neighborhood distributions.

Every sequence therefore has a subsequence converging locally weakly
to a rooted, finitely labelled graphing.  Any fixed-radius observable
depending continuously on those labels is approximable by finitely
many rooted types.

The explicit three-fibre family is a calibration, not an exception:
its correlation graph is a disjoint union of triangles, has degree
two, and has
\[
\Lambda_4\to2/81.
\]
It lies exactly inside the bounded-degree branch.

## 5. A one-subset mesoscopic version

For completeness, the spectral third-moment bound also gives a direct
\(n-o(n)\) statement, though with a growing degree bound.  The
mixed-sign Grothendieck inequality yields
\[
\operatorname{tr}(C^{3/2})=O(n),
\qquad
\|C\|_{\rm op}=O(\sqrt n).
\tag{16}
\]
Split \(C=C_{\le L}+C_{>L}\) spectrally.  Delete vertices with
\((C_{>L})_{ii}>\delta\).  The deleted fraction is
\[
O((\delta\sqrt L)^{-1}),
\tag{17}
\]
and on the remainder
\[
\frac1n\sum_{ij}(C_{>L})_{ij}^4
=O(\delta^2n^{1/4}).
\tag{18}
\]
Also
\[
\sum_j(C_{\le L})_{ij}^2\le L.
\tag{19}
\]
For example, taking
\[
L=n^{3/10},\qquad
\delta=n^{-7/50},\qquad
\tau=n^{-4/25}
\]
deletes \(O(n^{99/100})=o(n)\) vertices, discards \(o(n)\)
unnormalized fourth-moment mass, and leaves a threshold graph of
maximum degree
\[
O(n^{31/50}).
\tag{20}
\]
More generally the exponent can be brought down to
\(1/2+o(1)\).  This identifies the genuinely mesoscopic
\(\sqrt n\)-degree scale, but the fixed-\(\varepsilon\) reduction in
Sections 2--4 is the form useful for compactness.

## 6. What remains

The reduction above concerns the row-correlation kernel \(C_R=B_R^2\).
A depth-two cavity observable also remembers a compatible square-root
transport \(B_R\); the signs/phases of that transport are not
determined by \(C_R\) alone.

The remaining theorem needed for the universal \(c_2\) bound is:

> **Bounded-degree compatible-transport theorem.**  For every fixed
> \(\Delta\), finite edge-label set, and local weak limit of pairs
> \((B,C=B^2)\) whose \(C\)-correlation graph has degree at most
> \(\Delta\), the optimized dependent two-step Boolean response has
> energy at least \(c_2\).

The finite-fibre conditional-score conjecture in
`finite_type_conference_fibres.md` proves the correct candidate local
inequality when the graphing consists of finite components.  Extending
it to arbitrary bounded-degree rooted graphings is the precise next
step.

If that theorem holds, then (8)--(14) give
\[
\frac{Q(A)}{n^{3/2}}
\ge
(1-\varepsilon)^{3/2}(c_2-o_{\eta}(1))-o_n(1).
\]
Sending \(n\to\infty\), then \(\eta\downarrow0\), then
\(\varepsilon\downarrow0\), proves the universal depth-two lower
bound \(c_2\).

## 7. Compatible transport carried by cuts

The extra datum beyond \(C=B^2\) has an exact combinatorial form.
Let the scalar first residual have covariance kernel
\[
K_{ii}=s_t^2,\qquad K_{ij}=\kappa(q_{ij})\quad(i\ne j),
\]
where
\[
\kappa(q)=\sum_{\substack{\ell\ge3\\\ell\ {\rm odd}}}
u_\ell^2q^\ell.
\tag{21}
\]
For each row \(i\) of the flat square root \(B\), define a cut
\[
z^{(i)}_j=\sqrt{r-1}\,B_{ij}\in\{\pm1\}
\]
away from the negligible diagonal coordinate.  Then
\[
\boxed{
q_{jk}
=\frac1r\sum_i z^{(i)}_jz^{(i)}_k+o(1).
}
\tag{22}
\]
Thus the row-correlation graph comes with a representing probability
measure on its cuts.

The residual receiving field has covariance \(D=BKB\).  Its diagonal
at row \(i\) is
\[
\boxed{
D_{ii}
=s_t^2+
\frac2r\sum_{j<k}
\kappa(q_{jk})z^{(i)}_jz^{(i)}_k+o(1).
}
\tag{23}
\]
Averaging (23) over the representing cuts recovers the positive
orientation-even identity
\[
\frac1r\operatorname{tr}D
=s_t^2+
\frac2r\sum_{j<k}q_{jk}\kappa(q_{jk})+o(1)
\ge s_t^2.
\tag{24}
\]

Equations (22)--(23) are the exact compatible-transport information
that a rooted graphing must retain.  An arbitrary PSD weighted graph
with the same edge labels is not enough.

## 8. Graphing variational functional

For a finite approximation let
\[
S\sim\operatorname{Unif}\{\pm1\}^r,\qquad
G\sim N(0,C)
\]
with the usual local cavity coupling, put
\[
R_j=r_t(S_j,G_j),\qquad K=\mathbb ERR^\top,
\qquad W\sim N(0,D),\quad D=BKB,
\]
and use the scalar optimal moments
\[
(\alpha,\beta,\gamma,\delta)
=(0.5618430437,\ 0.5327670921,\
0.4238402492,\ 0.2179176762).
\]
The explicit compatible conditional score is
\[
\boxed{
L_B
=\beta CS+\alpha G
+\frac{\delta}{s_t}CR
+\gamma s_t\,C D^\dagger W,
\qquad
Y_B=\operatorname{sign}L_B.
}
\tag{25}
\]
The bounded-degree graphing conjecture is
\[
\boxed{
\liminf_{r\to\infty}
\frac1r\mathbb E Y_B^\top BY_B\ge c_2.
}
\tag{26}
\]
Fixed dither makes (25)--(26) continuous under local weak convergence
of the finitely labelled correlation graph together with its
cut-transport law.

## 9. Gauge invariance and a finite-girth lemma

Conjugating by a diagonal sign matrix sends
\[
C\mapsto\Sigma C\Sigma,\qquad
q_{ij}\mapsto\sigma_i\sigma_jq_{ij}
\]
without changing the Boolean optimization or the conditional-score
value.  Therefore every monomial in the fixed-dither Taylor expansion
must be Eulerian: each graph vertex has even total degree in the
multiset of correlation edges appearing in that monomial.

Consequences:

1. there is no linear correction at \(q=0\);
2. the quadratic Hessian is diagonal over edges and has the same
   coefficient as the two-vertex problem;
3. on a graph of girth \(g\), every orientation-sensitive correction
   of degree below \(g\) is a product of even powers of individual
   edges.

For the relaxed two-vertex covariance
\[
C(q)=
\begin{pmatrix}1&q\\q&1\end{pmatrix},
\qquad -1<q<1,
\]
direct high-precision evaluation of (25) gives an even kernel
\(\mathcal V_2(q)\) with
\[
\mathcal V_2(q)
=c_2+(0.170\ldots)q^2+O(q^4).
\tag{27}
\]
Representative computed values are:

| \(q\) | \(\mathcal V_2(q)\) |
|---:|---:|
| 0 | 0.783386 |
| 0.05 | 0.783830 |
| 0.10 | 0.785075 |
| 0.30 | 0.798681 |
| 0.50 | 0.82676 |
| 0.90 | 0.94950 |

These values are computed, not interval-certified.

The natural finite-girth statement is:

> **Finite-girth lemma.**  For every degree bound \(\Delta\), there
> are \(q_0(\Delta),c_0(\Delta)>0\) such that if the compatible
> correlation graph has girth tending to infinity and
> \(\max_e|q_e|\le q_0(\Delta)\), then
> \[
> \frac1r\mathbb E Y_B^\top BY_B
> \ge c_2+
> c_0(\Delta)\frac1r\sum_eq_e^2-o(1).
> \tag{28}
> \]

Gauge invariance proves the diagram support needed for (28), and
(27) identifies a positive edge Hessian.  What remains is a uniform
bound on the sum of higher doubled-edge diagrams.  Absolute
monotonicity in the variables \(q_e^2\) would prove (28) immediately,
but it has not yet been established.  Gaussian Anderson/Borell
comparison does not apply directly because the non-Gaussian source
pair \((S,R)\) is transported together with the receiving fields.

## 10. Audited fixed-dither Hessian

The first version of this section used the neighbor field
\(\sum_jq_{ij}(\beta S_j+\delta R_j)\).  Its asserted coefficient
\(2\mathbb E\psi_\tau'(U)\mathbb E(\beta S+\delta R)^2\), and the
uniform cluster bound that followed it, were not justified: covariance
derivatives involving the residual channel had been omitted.  Those
claims are withdrawn.

There is a simpler perturbation for which the complete quadratic
coefficient can be audited.  Fix \(\tau>0\), let
\[
\psi_\tau(u)=2\Phi(u/\tau)-1,
\tag{29}
\]
and take a *dithered stationary scalar response*
\[
y=\psi_\tau(U),\qquad
U=hS+aG+jR+mZ.
\tag{30}
\]
Here \(S\) is a sign, \(G,Z\) are independent standard Gaussians,
\(R\) is the normalized first residual, and
\[
a=\mathbb E Sy,\quad h=\mathbb E Gy,\quad
m=\mathbb E Ry,\quad j=\mathbb E Zy.
\tag{31}
\]
Thus the coefficients in (30) are precisely the reversed response
moments, as required by the scalar stationarity equation.  Put
\[
\ell_\tau=\mathbb E\psi_\tau'(U)>0
\tag{32}
\]
and
\[
k_\tau=-\mathbb E\!\left[U\psi_\tau''(U)\right]
=\frac{2}{\tau^3}
\mathbb E\!\left[U^2\phi(U/\tau)\right]>0.
\tag{33}
\]

Write \(C=I+Q\), where \(Q_{ii}=0\) and \(Q_{ij}=q_{ij}\).
For \(0<\lambda\le1\), use only neighboring *initial spins*:
\[
y_i^{(\lambda)}
=\psi_\tau\!\left(
U_i+\lambda h\sum_{j\ne i}q_{ij}S_j
\right).
\tag{34}
\]
This avoids every residual-covariance ambiguity because the \(S_j\)
are independent of one another and of the receiving variables at
site \(i\).

Let
\[
d_i=\sum_{j\ne i}q_{ij}^2.
\tag{35}
\]
At \(Q=0\), Taylor expansion of (34), followed by the sign averages,
gives the following response-matrix expansions through degree two:
\[
\begin{aligned}
\mathcal A_{ij}
&=\lambda h\ell_\tau q_{ij}+O_3
&& (i\ne j),\\
\mathcal A_{ii}
&=a+\frac{(\lambda h)^2d_i}{2}
  \mathbb E[S\psi_\tau''(U)]+O_3,\\
\mathcal H_{ii}
&=h+\frac{(\lambda h)^2d_i}{2}
  \mathbb E[G\psi_\tau''(U)]+O_3,\\
\mathcal M_{ii}
&=m+\frac{(\lambda h)^2d_i}{2}
  \mathbb E[R\psi_\tau''(U)]+O_3,\\
\mathcal J_{ii}
&=j+\frac{(\lambda h)^2d_i}{2}
  \mathbb E[Z\psi_\tau''(U)]+O_3.
\end{aligned}
\tag{36}
\]
Here \(O_3\) means terms of total degree at least three in the edge
variables.  There are no unlisted degree-two terms:

1. Gaussian integration by parts makes \(\mathcal H\) and
   \(\mathcal J\) diagonal through this order.
2. The residual has Hermite rank at least three across distinct
   sites, so \(\mathcal M_{ij}=O_3\) for \(i\ne j\).
3. Its covariance is
   \(K=s_t^2I+O_3\), and hence the normalized fresh field has
   covariance \(C+O_3\).
4. Mixed terms \(q_{ij}q_{ik}\) in the diagonal entries vanish
   because \(\mathbb E S_jS_k=0\) for \(j\ne k\).

Substitution into the exact paired-energy trace has two contributions.
The off-diagonal term in \(\mathcal A\), contracted with
\(C=I+Q\), contributes
\[
\frac{2\lambda h^2\ell_\tau}{r}
\sum_{i\ne j}q_{ij}^2.
\tag{37}
\]
The four diagonal shifts in (36) contribute
\[
\frac{(\lambda h)^2}{r}
\sum_i d_i
\left[
h\,\mathbb E(S\psi_\tau''(U))
a\,\mathbb E(G\psi_\tau''(U))
j\,\mathbb E(R\psi_\tau''(U))
m\,\mathbb E(Z\psi_\tau''(U))
\right].
\tag{38}
\]
By (30), the bracket is
\(\mathbb E[U\psi_\tau''(U)]=-k_\tau\).  Since
\(\sum_i d_i=\sum_{i\ne j}q_{ij}^2\), the full Hessian is therefore
\[
\boxed{
\mathcal E(Q)
=c_{2,\tau}
+h^2\left(2\lambda\ell_\tau-\lambda^2k_\tau\right)
\frac1r\sum_{i\ne j}q_{ij}^2
+O_3.
}
\tag{39}
\]
In particular, with
\[
\lambda_\tau=\min\{1,\ell_\tau/k_\tau\},
\tag{40}
\]
the coefficient has the explicit safe lower bound
\[
\boxed{
h^2\left(2\lambda_\tau\ell_\tau
-\lambda_\tau^2k_\tau\right)
\ge h^2\lambda_\tau\ell_\tau>0.
}
\tag{41}
\]
Equations (36)--(41), rather than the former (31)--(38), are the
certified perturbative conclusion.

## 11. Uniform operator-norm remainder

The dense transport does not in fact destroy the remainder estimate,
but it has to be controlled in operator norm rather than by declaring
it local.  The clean statement uses a fixed dither in *both* sign
steps.  Let the first sign be replaced by
\(\psi_\sigma(G+tS)\), regress away its \(S\)- and \(G\)-components,
and call the resulting smooth residual \(R_\sigma\).  At the end one
may send \(\sigma\downarrow0\), after the limit in \(r\).

Let
\[
\widehat Q=(|q_{ij}|)_{ij},\qquad
\varrho=\|\widehat Q\|_{\rm op}.
\tag{42}
\]
For a graph of maximum degree \(\Delta\),
\[
\varrho\le\Delta q_{\max}.
\tag{43}
\]
The following is the dimension-free Taylor estimate needed here.

> **Smooth covariance-response lemma.**  Fix
> \(\sigma,\tau>0\), the scalar coefficients in (30), and
> \(0<\varrho_*<1/4\).  There is a finite constant
> \(\mathfrak R_{\sigma,\tau,\varrho_*}\), independent of \(r\),
> of the compatible square root \(B\), and of the support graph, such
> that whenever
> \(\|\widehat Q\|_{\rm op}\le\varrho_*\),
> \[
> \left|
> \mathcal E(Q)-c_{2,\sigma,\tau}
> -h^2(2\lambda\ell_\tau-\lambda^2k_\tau)
>   \frac{\operatorname{tr}Q^2}{r}
> \right|
> \le
> \mathfrak R_{\sigma,\tau,\varrho_*}\,
> \|\widehat Q\|_{\rm op}
> \frac{\operatorname{tr}Q^2}{r}.
> \tag{44}
> \]

Here \(\mathcal E(Q)\) is the exact paired response functional, not a
truncated diagram series.  Notice that
\(\operatorname{tr}Q^2=\sum_{i\ne j}q_{ij}^2\).

### Proof of the smooth covariance-response lemma

Interpolate
\[
C_u=I+uQ,\qquad 0\le u\le1.
\tag{45}
\]
If \(B^2=C_1\), write
\[
B=O\,C_1^{1/2},\qquad O=\operatorname{sgn}(B),
\]
and use
\[
B_u=O\,C_u^{1/2}.
\tag{46}
\]
The matrices \(O,C_u,Q\) commute, \(O\) is orthogonal, and
\[
\|C_u^{\pm1/2}\|_{\rm op}\le(1-\varrho_*)^{-1/2},\qquad
\|B_u^{(k)}\|_{\rm op}
\le c_k(1-\varrho_*)^{1/2-k}\|Q\|_{\rm op}^k
\tag{47}
\]
for \(k\le3\).

The cross-site residual covariance has the entrywise form
\[
K_u=s_{\sigma}^2I+\kappa_\sigma[uQ],
\qquad
|\kappa_\sigma(q)|\le s_\sigma^2|q|^3.
\tag{48}
\]
The brackets mean entrywise application off the diagonal.  Hence
\[
\|\kappa_\sigma[uQ]\|_{\rm op}
\le s_\sigma^2 q_{\max}^2
\|\widehat Q\|_{\rm op}
\le s_\sigma^2\varrho_*^3.
\tag{49}
\]
In particular \(K_u\) is bounded above and below uniformly.  For the
fresh field,
\[
D_u=B_uK_uB_u,
\tag{50}
\]
so
\[
\|D_u-s_\sigma^2C_u\|_{\rm op}
\le
\|B_u\|_{\rm op}^2
\|\kappa_\sigma[uQ]\|_{\rm op}
\le
(1+\varrho_*)s_\sigma^2\varrho_*^3.
\tag{51}
\]
This is the step that controls the dense compatible transport.  No
sparsity of \(B\) is asserted or needed.

Differentiate the exact response traces three times along (45).
Gaussian covariance differentiation inserts one copy of \(Q\) and
two one-site derivatives.  Differentiating (46), the explicit
neighbor field in (34), or (48) also inserts edge weights from \(Q\).
Because \(\psi_\sigma,\psi_\tau\) have bounded derivatives of every
fixed order, all one-site factors produced by three differentiations
have a finite bound depending only on \(\sigma,\tau\) and the scalar
coefficients.  Equations (47)--(51) bound every intervening dense
matrix in operator norm.

It remains to check that the index sums are dimension-free.  Gauge
invariance makes each resulting edge multigraph Eulerian.  After two
edge factors are placed in Hilbert--Schmidt norm, every remaining
edge factor is bounded by \(\widehat Q\) in operator norm.  Thus each
third-or-higher connected contraction is bounded by
\[
\|\widehat Q\|_{\rm op}
\|Q\|_F^2.
\tag{52}
\]
The same estimate applies to (50): the two \(B_u\)'s cost only the
operator-norm factors in (47), while the first nonconstant term of
\(\kappa_\sigma\) already contains three edge factors.  Dividing the
trace by \(r\), the third directional derivative therefore obeys
\[
\sup_{0\le u\le1}
\left|\mathcal E'''(u)\right|
\le
6\mathfrak R_{\sigma,\tau,\varrho_*}\,
\|\widehat Q\|_{\rm op}
\frac{\|Q\|_F^2}{r}.
\tag{53}
\]
Taylor's formula with integral remainder, together with the audited
zeroth, first, and second derivatives in Section 10, proves (44).

For complete explicitness, the constant in (44) can be defined
without an asymptotic convention as
\[
\mathfrak R_{\sigma,\tau,\varrho_*}
=\frac16\sup
\frac{r|\mathcal E'''(u)|}
{\|\widehat Q\|_{\rm op}\|Q\|_F^2},
\tag{54}
\]
where the supremum is over \(r\), \(0\le u\le1\), compatible
instances satisfying (42), and the fixed smoothed response above.
The derivative calculation (47)--(53) proves that this supremum is
finite.  It is effectively computable from the derivatives through
order six of \(\psi_\sigma,\psi_\tau\); no compactness in the matrix
dimension is being assumed in (54).

## 12. Rigorous small-edge radius

Take \(\lambda=\lambda_\tau\) from (40) and set
\[
a_{\sigma,\tau}^{\rm safe}
=h^2\lambda_\tau\ell_\tau>0.
\tag{55}
\]
Choose
\[
\boxed{
\rho_{\sigma,\tau}
=\min\left\{
\frac18,\,
\frac{a_{\sigma,\tau}^{\rm safe}}
{2\mathfrak R_{\sigma,\tau,1/8}}
\right\}.
}
\tag{56}
\]
If
\[
\boxed{
\Delta q_{\max}\le\rho_{\sigma,\tau},
}
\tag{57}
\]
then (43)--(44) give the uniform estimate
\[
\boxed{
\mathcal E(Q)
\ge c_{2,\sigma,\tau}
+\frac{a_{\sigma,\tau}^{\rm safe}}2
\frac1r\sum_{i\ne j}q_{ij}^2.
}
\tag{58}
\]
This proves the small-edge theorem with a radius uniform in the
number of vertices.  The radius can shrink when the two dithers are
removed; that is harmless in the order
\[
r\to\infty,\qquad
\text{then the correlation cutoff}\downarrow0,\qquad
\text{then }\sigma,\tau\downarrow0.
\]

What remains in the bounded-degree compatible-transport theorem is
therefore the graph formed by edges above the fixed threshold
\(\rho_{\sigma,\tau}/\Delta\), together with their interaction with
the perturbative background.  Treating those large connected
components still requires the finite-type/graphing variational
inequality; (58) does not by itself justify conditioning them away.
