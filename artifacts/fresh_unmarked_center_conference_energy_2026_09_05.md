# Unmarked-center response: an exact conference energy law

Date: 2026-09-05. Status: independent derivation; root independently
reconstructed the full derivative bounds, endpoint terms, cavity covariance,
fourth-order zero term, and factor two, and provisionally passed the proof.
This is a conference-specific theorem. No universality assertion under only
the low quadratic-cap hypothesis is made.

Final scope update: the full conference-energy identity remains a
restricted theorem, independently audited by the director and variational
track. It was a precursor to, and is not a required premise of, the later
actual-sign weighted projection theorem
`fresh_full_weighted_unmarked_projection_2026_09_05.md` and its all-odd
extension. Those later proofs establish a different, specific weighted
identity under any fixed operator cap and support the unrestricted
normalized-gain theorem. They do not assert that the full nonlinear
conference formula (1) holds for every bounded-op signing. The unresolved
general-matrix discussion in Section 5 is historical in this precise sense.

## 1. Statement and normalization

Let `A=A^T` be a hollow sign matrix of order `n`, put `m=n-1` and
`B=A/sqrt(m)`, and assume the exact conference identity `B^2=I`.
Let `S` have independent uniform sign coordinates. Choose a fixed odd
`h in C_b^infinity(R)` such that

\[
 E h(N)^2=1,\qquad E[N h(N)]=E h'(N)=0,
 \qquad N\sim N(0,1).
\]

Define `G=BS`, `Z=B h(G)` coordinatewise. Let `F,H` be fixed smooth bounded
functions on `R^2` with all derivatives bounded, `F` jointly odd and `H`
jointly even. Set

\[
 f_i=F(G_i,Z_i)+S_iH(G_i,Z_i).
\]

For independent standard Gaussians `g,z`, the claimed limit is

\[
\boxed{
 \frac1{2n}E f^T Bf\ \longrightarrow\
 E H(g,z)\,E[gF(g,z)]
 +E[zF(g,z)]\,E[h(g)F(g,z)].
}                                                    \tag{1}
\]

In particular, if `|F|+|H|<=1`, the right side is a valid asymptotic lower
bound for `q(A)/(n sqrt(m))`. The second product is a genuinely unoriented
self-energy term; it is absent from the earlier all-marked odd-degree tree
basis. Taking `H=0` gives its standalone energy identity.

All functions and approximation parameters must be fixed before `n` tends
to infinity. This statement does not assert the same limit for arbitrary
discontinuous hard rules without an additional continuity argument.

## 2. Why the new field is outside the earlier tree family

Formally choosing `h=He_3/sqrt(6)` gives

\[
 Z=\frac1{\sqrt6} B\,\mathrm{He}_3(BS),\qquad
 BZ=\frac1{\sqrt6}\mathrm{He}_3(G).                    \tag{2}
\]

Its leading injective diagram has an unmarked internal center of degree
four and three spin-marked leaves. The earlier tree theorem marks every
free vertex and requires all free degrees to be odd. Equation (2) is an
exact matrix identity, not a Gaussian approximation. Under `A -> -A`, this
new cubic field is unchanged whereas the old marked fields change sign.

For the cubic Gaussian tensor `T_i=sum_j B_ij b_j^(tensor 3)`, orthogonality
of the rows `b_j` gives unit norm and every proper contraction has squared
norm `sum_j B_ij^4=1/m`. Its repeated-index coefficients satisfy

\[
 T_i(a,a,b)=\frac{\delta_{ib}-B_{ia}B_{ab}}m.
\]

Thus repeated-index mass is `O(1/n)` and the multilinear cubic has maximum
influence `O(1/n)`. On the sign cube the exact own-spin central derivative
of the full cubic field is

\[
 D_i Z_i=-\frac{3G_i^2-1}{\sqrt6\,m}.
\]

Consequently the own spin asymptotically decouples. These observations give
a separate route to its local Gaussian law, but a local law alone would
not prove (1): the energy sums `n^2` small correlations. Sections 3--4
instead compare the complete energy and then compute it.

## 3. Whole-energy Rademacher-to-Gaussian replacement

We supply quantitative bounds to avoid importing a generic AMP theorem.
Replace input coordinates one at a time by independent standard normals.
At an intermediate replacement fix coordinate `k` at a real value `t`;
all other coordinates are independent symmetric signs or normals. Write
`b=B e_k`. Differentiation with respect to this coordinate gives

\[
 G'=b,\quad Z'=B[h'(G)b],\quad
 Z''=B[h''(G)b^2],\quad Z'''=B[h'''(G)b^3].           \tag{3}
\]

Products and powers inside brackets are coordinatewise. Orthogonality and
flatness give deterministic bounds

\[
 \|Z'\|_2\le C,\quad \|Z''\|_2\le Cn^{-1/2},\quad
 \|Z'''\|_2\le Cn^{-1}.                             \tag{4}
\]

The needed coordinatewise moment improvements are, for every fixed finite
`p`, uniformly in all hybrid choices,

\[
 \max_i\|Z_i'\|_{L^p}\le C_p(1+t^2)n^{-1/2},\qquad
 \max_i\|Z_i''\|_{L^p}\le C_p(1+|t|)n^{-1}.         \tag{5}
\]

Here is an elementary verification. For `W=sum_j w_j r(G_j)`, its
continuous gradient has squared norm at most
`||r'||_infinity^2 ||w||_2^2`, by orthogonality. On sign coordinates the
central derivative differs from the corresponding gradient at the actual
field by at most `C sum_j |w_j| B_jl^2`, which is at most
`C ||w||_1/m` for each coordinate `l`. Therefore the total product-space
gradient square is bounded pointwise by

\[
 C\bigl(\|w\|_2^2+n\|w\|_1^2/m^2\bigr).           \tag{6}
\]

The tensorized Gaussian/sign Poincare inequality, followed by its usual
power induction, bounds every fixed centered `L^p` norm by a constant
times the square root of (6). One can perform the induction using even
powers and `|D W^r|<=r |DW| max(|W_+|,|W_-|)^(r-1)`; all fixed moments
exist. The gradient bound is pointwise, so its use introduces no
dimension-dependent constant.

For `Z_i'`, take `w_j=B_ij B_jk`, `r=h'`: `||w||_1<=1` and
`||w||_2^2<=1/m`. The mean of `h'(G_j)` is `O((1+t^2)/n)`.
Indeed, excluding coordinate `k`, ordinary fourth-order one-dimensional
Lindeberg replacement costs `O(sum_l B_jl^4)=O(1/n)`; the corresponding
Gaussian variance differs from one by `O(1/n)`. The fixed shift `B_jk t`
has zero first-order mean term since `h''` is odd and the unshifted sum is
symmetric. Finally `E h'(N)=0` by hypothesis. This proves the first bound
in (5).

For `Z_i''`, take `w_j=B_ij B_jk^2`, `r=h''`: `||w||_1=O(n^-1/2)` and
`||w||_2^2=O(n^-2)`. Symmetry and oddness of `h''` give
`|E h''(G_j)|<=C |t|/sqrt(n)`. Equation (6) now proves the second bound.

Let `f(s)` denote the response vector for this hybrid input. The chain
rule, (4)--(5), and Holder give, for every fixed finite `p`,

\[
 \|\,\|f\|_2\,\|_p=O((1+|t|)\sqrt n),\quad
 \|\,\|f'\|_2\,\|_p=O((1+|t|)^C\sqrt{\log(n+2)}),
\]
\[
 \|\,\|f''\|_2\,\|_p
 =O((1+|t|)^C n^{-1/2}\sqrt{\log(n+2)}),\quad
 \|\,\|f'''\|_2\,\|_p
 =O((1+|t|)^C n^{-1}\sqrt{\log(n+2)}).               \tag{7}
\]

The harmless logarithm bounds `max_i |s_i|` in `L^p`. For example
`sum_i E|Z_i'|^6=O(n^-2)` and
`sum_i E|Z_i' Z_i''|^2=O(n^-2)`. The explicit derivative of the factor
`s_i H` is confined to `i=k`; since `B_kk=0`, its contributions to the
second and third derivatives are respectively `O_Lp(n^-1/2)` and
`O_Lp(n^-1)`. This endpoint check is essential.

For `Ecal(s)=f(s)^T Bf(s)/n`, symmetry gives

\[
 Ecal'''=\frac2n\{f'''{}^T Bf+3f''{}^T Bf'\}.
\]

Since `||B||op=1`, (7) yields
`E|Ecal'''|<=C(1+|t|)^C n^-3/2 log(n+2)`.
A second-order Taylor replacement uses only equality of the first two
input moments. Integrating the remainder along the segment from zero to
the input coordinate costs the same bound with a fixed finite Gaussian
moment. Summing over all coordinates gives

\[
 |E Ecal(S)-E Ecal(N^{(n)})|
       =O(n^{-1/2}\log(n+2))=o(1).                  \tag{8}
\]

Thus the full energy, not merely individual fields, has been replaced.

## 4. Exact Gaussian-input cavity calculation

For Gaussian input, `G=BS` has independent standard Gaussian coordinates
and `S=BG`. Put `R(g)=(g,h(g))` and

\[
 U_i=(S_i,Z_i)=\sum_k B_{ik}R(G_k),\qquad
 \mathfrak f(g,s,z)=F(g,z)+sH(g,z).
\]

The vector `R(N)` is centered, has identity covariance, and is centrally
symmetric. Also `U_i` is independent of `G_i` because `B_ii=0`.

For `i!=j`, remove the two endpoints and write

\[
 U_i=U_i^{ij}+B_{ij}R(G_j),\qquad
 U_j=U_j^{ij}+B_{ij}R(G_i).
\]

The cavity vectors have covariance `(1-B_ij^2) I_2` separately and zero
cross-covariance, exactly: `sum_(k!=i,j) B_ik B_jk=0`.
Taylor expansion in the two endpoint increments has expected remainder
`O(B_ij^2)=O(1/n)`. Constants are uniform because all cavity moments of
fixed order are bounded and `mathfrak f` grows only linearly in `s`.

The zero-order term is `O(1/n)`, not merely `o(1)`. After integrating the
two endpoint Gaussians it becomes
`E barf(U_i^{ij}) barf(U_j^{ij})`, where `barf` is jointly odd. Replace
each independent vector summand `R(G_k)` by a two-dimensional standard
Gaussian. Their first three joint moments agree (central symmetry kills
the third moments). Fourth-order Lindeberg replacement costs
`C sum_k (|B_ik|+|B_jk|)^4=O(1/n)`. The replacing cavity vectors are
independent centered Gaussians, so the zero-order expectation is exactly
zero after replacement.

For the first-order coefficient the same replacement, or its elementary
central limit consequence, gives

\[
 E f_i f_j
 =2B_{ij}\sum_{a=0}^1
  E[\partial_{u_a}\mathfrak f(g,U)]
  E[R_a(g)\mathfrak f(g,U)]+O(1/n),                 \tag{9}
\]

where `g` and the two coordinates of `U` are independent standard
Gaussians. Replacing covariance `1-B_ij^2` by one costs another `O(1/n)`.
The four expectations in (9) are respectively

\[
 EH,\quad E[gF],\quad E[\partial_z F]=E[zF],\quad E[h(g)F].
\]

Finally `sum_(i,j) B_ij^2=n`. The accumulated error after multiplying
(9) by `B_ij/(2n)` is `O(n^-1/2)`. This proves (1), in conjunction with
(8).

## 5. Exact general-matrix defects: no universal conclusion yet

This section concerns extending the full formula (1), not the later
weighted projection theorem, whose broader hypotheses were proved
separately after this note.

For arbitrary hollow symmetric sign `A`, put `Q=B^2`; then `Q_ii=1` but
typically `Q!=I`. For Gaussian input, `G` has covariance `Q`. If the odd
function `h` has normalized Hermite expansion `h=sum_(r odd>=3) u_r h_r`,
then the following are exact:

\[
 K_h(Q)_{ij}=E h(G_i)h(G_j)=\sum_{r\ge3\text{ odd}}u_r^2 Q_{ij}^r,
\]
\[
 \operatorname{Cov}(Z)=B K_h(Q)B,\qquad
 \operatorname{Cov}(G,Z)=0,\qquad BZ=Qh(G),\qquad BG=QS.
                                                               \tag{10}
\]

In particular

\[
 \frac1n E\|Z\|_2^2
 =\sum_{r\ge3\text{ odd}}u_r^2\frac1n\sum_{i,j}Q_{ij}^{r+1}
 \ge1.                                                        \tag{11}
\]

For the cubic this is exactly `1+n^-1 sum_(i!=j) Q_ij^4`.
The nonnegative variance excess in (11) is not an energy inequality.
It neither guarantees asymptotic Gaussianity of `Z`, nor gives a sign to
the extra terms involving `(Q-I)h(G)`.

The proof above uses `B^2=I` in three genuinely different places: the
uniform operator bound in whole-energy replacement; independence of the
Gaussian coordinates `G`; and the exact independent cavity covariance.
The low-cap condition `q(A)=O(n^1.5)` only supplies operator bounds growing
with `n` and does not by itself validate any of these replacements.
Therefore (1) must not be advertised as a universal lower-bound theorem.

The whole-energy replacement in Section 3 does extend, with constants
depending on `L`, to any flat sign matrix with `||B||op<=L` uniformly.
Indeed, the deterministic derivative bounds and the product-gradient
bound acquire only fixed powers of `L`; the scalar mean cancellations use
flat rows and parity, not orthogonality. Thus under a bounded operator norm
the unresolved issue is the Gaussian-input calculation alone.

If additionally `rho=max_(i!=j)|Q_ij| -> 0`, then

\[
 \|K_h(Q)-I\|_{op}
 \le \max_i\sum_{j\ne i}\sum_{r\ge3\text{ odd}}u_r^2|Q_{ij}|^r
 \le \rho\max_i\sum_j Q_{ij}^2\le L^4\rho=o(1).      \tag{12}
\]

Consequently the following linear-feature energies are proved, with
`tau_3=Tr(B^3)/n` (which need not itself converge):

\[
 \frac1{2n}E G^T BG=\tau_3/2,\qquad
 \frac1{2n}E Z^T BZ=\tau_3/2+o(1),\qquad
 \frac1{2n}E h(G)^T B h(G)=o(1),
\]
\[
 \frac1n E Z^T B h(G)=1+o(1),\qquad
 E G^T BZ=E G^T B h(G)=0.                            \tag{13}
\]

These are Gaussian-input identities followed by (12); (8) transfers the
corresponding fixed smooth bounded-response energies when applicable.
Equations (12)--(13) suggest, but do not prove, a nonlinear formula with
an additional term `tau_3 ((E gF)^2+(E zF)^2)/2`. A pointwise multivariate
central limit theorem is not sufficient to justify that formula. No such
formula is banked here.
