# General masks after a diagonal gauge: independent proof audit

Date: 2026-09-05. Auditor: fresh-limit literature subagent.

## Verdict

**Passed.** A preselected diagonal gauge removes the constant Schur-kernel
obstruction to the uniform rooted Gaussian theorem for arbitrary, rather
than even, smooth masks. The extension also supports finitely many masks and
smooth final functions of all resulting fields. The direct-energy identity
below follows by the same genuinely collective indirect-path cancellation
as in the independently audited second-response proof.

The associated finite-first-layer **paired/direct-energy certificate** has
a universal ceiling strictly below `1/2`. This does not bound every actual
quadratic energy produced by the rounding algorithm, and it does not solve
convergence of the original minima.

All functions and the number of fields below are fixed before dimension
tends to infinity. Smooth masks have bounded derivatives of all needed
orders. Any removal of smoothing is performed after the dimension limit.

## 1. Choosing the gauge, and what it preserves

Start with a hollow symmetric signing `A`, put `m=n-1`, `B=A/sqrt(m)`, and
assume the fixed cap `q(A)<=C n^(3/2)`. The existing bootstrap gives
`||B²||op=O_C(sqrt(n))` and `Tr B⁴=O_C(n^(3/2))`.

For independent gauge signs `epsilon`, each `(B epsilon)_i` is subgaussian
with variance proxy `sum_j Bij²=1`. Hence

\[
P\{\max_i|(B\epsilon)_i|>\sqrt{2\log(4n)}\}
\le 2n e^{-\log(4n)}=1/2.
\]

Choose any successful gauge deterministically and set
`D=diag(epsilon)`, `B'=DBD`. Then

\[
\max_i |(B'\mathbf1)_i|^2\le2\log(4n).
\]

The Boolean quadratic optimum, symmetry, hollowness, entry magnitudes,
operator norm, and traces of even powers are preserved. The Gram matrix is
**conjugated**, `Q'=DQD`, not literally unchanged; its diagonal, spectrum,
and absolute entry values are preserved. Below rename `B'` as `B`.

The probe signs `S` are fresh and independent of this gauge choice. Choosing
the gauge using the same probe signs would not justify the subsequent
independence arguments. A single successful gauge works simultaneously for
every fixed finite family of masks, because its defining event contains no
response function.

## 2. Transport of arbitrary smooth row functions

For any fixed smooth `g`, endpoint extraction remains exact:

\[
E[S_jg(G_j)S_kg(G_k)]
=m^{-1}K_g(Q_{jk})+O_g(m^{-3/2}),\quad j\ne k,
\]

where `G=BS` and `K_g(q)=E[g'(Z)g'(Z')]` for unit Gaussian correlation `q`.
The earlier factored Taylor/Lindeberg argument has no parity hypothesis.
The only change is the Hermite expansion

\[
K_g(Q)=c_0J+\sum_{r\ge1}c_r Q^{\circ r},\quad
c_0=(Eg'(Z))^2,\quad c_r\ge0,\quad\sum_{r\ge0}c_r=Eg'(Z)^2.
\]

All nonconstant Schur powers have operator norm at most `||Q||op`. The
constant piece must be handled explicitly, rather than included in that
bound. Sandwiching by row `B_i` gives its contribution
`c0 (B1)_i²/m=O_g(log(n)/n)`. Thus, uniformly in `i` and the permitted
matrices,

\[
E\{B[Sg(BS)]\}_i^2=Eg(Z)^2+o(1).
\]

This applies to arbitrary fixed polynomial remainders, and also to `h''`
when `h` is not even. In fact, for a Hermite remainder after retaining degree
one, `E(g')=E[Zg]=0` exactly, so its constant-kernel obstruction disappears
even without a gauge. No derivative-norm approximation as the truncation
degree grows is needed: dimension tends to infinity first for each fixed
remainder.

## 3. Every fixed rooted Hermite order, including order one

For every integer `r>=1`, define

\[
F_{i,r}=\sum_jB_{ij}Z_jH_r(B_jZ)
=I_{r+1}\!\left(\operatorname{Sym}
       \sum_jB_{ij}e_j\otimes B_j^{\otimes r}\right).
\]

Hollowness makes this a pure chaos of degree `r+1`, regardless of the parity
of `r`. Its exact variance is

\[
r!\left[1+\frac r m
 \left(B_iQ^{\circ(r-1)}B_i^\top-1\right)\right].
\]

For `r>=2`, the nonconstant Schur bound gives limit `r!`. For `r=1`, interpret
`Q^(circ 0)=J`; the exact variance is

\[
1+\frac{(B\mathbf1)_i^2-1}{m}\longrightarrow1.
\]

The four previously audited contraction cases do not use parity. For
`1<=ell<=r`, their squared norms are bounded respectively by
`Tr Q²/m²`, `Tr Q²/m²`, `||Q||op/m`, and `||Q||op²/m²`.
At `r=1` only the first three cases occur, with residual branch exponent
one; all still vanish. Thus the same primary fixed-chaos Gaussian theorem
applies at every degree `r+1>=2`.

The imported statements are Theorem 1 of
[Nualart--Peccati](https://arxiv.org/pdf/math/0503598) and Theorem 7, with
setup (9), of
[Nualart--Ortiz-Latorre](https://arxiv.org/pdf/math/0703240). Their actual
hypotheses were read in the preceding audit: neither imposes a parity
restriction, and the multivariate theorem allows first chaos and repeated
chaos orders.

The exact discrete Hermite recurrence, bounded coefficient/root count,
low-influence Rademacher replacement, and repeated-index Gaussian Wick
comparison likewise have no parity restriction. For example,
`H3(Gj)=Wj,3-(2/m)Wj,1`; the transported lower-order term is uniformly L²
bounded. The linear field `Gi` is included as a first-chaos component.

Finite collections of rooted Hermite orders therefore converge jointly to
independent Gaussians of variances `r!`, independent of the linear field.
Their distinct chaos degrees are orthogonal. Tail transport then proves the
uniform arbitrary-mask rooted Gaussian law.

## 4. Finite families and coordinate-pair limits

For fixed masks `h_1,...,h_k`, put

\[
Y_i^a=\{B[S h_a(BS)]\}_i,\qquad
U_i=(G_i,Y_i^1,\ldots,Y_i^k).
\]

Uniformly in the coordinate, `U_i` converges to a centered Gaussian vector
`U=(G,Y_1,...,Y_k)` with

\[
EG^2=1,\quad EGY_a=Eh_a(Z),\quad
EY_aY_b=E[h_a(Z)h_b(Z)].
\]

Singular covariance matrices are allowed: use the independent normalized
Hermite-chaos coordinates first and then take the required linear
combinations. No covariance inversion is needed for this limit theorem.

At two coordinates, the order-one rooted covariance is exactly

\[
\operatorname{Cov}(F_{i,1},F_{j,1})
=Q_{ij}+\frac{(B\mathbf1)_i(B\mathbf1)_j-Q_{ij}}m.
\]

Its extra term is uniformly `o(1)` after the gauge. Higher rooted orders have
covariance `r!Qij+O_r(||Q||op/m)`. Hence whenever `Qij->0`, the entire
coordinate blocks `U_i,U_j` converge to independent copies of `U`.
The existing averaged-pair argument applies because
`Tr Q²/n²->0`; conjugating the Gram matrix did not change its absolute
entry values.

## 5. General joint-final-function energy identity

Let `F,H` be fixed bounded C³ functions on `R^(k+1)` with bounded derivatives.
No parity hypothesis on these functions is necessary. Write `partial_0`
for the derivative in `G`, and `partial_a` for the derivative in `Y_a`.
Then

\[
\boxed{\begin{aligned}
\frac1n E F(U)^\top B[S H(U)]\ \longrightarrow\quad
&E[\partial_0F(U)]\,E[H(U)]\\
&+\sum_{a=1}^k E[\partial_aF(U)]\,
                       E[h_a(G)H(U)].
\end{aligned}}
\]

On the left, `F(U)` means the vector with entries `F(U_i)`; on the right,
`U` denotes the limiting single-coordinate Gaussian vector.

Here is the new bookkeeping needed beyond the scalar proof. The direct
increment of `Gi` under `Sj` is exactly `Bij`. The increment of `Y_i^a` is
`Bij h_a(Gj)+Tij^a+O_L2(1/m)`. The restoration estimate is valid because
arbitrary-function transport now applies to `h_a''`. The cube-gradient and
eighth-moment estimates for every `Tij^a` did not require parity, so

\[
\sup_j E\max_i\|D_jU_i\|=o(1),\qquad
\sup_{i,j} E\|D_jU_i\|^2=O(1/m).
\]

These statements use fixed `k`. The multivariate endpoint Taylor error is
still `O(1/m)` per pair and therefore `o(n)` in total. Every off-diagonal
indirect contribution has exactly the old form with `F'` replaced by
`partial_a F`; the bilinear cap bounds its endpoint variation by
`C beta(B)/m` times the maximal increment of `U`. The own-coordinate
terms split into one small common `Z0^a` and one removed diagonal term,
exactly as before. There is no derivative of `H` in the own `G` direction,
because `D_jG_j=Bjj=0`.

The surviving direct term factors by the coordinate-pair Gaussian limit.
This produces precisely the boxed formula; no generic AMP or growing-depth
universality claim has been imported.

For Boolean conversion it is sufficient that `|F(u)|+|H(u)|<=1` everywhere:
the two means `sigma F(U_i)+S_i H(U_i)` then lie in the cube. Conditional
independent rounding and the half difference of the two oriented energies
give a valid lower bound on `q(B)`. Alternatively the original scalar
threshold-channel construction enforces this cube condition automatically.
The resulting signs transfer back to the original matrix by multiplication
by `D`, preserving their energies.

## 6. A bounded counterexample if the gauge is omitted from the uniform law

Take a symmetric Sylvester Hadamard matrix with all-positive first row,
remove its diagonal, and normalize by `sqrt(n-1)`. This is a genuine low-cap
family with bounded operator norm, but its first normalized row sums to
`sqrt(m)`. Therefore its degree-one rooted chaos has variance tending to
two, not one.

This defect also occurs for an admissible bounded mask: take `h(z)=sin z`.
The constant Schur-kernel coefficient is `(E cos Z)²=e^-1`. At the first
coordinate the rooted field has asymptotic variance

\[
E\sin^2Z+e^{-1}=\frac{1-e^{-2}}2+e^{-1},
\]

not `E sin²Z`. Indeed the same finite-chaos contraction proof gives a
Gaussian limit with this augmented variance; only the quadratic chaos has
its variance doubled. The Hermite tails after degree one have zero mean
derivative, so the transported L² approximation remains valid there.
Thus the gauge addresses an actual obstruction to a uniform-coordinate
arbitrary-mask statement, not merely a proof convenience.

## 7. Ceiling for the finite-first-layer paired certificate

Adjoin `h_0=1`, put `Sigma=(<h_a,h_b>)_(a,b=0)^k` in Gaussian L², and let
`V=(G,Y_1,...,Y_k)` have covariance `Sigma`. Define

\[
k_0(g)=E[H(V)\mid G=g],\qquad b_a=\langle h_a,k_0\rangle.
\]

Then `b` lies in the range of `Sigma`: every null relation among the `h_a`
is also orthogonal to `k_0`. Set

\[
W=(\Sigma^+b)\cdot V.
\]

Gaussian integration by parts, valid also for this singular covariance,
identifies the direct-energy certificate as

\[
(E\nabla F)\cdot b=E[WF(V)].
\]

Moreover

\[
EW^2=b^\top\Sigma^+b
=\|\operatorname{Proj}_{\operatorname{span}\{h_a\}}k_0\|_2^2
\le Ek_0(G)^2\le EH(V)^2\le u:=E|H(V)|.
\]

This argument allows signed `H`. Since `W` is centered Gaussian and
`|F|<=1-|H|`, central rearrangement at fixed `E|H|=u` gives

\[
E[WF]\le E[|W|(1-|H|)]
\le2\sqrt u\,\phi\!\left(\Phi^{-1}((1+u)/2)\right).
\]

The endpoint cases `u=0,1` have value zero. Writing
`u=2Phi(t)-1<=sqrt(2/pi)t`, squaring, and maximizing `t exp(-t²)` gives the
fully analytic upper bound

\[
\sup \text{certificate}
\le\sqrt{\frac{2}{\pi\sqrt{\pi e}}}
\approx0.4667442461530095<\frac12.
\]

Even the strict comparison needs no numerical optimization: `pi>3` and
`e>8/3` imply `pi³e>72>64`. The exact sharper scalar ceiling is the maximum
of the displayed one-variable function. Its unique positive maximizer
satisfies `phi(t)=t(2Phi(t)-1)`, because the difference has derivative
`-(2Phi(t)-1)-3t phi(t)<0`. Independent numerical diagnostics give
`t≈.657306552578`, `u≈.489016175949`, and maximum `.449554961564064`.
Only the symbolic supremum and the analytic `.466744...` expression are
being used as proved bounds here.

### Necessary scope restriction

This bounds the **paired direct-energy certificate**, not every actual
energy output of such a rounding procedure. The two unoriented quadratic
terms cancel in the half difference and have not been bounded by this
argument. For example, on a regular symmetric Hadamard family, `F=1`, `H=0`
has paired certificate zero while its two constant sign outputs can both
have normalized energy tending to `1/2`. The regular family already has
bounded normalized row sums, so it also satisfies the gauged hypotheses.

Thus the valid conclusion is that optimizing this finite-first-layer
direct-energy identity alone cannot reach `1/2`. Additional control of the
unoriented energies, a genuinely deeper response, or other structural
information is outside this ceiling theorem.
