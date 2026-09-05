# Independent audit: the second rooted response and smooth-mask certificate

Date: 2026-09-05. Auditor: fresh-limit literature subagent.

## Verdict

**Passed.** I independently checked the new proof in
`fresh_limit_custom_response_2026_09_05.md` and the exact-function extension in
`fresh_limit_response_variational_2026_09_05.md`. The cube moment induction,
endpoint errors, indirect-path cancellation, coordinate-pair factorization,
and scalar certificate all withstand the audit. No mathematical repair was
identified.

Together with the previously audited rooted Gaussian theorem, they imply

\[
\liminf_n M_n/n^{3/2}
\ge 0.385785876908778466066127787300791530025224664775532506435400
>0.38578.
\]

This is an asymptotic lower bound. Neither convergence nor an upper-bound
improvement follows. The decimal is a rigorous rational lower endpoint for a
fixed explicit smooth mask, not a certified global numerical optimum.

## 1. Pointwise gradient and high-moment estimate

Keep `m=n-1`, `B=A/sqrt(m)`, `Q=B²`, and a fixed cap `q(A)<=C n^(3/2)`.
The earlier elementary bootstrap gives `beta(B)=O_C(n)` and
`||B||op²=||Q||op=O_C(sqrt(n))`.

For a fixed pair `i,j`, put `a_k=Bik Bkj` and
`T=sum_k a_k S_k h'(Gk)`. Hollowness gives `ET=0`. Also
`sum a_k²<=1/m` and `sum |a_k|<=1`.

Under a flip of `S_l`, Taylor expansion of `h'` gives, for the half-difference,
the direct term `-a_l S_l h'(Gl)`, the indirect term
`-S_l sum_k Blk a_k S_k h''(Gk)`, and a remainder bounded by `C_h/m`.
The last estimate follows by summing `|a_k| Blk²`, not by summing an
unweighted error over `k`. Therefore, pointwise on the cube,

\[
\sum_l|\delta_lT|^2
\le C_h\left(\|a\|_2^2+
       \|B\|_{op}^2\|a\|_2^2+n/m^2\right)
\le C_h(1+\|B\|_{op}^2)/m=O_{h,C}(n^{-1/2}).
\]

The claimed moment induction is valid. To make it explicit, the cube
Poincare inequality is `Var f<=E sum_l(delta_l f)²`; it follows directly by
expanding `f` in the orthogonal Walsh basis, since a nonconstant coefficient
is counted once in the variance and at least once in the gradient sum.
If the pointwise sum is at most `w`, then `ET²<=w` and
`|T(S^l)-T(S)|<=2sqrt(w)`. The power-difference identity gives

\[
\operatorname{Var}(T^k)
\le k^2w E(|T|+2\sqrt w)^{2k-2}.
\]

Moreover `|ET^k|²<=ET² E|T|^(2k-2)`. Combining these two inequalities and
inducting yields `E|T|^(2k)<=C_k w^k`. Thus
`E|Tij|^8=O(n^-2)`, and for each fixed column

\[
E\max_i|T_{ij}|\le
 \left(\sum_i E|T_{ij}|^8\right)^{1/8}
=O(n^{-1/8})=o(1).
\]

This uses no Gaussian concentration transfer. The moment constants may
depend on the fixed smooth mask, which is sufficient for the ordered limits.

## 2. The endpoint increment and two-site Taylor error

Let `D_jYi` be the half difference obtained by fixing `Sj=+1` and `Sj=-1`.
The exact formula in the source correctly separates `Bij h(Gj)` from all
indirect paths. Its second moment is `O_h(1/m)`.

The stronger restoration estimate

\[
D_jY_i=B_{ij}h(G_j)+T_{ij}(S)+R_{ij},
\qquad \sup_{i,j}\|R_{ij}\|_2=O_h(1/m)
\]

is also valid. The central derivative differs from `h'(Gk)` by a term whose
principal contribution, up to sign, is

\[
\frac{S_j}{m}\sum_{k\ne j}B_{ik}S_k h''(G_k-Bkj S_j).
\]

The sum is bounded in L²: `h''` is even, so the prior even-function
transport lemma applies after restoring `Gk`, and the restoration itself is
pointwise bounded by `||h'''||infinity sum_k|Bik Bkj|=O_h(1)`.
The remaining third-order increments, including their row weights, are
`O_h(1/m)` pointwise. This is the important place where evenness of `h`
enters the sharper error estimate.

The same exact increment is pointwise bounded by
`||h||infinity/sqrt(m)+||h'||infinity sum_k|Bik Bkj|=O_h(1)`.
Consequently its third absolute moment is `O_h(1/m)`. Its expected maximum
over the output coordinate is `o(1)` by the preceding high-moment estimate
and the coarser pointwise restoration error `O(m^-1/2)`.

For clarity, write the endpoint-independent midpoints as `U_i,U_j` and the
increments as `D_i,D_j`. Averaging `Sj` in
`Sj F(U_i+Sj D_i) H(U_j+Sj D_j)` gives the first differential
`D_i F'(U_i)H(U_j)+D_j F(U_i)H'(U_j)` with remainder bounded by
`C(|D_i|+|D_j|)^3`. Evaluating the derivatives at the actual fields adds
only `C(D_i²+D_j²)`. Substituting the sharper increment formula adds
`O(1/m)` in expected absolute error. Thus the claimed total expected error
is `O(1/m)` per ordered pair before multiplication by `Bij`, and

\[
\sum_{i,j}|B_{ij}|O(1/m)=O(\sqrt n)=o(n).
\]

No unjustified independence between the actual fields and `Sj` is used.

## 3. Both indirect contributions really vanish

For the off-diagonal indirect contribution, reordering gives

\[
D_1=\sum_k E[S_k h'(G_k)C_k],\quad
C_k=m^{-1}p_k^\top Bq_k,
\]

with `(p_k)i=sqrt(m)Bki F'(Yi)` and
`(q_k)i=sqrt(m)Bki H(Yi)`. All coordinates are bounded by fixed constants.
The exact identity for the difference of two bilinear products gives

\[
|D_kC_k|\le C_{F,H}\frac{\beta(B)}m
                       \max_i|D_kY_i|.
\]

Here the bilinear cap applies to arbitrary vectors in the corresponding
bounded cubes by multilinearity. It applies separately to the two endpoint
differences, even though every vector depends on the same input signs.
Since `h'(Gk)` is exactly independent of `Sk`, endpoint averaging gives
`E Sk h'(Gk) Ck=E h'(Gk) DkCk`. As `beta(B)/m=O(1)` and the expected
maximal increment is `o(1)`, this proves `D1=o(n)`.

This is a collective bound on all indirect paths. Bounding individual
`Tij` in L² and summing their magnitudes would not suffice, but that is not
the argument being used.

For the own-coordinate term, flatness gives exactly

\[
T_{jj}=Z_0-m^{-1}S_jh'(G_j),\quad
Z_0=m^{-1}\sum_kS_kh'(G_k).
\]

Endpoint extraction gives `EZ0²=O(1/m)`. The common-factor contribution is
bounded by `C beta(B) E|Z0|=O(sqrt(n))`. The removed contribution is a
single bounded-cube bilinear form divided by `m`, hence is
`O(beta(B)/m)=O(1)`. The asserted `D2=o(n)` follows with the stated
normalization and signs.

## 4. Coordinate-pair Gaussian factorization and final scalar formula

For any fixed even rooted Hermite order `r`, the exact Gaussian summand
covariance yields

\[
\operatorname{Cov}(F_{i,r},F_{j,r})
=r!Q_{ij}+O_r(\|Q\|_{op}/m).
\]

The error is uniform because both sandwiching rows have Euclidean norm one.
Different chaos degrees are orthogonal, and their covariances with the
linear fields are zero. Thus along any coordinate-pair sequence with
`Qij->0`, the multivariate fixed-chaos theorem applies to both coordinates
at once, including repeated chaos orders, and gives independent coordinate
blocks. The existing joint Rademacher replacement and L² tail transport
remain valid for this finite enlarged vector. Therefore
`(Yi,Gj,Yj)` has the claimed limit, with `Yi` independent of `(Gj,Yj)`.

This is enough for averaging. For example, the fraction of pairs with
`|Qij|>n^-1/8` is at most `Tr(Q²)/n^(7/4)=O(n^-1/4)`.
If the average bounded-test error on the other pairs failed to vanish,
one could choose a violating coordinate-pair sequence there; its Gram entry
tends to zero, contradicting the just-proved joint convergence. Thus no
unproved uniform mixing estimate is needed.

Combining the direct term with the two vanishing terms proves exactly

\[
\frac1n E F(Y)^\top B[S H(Y)]
\longrightarrow E F'(Y_*)\,E[h(G_*)H(Y_*)],
\]

where `Var G*=1`, `Var Y*=Eh²`, and `Cov(G*,Y*)=Eh`.
For the smoothed outer threshold, conditional independent rounding realizes
the two means as genuine Boolean vectors. The half difference of their
energies is the numerator above, with no additional factor of two, and is
at most `q(B)`.

For fixed masks and fixed positive smoothing parameters all estimates hold
uniformly along a fixed-cap sequence. Dimension tends to infinity first.
For a hard outer threshold, the Gaussian identity
`EF_epsilon'(Y*)=2phi_(v+epsilon²)(alpha)` permits the subsequent smoothing
limit directly. Bounded convergence gives the accompanying indicator
expectation. For a hard first mask one may subsequently let its smoothing
vanish too. This proves the stated lower bounds; it does **not**, without
another argument, assert convergence of the original energy of a fixed
discontinuous hard-rule algorithm.

## 5. Explicit smooth-mask extension and exact certificate

At `Eh²=1` and `Eh=rho`, let
`Y=rho G+sqrt(1-rho²)Z` and
`k(g)=P(|Y|<=alpha | G=g)`. Then `Ek=p=2Phi(alpha)-1` and
`s²=Var k>0`. The proposed function

\[
h_*(g)=\rho+\frac{\sqrt{1-\rho^2}}s(k(g)-p)
\]

is genuinely bounded, smooth, and even, with bounded derivatives of every
fixed order. Directly `Eh*=rho`, `Eh*²=1`, and
`Eh*k=rho p+sqrt(1-rho²)s`. Thus the Hilbert-space Cauchy--Schwarz
optimizer is achieved inside the actual admissible class. This is not a
relaxation requiring an unproved realization theorem.

The intermediate vector `S h*(BS)` need not be Boolean or have magnitude
at most one. Only the final smoothed-response means need lie in the cube,
and they do. Alternatively multiplying `h*` and the outer threshold by
the same positive constant leaves the scalar energy formula unchanged.
There is no circular normalization: `k` is first defined at the chosen
fixed `rho,alpha`, and the displayed formula then verifies the exact
moments needed for its own Gaussian limit.

I independently checked the positive Hermite series

\[
s^2=4\phi(\alpha)^2\sum_{j\ge1}
 \frac{\rho^{4j}H_{2j-1}(\alpha)^2}{(2j)!}.
\]

The indicator coefficient is obtained by integrating
`(phi H_(2j-1))'=-phi H_(2j)`, conditional Gaussian expectation multiplies
degree `2j` by `rho^(2j)`, and Parseval squares that factor. For a cutoff
`K`, the omitted variance is bounded above by
`rho^(4(K+1)) p(1-p)`, because all terms are nonnegative and the full
unweighted indicator variance is `p(1-p)`.

I inspected and executed
`computations/fresh_limit_response_variational_certificate.py`, including
its previously audited outward-interval dependency. The Hermite recurrence
and factorial indices are correct. At `rho=47/50`, `alpha=81/100`, `K=100`,
the exact rational run reproduces

\[
0.385785876908778466066127787300791530025224664775532506435400
\le D\le
0.385785876909691210517840971183287938897036818354614542045282.
\]

The lower bound uses no numerical optimization, quadrature, or numerical
differential equation. The optimizer's exact defining integral need not be
evaluated to construct the mathematical witness; the positive series and
tail bound suffice to certify its achieved value.

## 6. Independent numerical stress checks (not proof)

I also tested a fixed smoothed first indicator and fixed smoothed outer
threshold, with `t=1.3`, `alpha=.76`, both smoothing parameters `.1`, and
random seed `99291`. Independent one-dimensional Gaussian quadrature gives
the predicted value `.3661345415705465`. With 4096 random inputs at order
4096, a hollow Sylvester matrix gave `.3661848274`, standard error
`.0000292597`. Attaching 64 all-positive twins, a genuine low-cap family
with `||Q||op` of order `sqrt(n)`, gave `.3660396212`, standard error
`.0001393177`. Smaller orders also trend to the same formula.

These checks do not certify the theorem or rule out all counterexamples.
Their role was to test a structurally difficult allowed family in addition
to the independent proof audit. The variational auditor reported a negative
finite hard-rule example at order six; that does not contradict the ordered
smooth asymptotic theorem or the liminf certificate.
