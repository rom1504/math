# Custom second response: independent falsification and proof audit

Date: 2026-09-05. Auditor: fresh-limit variational agent.

## 1. Exact finite counterexample to unsmoothed positivity

Use `t=7/8`, `a=2phi(t)`, `B=A/sqrt(n−1)`,
`v=S 1{|BS|≤t}`, `Y=Bv`, and `y±=sign(aS±Y)`.
The oriented response energy is

\[
e_y=\frac1{4n}E[(y^+)^TBy^+-(y^-)^TBy^-].
\]

The matrix of order six whose off-diagonal entries are all +1 except the
three negative edges inside vertices `{1,2,3}` (zero-based indexing) gives

\[
e_y=-\frac3{32\sqrt5}<0.
\]

This was found by exhaustive integer enumeration of every signing modulo
diagonal switching and every spin modulo global reversal, and independently
checked in the explicit seed. The raw energy-difference sum over the 32
projective spins is −72. The first mask is evaluated by the exact integer
comparison `64 Graw²≤49(n−1)`. At order six, `1<a sqrt5<2`, so the second
mask is exactly `|Yraw|≤1`. No floating-point optimizer is used.

`computations/fresh_custom_response_small_exact.py` covers all orders 3–7.
The minimum raw sums are respectively `0,0,0,−72,920`. Results, counts, and
matrices are saved in
`computations/results/fresh_custom_response_small_exact_2026_09_05.json`.
This disproves a finite-dimensional positivity identity, **not** asymptotic
positivity for low-cap sequences.

## 2. Structured asymptotic counterexample search

`computations/fresh_custom_response_falsify.py` tested the required response
and the separate cross term. The 38-case saved run uses orders 128 and 256,
8192 independent whole spin vectors per matrix, and includes random signs,
signs of Gram matrices with rank ratios 0.03 through 4, added Gaussian noise,
positive and negative planted cliques of size `n^(3/4)`, regular Hadamards,
and full-seed Hadamard tensors. Standard errors use independent whole-spin
samples, not correlated individual coordinates.

Every saved estimate of `e_y` is positive: approximately 0.307–0.331 at
order 128 and 0.320–0.323 at order 256. Spectral third moments range to
about 3.3 in the Gram tests. This asymmetry changes the *unoriented* energy
substantially but not the sign of the oriented difference. Data are in
`computations/results/fresh_custom_response_falsify_2026_09_05.json`.

The negative six-point seed was separately lifted with `H4` tensors and
three diagonal completions. Its responses became positive, roughly 0.23–0.29
at order 24, 0.30 at order 96, and 0.316–0.320 at order 384. No scalable
negative example was found. These simulations do not certify a limiting law.

## 3. Audit of the smoothed positive-energy formula

The parent/algebra proof in `fresh_limit_custom_response_2026_09_05.md`
uses fixed smooth even `h`, fixed bounded smooth odd `F` and even `H`, and
claims

\[
\frac1n E F(Y)^TB[S H(Y)]
\longrightarrow E F'(Y_*)\,E[h(G_*)H(Y_*)],
\]

where `Y=B[S h(BS)]` and the Gaussian pair has
`Var G*=1`, `Var Y*=E h(Z)²`, `Cov(G*,Y*)=E h(Z)`.
I independently audited its new steps rather than importing a second-AMP
universality assumption.

### Maximal influences

For `Tij=sum_k Bik Bkj Sk h'(Gk)`, the coefficient vector has squared norm
at most `1/m` and absolute sum at most one. Endpoint-spin extraction gives
`E Tij²=O(1/m)` and the mean is exactly zero.
The pointwise cube gradient decomposes into a diagonal term, a term of the
form `B[coefficients·S·h''(G)]`, and remainders bounded by `C/m` in each
coordinate. Therefore its squared gradient sum is at most
`C(1+||B||op²)/m=O(n^(-1/2))`.

The Efron--Stein moment iteration is valid. If a centered cube function has
gradient square sum at most `w`, apply the variance inequality to `T^k` and
use `|ET^k|²≤ET² E|T|^(2k−2)`. Induction gives `E|T|^(2k)≤C_k w^k`.
Thus `||Tij||8=O(n^(-1/4))`, and
`E max_i |Tij|≤(sum_i E|Tij|8)^(1/8)=O(n^(-1/8))`.
The exact increment of `Yi` under one input spin differs from its direct
term plus `Tij` by a pointwise `O(m^(-1/2))`, so the same maximal-influence
conclusion holds for `D_j Yi`.

For the finer two-site expansion, restoration of the endpoint in `h'`
has leading error `(Sj/m)sum_k Bik Sk h''(Gk^(j))`. Evenness of `h''`
allows the previously proved transport bound, giving `L²` error `O(1/m)`.
All remaining Taylor errors have the claimed order. The exact increments
also have second moments `O(1/m)`.

### Two-site and indirect-path cancellation

Central Taylor expansion in `Sj` produces the direct term
`m^(-1)sum_ij E F'(Yi)H(Yj)h(Gj)`, plus the two indirect terms in the
candidate. Per-pair error `O(1/m)` is sufficient, since multiplication by
`sum_ij |Bij|=n sqrt(m)` makes its total `O(sqrt(n))=o(n)`.

For the off-diagonal indirect term, the exact reordered coefficient is
`Ck=m^(-1) pk^T B qk`, with uniformly bounded vectors depending on `Y`.
The two-endpoint difference of this bilinear expression gives

\[
|D_k C_k|\le C\frac{\beta(B)}m\max_i|D_kY_i|.
\]

Since `h'(Gk)` is exactly independent of `Sk`, discrete integration by parts
and `beta(B)=O(n)` make the sum over `k` equal to `o(n)`. This controls
the entire indirect term together; pointwise-small entries alone would not
have justified dropping it.

For the own-coordinate indirect term, flatness gives
`Tjj=Z0−Sj h'(Gj)/m`, with `Z0=m^(-1)sum_k Sk h'(Gk)` and
`E Z0²=O(1/m)`. The bilinear cap bounds its accompanying sum by `O(n)`
pointwise, giving total `O(sqrt(n))`. The removed diagonal term is bounded
by `C beta(B)/m=O(1)`. Both controls check out.

### Direct-term factorization

At fixed even rooted degree, output covariance between roots `i,j` equals
`r! Qij+O_r(||Q||op/m)`, by sandwiching the exact rooted covariance kernel.
The already audited finite-chaos multivariate theorem therefore gives
independent coordinate blocks whenever `Qij→0`. This includes the mixed
triple `(Yi,Gj,Yj)` needed here. Since `Tr Q²/n²→0`, only a vanishing
fraction of ordered pairs have a correlation bounded away from zero.
Boundedness of `F',H,h` then proves the required average factorization.

## 4. Conclusion and scope

No fatal gap was found in the new smoothed positive-energy proof. In the
intended nonnegative-threshold channels its Gaussian limit is positive.
The dimension limit must be taken with smoothing fixed, and smoothing
removed only afterwards. This is sufficient for the claimed new lower
bound on `q(A)`; it does **not** automatically establish convergence of
the exact unsmoothed fixed-rule energy without an additional uniform energy
continuity argument. The finite order-six counterexample is therefore
compatible with the audited asymptotic smoothed theorem.

No original-limit convergence claim is made.

## 5. Independent audit of the optimized smooth first mask

The later variational calculation in
`fresh_limit_response_variational_2026_09_05.md` also passes an independent
audit. Normalize `E h(G)^2=1` and `E h(G)=rho`, with `0<rho<1`.
The limiting Gaussian pair is then
`Y=rho G+sqrt(1-rho²) Z`, independently of the remaining shape of `h`.
For the outer central threshold `H(y)=1_{|y|≤alpha}`, put

\[
k(g)=\Pr\{|\rho g+\sqrt{1-\rho^2}Z|\le\alpha\},\quad
p=Ek(G),\quad s^2=E(k(G)-p)^2.
\]

Both constraints on `h` determine the norm of its centered component:
`||h-rho||₂=sqrt(1-rho²)`. Cauchy--Schwarz consequently gives

\[
E[h(G)k(G)]\le \rho p+\sqrt{1-\rho^2}\,s,
\]

with equality for

\[
h_*(g)=\rho+\frac{\sqrt{1-\rho^2}}s(k(g)-p).
\]

This is self-consistent: its mean is exactly `rho`, its squared norm is
exactly one, and these are precisely the two quantities determining the
Gaussian law used to define `k`. Furthermore, `h*` is bounded smooth even
with bounded derivatives of every fixed order. It belongs to the proved
response class, not just a larger Hilbert-space relaxation. The intermediate
mask need not lie in `[-1,1]`; positive rescaling of both it and the outer
threshold would impose such a bound without changing the energy formula.
The optimization is exact at fixed `(rho,alpha)` and is not a claim that the
central threshold is the globally optimal outer response.

### Exact scalar certificate

The independently checked Hermite identity is

\[
s^2=4\phi(\alpha)^2\sum_{j\ge1}
 \frac{\rho^{4j}H_{2j-1}(\alpha)^2}{(2j)!}.
\]

The factor `rho^(4j)` is correct: conditional Gaussian expectation first
multiplies degree `2j` by `rho^(2j)`, and Parseval squares the coefficient.
All terms are nonnegative. After truncation at `K`, the remaining sum is
at most `rho^(4(K+1)) p(1-p)`, using the entire unweighted indicator variance
as an upper bound. This avoids any quadrature-error assumption.

I read and reran
`computations/fresh_limit_response_variational_certificate.py`, including
its imported exact interval primitives, previously audited separately.
It uses rational `rho=47/50`, `alpha=81/100`, the exact probabilists'
Hermite recurrence, and `K=100`. The returned interval is

\[
\begin{split}
0.385785876908778466066127787300791530025224664775532506435400
\le D(47/50,81/100)\\
\le
0.385785876909691210517840971183287938897036818354614542045282.
\end{split}
\]

In particular, `D>19289/50000=0.38578`. The scalar certification is exact;
its implication for the original problem uses the independently audited
smooth two-site theorem and then the ordered smoothing limit.

### Targeted tests of the actual optimized function

`computations/fresh_hstar_response_falsify.py` evaluates the explicit `h*`
at the above rational parameters, then uses the finite hard responses
`y+=sign(alpha S+Y)`, `y-=sign(alpha S-Y)`. It records oriented energy
`E(y+^T B y+-y-^T B y-)/(4n)`, unoriented energy, mask moments, and field
moments. Results are saved in
`computations/results/fresh_hstar_response_falsify_2026_09_05.json`.

The completed run comprises 26 models at orders 128, 256, and 512, with
8192 whole-spin samples each. Models include random signs, Gram-sign
matrices with rank ratios 0.1, 0.25, 1, and 4, positive and negative planted
cliques of size approximately `n^(3/4)`, and symmetric Hadamard/tensor
examples. There is no sign of a counterexample. Oriented energies range
from about 0.38570 to 0.38810; the order-512 range is 0.38588 to 0.38654.
Individual whole-spin Monte Carlo standard errors range from about
0.000086 to 0.00065. Transported-field variances range from 0.9986 to
1.0067, consistent with the theoretical normalization one. Asymmetric
Gram spectra and signed planted structures do not produce the suspected
negative oriented-energy defect.

These finite hard-rule experiments are corroborating tests, not numerical
proofs of convergence. In particular they do not remove the ordered-limit
caveat in Section 4.

## 6. Two-field extension and exact centering reduction

The multivariate endpoint formula in
`fresh_limit_two_field_variational_2026_09_05.md` passes an independent
audit. For fixed admissible smooth functions, it is

\[
\frac1n E F(G,Y)^T B[S H(G,Y)]
\longrightarrow E F_g(G_*,Y_*)\,E H(G_*,Y_*)
 +E F_y(G_*,Y_*)\,E[h(G_*)H(G_*,Y_*)].
\]

Indeed `D_j G_i=B_ij` exactly, while `D_j G_j=0`; the new direct term is
therefore exactly the asserted `F_g` term. All central Taylor cross errors
are controlled by `B_ij²=1/m` and the existing second moment of `D_jY_i`.
In the indirect bilinear derivative estimate, the only new contribution is
`max_i|D_kG_i|≤m^(-1/2)`. The own-coordinate indirect term still involves
only `H_y`. The joint rooted-chaos limit for four fields gives the two direct
factorizations outside a vanishing fraction of correlated coordinate pairs.
No new depth-universality assertion is needed.

For jointly even `0≤H(g,z)≤1`, its conditional average `k(g)=E_Z H(g,Z)`
is even. Taking `h=(k-p)/c`, where `p=EH` and `c²=Var k`, makes the two
limiting fields independent standard normals and gives `E hH=c`.
Gaussian integration by parts and pointwise optimization of `F` yield

\[
J(H)=E\bigl[|pG+cZ|(1-H(G,Z))\bigr].
\]

The measurable-mask extension also checks: Gaussian OU smoothing preserves
the interval `[0,1]`, joint evenness and mean, gives bounded derivatives at
each fixed smoothing parameter, and converges in Gaussian `L²`. Conditional
expectation contracts `L²`, so the conditional variance converges as well.
This justifies taking dimension first, sign smoothing second, and mask
smoothing last. If `c=0`, every OU-smoothed conditional average is constant,
and an arbitrary centered normalized smooth even first mask yields the
zero second coefficient required by the formula.

### Centering loses no generality for nonnegative outer masks

There is an additional reduction not needed for validity of the lower
certificate. It shows that the preceding centered construction already
captures the supremum over the whole smooth two-field response class with
nonnegative outer masks.

For a general first mask write `b=Eh`, `d²=Var h>0`, and
`Y=bG+dZ`. Pull the outer functions back to independent coordinates:
`K(g,z)=H(g,bg+dz)`, and similarly write the pulled-back odd channel as `F`.
Put `p=EK`, `k(g)=E_ZK(g,Z)`, and `s²=Var k`. The endpoint formula, followed
by Gaussian integration by parts, becomes

\[
E[(pG+rZ)F(G,Z)],\qquad
r=\frac{E[h(G)K(G,Z)]-bp}{d}
 =E\left[\frac{h(G)-b}{d}(k(G)-p)\right].
\]

Cauchy--Schwarz gives `|r|≤s`. Pointwise optimization over
`|F|≤1-K` gives

\[
J_K(r)=E|pG+rZ|(1-K).
\]

This is a convex function of `r`, so
`J_K(r)≤max(J_K(s),J_K(-s))`. The negative endpoint equals the positive
endpoint for the reflected admissible mask `K(g,-z)`. Each endpoint is
attained by an appropriately signed centered normalized first mask.
Consequently arbitrary first-mask means and variances do not enlarge the
supremum over `J(H)`. The degenerate case `d=0` has only a first-coordinate
response, corresponding to `r=0`, and obeys the same endpoint upper bound.

This statement is restricted to nonnegative outer masks. The full cube
constraint allows signed `H` with `|F|+|H|≤1`; the same reasoning produces
the potentially larger envelope

\[
J_{\rm signed}(H)
 =E|pG+cZ|(1-|H|),\quad
p=EH,\quad c^2=\operatorname{Var}(E_ZH).
\]

No reduction of this signed-mask envelope to nonnegative masks is asserted.
