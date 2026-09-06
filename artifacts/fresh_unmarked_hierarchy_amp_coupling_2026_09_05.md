# Unmarked transport coupled to the old hierarchy via a precise AMP theorem

Date: 2026-09-05. Status: algebra independently checked the cited primary
theorems, `Phi^2=0`, the final covariance row, the Schur complement, and the
query regularization, and passed the calculation. Exact for conference sequences; also valid under
the explicit all-power conditions below. Not a universal low-cap theorem.

Final scope update: this is an independently audited restricted AMP
route. Its all-power delocalization hypotheses remain mandatory. The
subsequent actual-sign weighted projection and normalized-gain proofs
bypass AMP for their specific identities and universal lower bounds;
they do not extend the entire formula (9) to every bounded-op signing.
The exact Section 6 improvement is over its explicitly named older mask,
not over the later best universal certificate. Its larger optimization
decimal is exploratory, as marked there. No convergence result for the
original minima follows from this note.

## 1. Imported results and exact matrix hypotheses

[Wang--Zhong--Fan, Proposition 2.7(b2) and Theorem 2.8](https://arxiv.org/pdf/2206.13037)
give finite-history AMP state evolution with independent side information
after uniform signed-permutation conjugation, provided the operator norm is
bounded, the spectral law converges, and, for every fixed integer `k>=1`
and every `epsilon>0`, eventually

\[
 \max_i |(B^k)_{ii}-n^{-1}\operatorname{Tr}B^k|
 +\max_{i\ne j}|(B^k)_{ij}|=O(n^{-1/2+\epsilon}).       \tag{1}
\]

The responses must be continuous with polynomial growth and Lipschitz in
the Gaussian-history arguments. Their auxiliary-variable law must have
all moments and polynomial density in `L^2`. The Gaussian covariance at
each queried step must be nonsingular. Independent signs (and extra
independent signs used for regularization) satisfy the auxiliary-law
requirements.

Every conference sequence satisfies (1) exactly: `B^(2k)=I` and
`B^(2k+1)=B`, with zero odd diagonal. The spectral law is
`(delta_-1+delta_1)/2` because `Tr B=0`; the operator norm is one. Signed
permutation conjugation preserves the original same-spin quadratic norm.

The state-evolution prescriptions are imported from
[Fan, Section 4.1, equations (4.4)--(4.7)](https://arxiv.org/pdf/2008.11892).
In the notation used below, if `Delta=E UU^T` and
`Phi_rs=E partial_s U_r`, they are

\[
 \text{Onsager}=\sum_{j\ge0}\kappa_{j+1}\Phi^j,
 \qquad
 \Sigma=\sum_{j\ge0}\kappa_{j+2}
       \sum_{a+b=j}\Phi^a\Delta(\Phi^T)^b,             \tag{2}
\]

where `kappa_j` are free cumulants of the limiting spectrum. The matrix
transpose convention in Fan's displayed Onsager matrix merely places the
coefficient of an earlier input in the appropriate later output; (2) is
written with the response Jacobian lower triangular.

## 2. Finite triangular old hierarchy and appended unmarked features

Let `X=(X_1,...,X_d)` be independent standard Gaussians and let `S` be an
independent sign. Choose a finite triangular family

\[
 a_1=1,\qquad a_j=a_j(X_1,...,X_{j-1}),\qquad
 E[a_j a_k]=\delta_{jk}.                              \tag{3}
\]

Initially take the functions smooth and globally Lipschitz, with the
required polynomial-growth extensions in the auxiliary sign. At matrix
level build fields sequentially from inputs

\[
 u_j=S\,a_j(x_1,...,x_{j-1}),\qquad x_j=Wu_j,
 \qquad W=\Pi B\Pi^T.                                \tag{4}
\]

All limiting mean Jacobians vanish because `E S=0` and the Gaussian
history is independent of `S`. Hence all Onsager corrections in (4) are
zero. Equations (2)--(3) give `Sigma=Delta=I_d`: these fields have exactly
the desired independent Gaussian empirical limit, jointly with `S`.

For the literal uniform Lipschitz hypothesis in the side variable's full
ambient domain, use a bounded continuous clipping extension of `S` (and
of the later auxiliary `T`) away from its support `{ -1,1 }`. The extension
equals the sign on that support and changes no actual iterate. This avoids
mistakenly using an unbounded Lipschitz constant proportional to a formal
real-valued auxiliary argument.

This includes smooth finite approximations of the previously constructed
marked-tree hierarchy: its input functions are the normalized even
Hermite products of earlier child fields. Smooth approximation and finite
Gram--Schmidt preserve triangularity; every approximation is fixed before
the matrix-order limit. The old finite certificate is therefore recovered
arbitrarily closely without assuming an unbounded-depth AMP theorem.

Choose further smooth globally Lipschitz features `h_1,...,h_r` of `X`
such that

\[
 E[h_\ell h_k]=\delta_{\ell k},\qquad
 E[X_j h_\ell(X)]=E[\partial_j h_\ell(X)]=0
 \quad\text{for all }j,\ell.                         \tag{5}
\]

There is no parity requirement. Append the matrix fields

\[
 z_\ell=W h_\ell(x).
\]

Their inputs also have zero mean Jacobians. Their correlations with the
old inputs vanish due to the independent factor `S`; their own input Gram
matrix is the identity by (5). Thus the joint empirical limit is

\[
 (S,x,z)\ \longrightarrow\ (S,X,Z),\qquad
 (X,Z)\sim N(0,I_{d+r}),\quad S\perp(X,Z).            \tag{6}
\]

This is a rigorous finite theorem for unmarked transport of arbitrary
old-hierarchy statistics, not a guess based on the first Gaussian field.

## 3. Exact final-energy calculation

Take smooth bounded responses `F(X,Z), H(X,Z)` and write

\[
 f=F(x,z)+S H(x,z),\qquad
 \alpha=E\nabla F(X,Z),\quad
 c_j=E[a_j(X)H(X,Z)],\quad
 c_{d+\ell}=E[h_\ell(X)F(X,Z)].                       \tag{7}
\]

Assume `|F|+|H|<=1` when using `f` for Boolean rounding. This constraint
is not needed just to compute the energy.

Append `f` as one final AMP input. Its mean Jacobian is `alpha^T`, since
the terms carrying `S` again average to zero. All earlier rows of `Phi`
are zero. Therefore `Phi^2=0`, and its prior input Gram block is `I`.
The sign matrix has limiting spectral mean zero and second moment one,
so `kappa_1=0`, `kappa_2=1`. Equation (2) then gives

\[
 Wf=y+\sum_{j=1}^{d+r}\alpha_j u_j,\qquad
 \operatorname{Cov}(y,(X,Z))=c+\kappa_3\alpha.        \tag{8}
\]

Here `y` has a joint Gaussian empirical limit with `(X,Z)`, independent
of the auxiliary sign. Gaussian integration by parts and (7)--(8) yield

\[
 \frac1n E[f^T y]\to\alpha\cdot(c+\kappa_3\alpha),
 \qquad
 \frac1n E\Big[f^T\sum_j\alpha_j u_j\Big]
       \to\alpha\cdot c.
\]

Consequently

\[
\boxed{
 \frac1{2n}E f^T Wf\to
 \sum_{j=1}^d E[\partial_{X_j}F]E[a_jH]
 +\sum_{\ell=1}^r E[\partial_{Z_\ell}F]E[h_\ell F]
 +\frac{\kappa_3}{2}\|E\nabla F\|_2^2.
}                                                       \tag{9}
\]

The first sum is the old paired certificate; the second is the new
unoriented energy. For conference matrices `kappa_3=0`. Under (1),
`kappa_3=lim Tr(B^3)/n`; applying the construction to either `B` or `-B`
allows the favorable sign of this spectral correction for a quadratic
absolute-norm lower bound.

The finite-dimensional covariance condition has not been ignored. Its
last Schur complement is

\[
 E(F^2+H^2)-\|c\|_2^2
       +(\kappa_4-\kappa_3^2)\|\alpha\|_2^2.         \tag{10}
\]

The old inputs, new feature inputs, and linear Gaussian coordinates form
an orthonormal collection in the scalar `L^2` model. Bessel's inequality
gives `E(F^2+H^2)>=||c||^2+||alpha||^2`. Also the scalar moment inequality
`m_4>=1+m_3^2` gives `kappa_4-kappa_3^2>=-1`. Thus (10) is nonnegative.
If it vanishes, replace only the final query input by `f+epsilon T`, where
`T` is an extra independent auxiliary sign. This makes (10) strictly
positive. The operator bound controls the energy difference by
`O(epsilon+epsilon^2)`; take `n -> infinity` and then `epsilon -> 0`.
This query regularization need not obey the rounding constraint because
the final lower bound is applied to the original cube-valued `f`.

## 4. Why this can improve an already good paired certificate

Suppose an old finite certificate uses a decision field `V=v dot X` and
`F_0=sign(V)(1-H_0(X))`, with `H_0` even. If the residual

\[
 r(X)=F_0(X)-\sum_j E[X_jF_0]X_j
\]

is nonzero, its normalized version is a permissible limiting feature `h`
after smooth Lipschitz approximation and exact linear projection. It has
`E[h F_0]=||r||_2>0`. Replacing the decision field by `V+delta Z` makes
`E[Z F_delta]` positive to first order in `delta`, whereas the loss in the
old pointwise-optimal decision objective is second order, provided the
unmasked conditional density on `V=0` is positive. Hence the new second
sum can strictly improve that fixed paired certificate. This paragraph
is a variational mechanism, not a certified numerical gain: the exact
density condition and approximation margins must be checked for a chosen
certificate.

## 5. Scope boundary

Theorem 2.8 is applicable here because conference matrices satisfy every
fixed-power delocalization condition, not merely entrywise flatness. A
general bounded-operator sign matrix need not satisfy (1); a low quadratic
cap gives still less. Thus (9) is presently a conference/full-power-
delocalized result. It is not a proof that the original universal lower
bound has increased, nor a proof of convergence of the original minima.

## 6. A fully quantitative strict gain over the existing two-field mask

Use exactly the compact-support mask `H=Htilde(V,W)` in Section 4 of
`fresh_limit_gaussian_mask_ascent_2026_09_05.md`; put `M=1-H`. Here
`Var V=1`, `Var W=p`, `Cov(V,W)=w`, and `s^2=p-w^2`. Let

\[
 F_0=M\operatorname{sign}(W),\quad
 E[F_0\mid\text{Gaussian first chaos}]=aV+bW,
 \quad r^2=E F_0^2-E(aV+bW)^2.
\]

Take the appended feature
`h=(F_0-aV-bW)/r`. Its use is justified by fixed smooth approximation,
orthogonal projection off the Gaussian linear coordinates, normalization,
and then taking the approximation limit in (9). Only Gaussian `L^2`
continuity is needed for the limiting energy formula after writing its
derivative expectations as first Gaussian moments.

Let `Z` be the resulting independent standard Gaussian. Keep `H` fixed
and replace the decision by

\[
 F_\tau=M\operatorname{sign}(W+\tau Z),\qquad \tau=1/125.
\]

The regression `E[UH\mid V,W]=cV+dW` gives the old paired value
`J_0=E[(cV+dW)F_0]`. Integrating only the new standard Gaussian,
with `e_tau(w)=2 Phi(w/tau)-1`, gives the exact coupled expression

\[
 J_\tau=E[(cV+dW)M e_\tau(W)] +2 D_\tau B_\tau,
\]
\[
 D_\tau=E[M\phi(W/\tau)],\qquad
 B_\tau=E[hM e_\tau(W)].                              \tag{11}
\]

The previously certified scalar intervals imply

\[
 0<c<.18,\quad0<d<.75,\quad |a|<.5,\quad|b|<.2,
 \quad r>.3.
\]

The new exact script checks these inequalities directly from that saved
certificate and independently verifies that the mask contains the
rectangle `|V|<=17/20, |W|<=1/6`. For the latter it encloses the degree-200
polynomial `g` on 34 rational Taylor bins and proves `g>13/20`; the
already certified coefficient `a_old>149/250` then suffices. Consequently

\[
 1_{\{|V|>1\}}\le M
 \le1_{\{|V|>17/20\}}+1_{\{|W|>1/6\}}.              \tag{12}
\]

Define `D_V=E[M |V| phi(W/tau)]`,
`D_VV=E[M V^2 phi(W/tau)]`, and
`D_W=E[M |W| phi(W/tau)]`. Symmetry of `M(v,w)` in `w`, the identity
`(p(w|v)-p(-w|v))/(p(w|v)+p(-w|v))=tanh(w v w_cov/s^2)`, and Mills'
inequality `x Phi(-x)<=phi(x)` prove

\[
 J_0-E[(cV+dW)M e_\tau(W)]
 \le 2\tau\{.18\,w D_{VV}/s^2+.75 D_\tau\}.         \tag{13}
\]

Also `2 Phi(-x)<=exp(-x^2/2)` and the residual formula give

\[
 B_\tau\ge .3-
 \frac{D_\tau+.5D_V+.2D_W}{.3\phi(0)}.                \tag{14}
\]

All quantities needed to bound (13)--(14) have elementary closed forms.
Multiplication of the `(V,W)` Gaussian density by `phi(W/tau)` gives total
mass `D_*=phi(0) tau/sqrt(p+tau^2)` and a centered tilted Gaussian with

\[
 \sigma_V^2=(s^2+\tau^2)/(p+\tau^2),\qquad
 \sigma_W^2=p\tau^2/(p+\tau^2).
\]

Use (12), exact one-dimensional Gaussian tail moments in `V`, and
Cauchy--Schwarz for the negligible `W` tail. Its probability is below
`10^-50`, verified by a finite positive rational Taylor lower bound for
`exp((1/(6 tau))^2/2)`. Thus no new numerical quadrature is required.

The exact Fraction computation in
`computations/fresh_unmarked_mask_tilt_certificate.py` proves

\[
\boxed{J_{1/125}-J_0
 >0.0000010923810129806398248858>10^{-6}.}               \tag{15}
\]

Combining with the prior certified lower endpoint gives the
conference/full-power-delocalized lower bound

\[
 J_{1/125}>0.4297875431186410219952.
\]

The numerical-only script `computations/fresh_unmarked_mask_tilt.py`
suggests `J_tau=0.4308115641...` near `tau=.25361`. That larger decimal
is not claimed as an interval certificate. The substantive banked result
here is the exact finite new interaction and a quantified positive gain;
other subsequently improved universal masks may already exceed this
particular numerical benchmark.
