# A positive second-return innovation with the actual remaining-slack mask

Date: 2026-09-06. Proved at the scalar zero-first scope below; independently
submitted for audit. This extends the constant-mask construction in
`continued_audit_actual_sine_innovation_floor_2026_09_06.md`. It is not a
floor for arbitrary marked responses, a larger universal lower bound, or
a convergence theorem.

## 1. Statement with all parameters fixed before matrix order

Let `B=A/sqrt(n-1)` be an actual symmetric hollow signing with
`||B||op<=L`, for fixed `L>=1`. Let f be a fixed bounded odd,
Gaussian-a.e.-continuous scalar function satisfying

```math
 |f|\le1,\qquad E[Nf(N)]=0,\qquad \tau^2=E f(N)^2>0.
```

Let H be fixed, bounded, even, and Gaussian-a.e.-continuous, with

```math
 0\le H(g)\le1-|f(g)|,\qquad
 \mu=E H(N)>0,\qquad h=\|H\|_\infty.
```

In particular `H=1-|f|` is allowed whenever its mean is positive. There
is NO assumption of strictly positive pointwise slack or a pointwise
variance floor. Choose fixed `w>0` and set `t=sqrt(w)/(L tau)`.
For independent Boolean seeds S define the ACTUAL fields

```math
 G=BS,\qquad Z=Bf(G),\qquad C=H(G)\sin(tZ),\qquad
 u=f(G)+C,\qquad J=1-|u|.
```

Then

```math
 \boxed{
 \liminf_{n\to\infty}\frac1n E\sum_i J_i |(Bu)_i|
 \ge \sqrt{\frac2\pi}\left(1-\frac2\pi\right)
 \frac{\mu^3}{h}
 \frac{w^3 e^{-w}}{6L^9\sqrt{\sinh w}}>0.}
 \tag{1}
```

The right side is deliberately conservative. The content is that an
entire feasible scalar-source/mask class has a genuinely positive next
innovation on arbitrary bounded-operator hollow signings. No involution
or spectral-flatness identity is used.

## 2. Actual coherent coefficient and the sufficient covariance statement

Put `Q=B^2` and, in normalized Hermite notation,

```math
 R=\sum_{p\ge3,\ p\text{ odd}} f_p^2 Q^{\circ p},\qquad T=BRB.
```

Every row of G has the same distribution, because all its nonzero
coefficients have magnitude `1/sqrt(n-1)`. Thus
`mu_n=E H(G_i)` is independent of i and tends to mu. The exact first
history coefficient from the marked-energy theorem is

```math
 e_i=e^{-t^2T_{ii}/2},\qquad a_i=\mu_n t e_i,\qquad
 L_0=B D_a Z,\qquad \eta=B(C-D_aZ).
```

The conditional mean is exactly zero and the old coherent source
coefficients vanish. In particular `Bu=Z+L_0+eta` is an identity.

Let U and N_T be independent Gaussian vectors of covariance Q and T.
Set `K_H(i,j)=E[H(U_i)H(U_j)]`. The normalized-nuclear covariance theorem
for a zero-first scalar source gives the ideal residual covariance

```math
 K_{\rm res}
 =K_H\circ[D_e\sinh(t^2T)D_e]
                -\mu^2t^2D_e T D_e,\qquad
 \frac1n\|\operatorname{Cov}(\eta)-BK_{\rm res}B\|_*\to0.
 \tag{2}
```

It also gives the SUFFICIENT symmetric cross statement

```math
 \frac1n\sum_i|E[\eta_i(L_0)_i]|\longrightarrow0.       \tag{3}
```

Here is the precise extension needed to justify (2)--(3). In the proof
of `continued_feedback_zero_first_masked_covariance_2026_09_06.md`, permit
a finite sum of products of even old responses and odd new responses.
Expand each covariance bilinearly before applying its forest estimates.
On both sides every old branch still has even total count and every
new branch has odd total count, with individual new source degree at
least three. The two-small-contraction and lone-star estimates therefore
apply term by term, with their original bounded coefficients. Whole
branch matches give precisely the covariance of the same finite sum
evaluated on independent U and N_T. No new branch pattern is introduced.

Bounded deterministic row-dependent coefficients are already allowed by
that proof. First approximate a_i by a finite coefficient catalog. To
include the linear response, use clipped linear functions, take the n
limit first, and remove the clipping afterwards. Fixed-polynomial
approximations of Z have bounded row moments and arbitrarily small
limiting averaged L2 error, so Z has uniformly integrable averaged
squares. The ideal Gaussian variances are bounded by `L^2 tau^2`.
Covariance Cauchy--Schwarz in nuclear norm controls every clipping and
coefficient error. These facts justify the finite-sum extension for
`C`, `D_a Z`, and `C-D_aZ`.

For the ideal variables,

```math
 E[H(U_i)\sin(t(N_T)_i)(N_T)_j]=\mu t e_iT_{ij}.
```

Thus the ideal residual has zero cross with its mean linear channel.
Replacing mu_n by mu has vanishing normalized-nuclear cost. Covariance
polarization of the row-weighted sums gives the SYMMETRIC cross matrix;
after multiplication on both sides by B this suffices for every diagonal
entry in (3). One must not infer an unrestricted nonsymmetric cross
matrix merely by scalar polarization. Diagonal nuclear-norm duality gives
the average ABSOLUTE error in (3), not only its signed average.

## 3. Positive cubic mass and a bounded ideal row variance

The Gaussian covariance of the even mask has the PSD decomposition

```math
 K_H=\mu^2\boldsymbol1\boldsymbol1^{\mathsf T}
                 +\sum_{q\ge2,\ q\text{ even}} H_q^2 Q^{\circ q}.
```

Consequently (2) is positive semidefinite, and its constant-mask cubic
piece remains positive even in the presence of all the other mask terms:

```math
 K_{\rm res}\succeq\frac{\mu^2t^6}{6}D_eT^{\circ3}D_e. \tag{4}
```

Indeed subtract only the constant-mask linear term from the entrywise
sinh series; all remaining Schur products have positive coefficients.
This positivity is particular to the present response and must not be
asserted for arbitrary row-dependent Hermite coefficients.

The correlation Schur multiplier estimate and row norm one give

```math
 R\preceq L^2\tau^2 I,\qquad
 T\preceq L^2\tau^2 Q,\qquad T_{ii}\le L^2\tau^2,
 \qquad \frac1n\operatorname{Tr}T\ge\tau^2.              \tag{5}
```

For the last inequality, expand `Tr(QR)` as
`sum_p f_p^2 sum_(i,j) Q_ij^(p+1)`: p+1 is even and Q_ii=1.
Using (4)--(5), then retaining the diagonal fourth powers, proves

```math
 \frac1n\operatorname{Tr}(QK_{\rm res})
 \ge\frac{\mu^2t^6}{6nL^2\tau^2}
                   \sum_{i,j}e_i e_j T_{ij}^4
 \ge\frac{\mu^2t^6\tau^6 e^{-w}}{6L^2}
 =\frac{\mu^2w^3e^{-w}}{6L^8}=:a_*.                    \tag{6}
```

The Schur multiplier of K_H has norm at most `E H(N)^2<=h^2`.
The bounds `||T||<=L^4 tau^2` and `max Tii<=L^2 tau^2`, applied termwise
to the positive sinh series, yield

```math
 0\preceq K_{\rm res}\preceq K_C,
 \qquad \|K_C\|_{op}\le h^2L^2\sinh(w)=:M_*.
```

Since B has unit row norms, each ideal variance
`v_i^*=(BK_resB)_ii` is at most M_*. If `v_i=E eta_i^2`, (2) gives
`average |v_i-v_i^*| ->0`, hence also convergence to zero of the average
square-root difference. It follows that

```math
 \liminf_n\frac1n\sum_i\sqrt{v_i}
 \ge a_*/\sqrt{M_*}.                                  \tag{7}
```

This uses an ideal row cap as well as an averaged variance lower bound.
An averaged variance lower bound alone would not justify (7).

## 4. The actual mask keeps enough mean slack

Feasibility gives the pointwise, actual inequality

```math
 J_i\ge H(G_i)[1-|\sin(tZ_i)|].
```

The established local comparison makes G_i and its new scalar Gaussian
noise independent, with their literal row variance T_ii retained. Thus
there are errors with average absolute value tending to zero such that

```math
 E J_i\ge\mu[1-E|\sin(t\sqrt{T_{ii}}N)|]-\epsilon_i
 \ge\mu\left(1-2/\pi\right)-\epsilon_i.                \tag{8}
```

For the last inequality use the absolutely convergent Fourier series
`|sin x|=2/pi-(4/pi) sum_(k>=1) cos(2kx)/(4k^2-1)`.
Gaussian averaging replaces each cosine by `exp(-2k^2 v)`, a nonnegative
number, so `E|sin(sqrt(v)N)|<=2/pi` for every v>=0. This proves (8)
for every fixed w>0 without a small-response restriction. Alternatively
`|sin z|<=|z|` gives the sharper small-w factor `1-sqrt(2w/pi)` whenever
useful. Neither argument is a hard-threshold approximation. If local
approximation is only available in averaged form, replace sqrt(v_i) by
sqrt(v_i*) before multiplying (8). The latter is uniformly bounded by
sqrt(M_*); its average absolute difference from sqrt(v_i) vanishes.
Thus no uncontrolled product of two merely averaged errors is used.
Combining (6)--(8),

```math
 \liminf_n\frac1n\sum_i E J_i\sqrt{v_i}
 \ge\mu(1-2/\pi)a_*/\sqrt{M_*}.                        \tag{9}
```

## 5. Ordered second-query comparison and extremal normalization

Apply Section 7 of
`continued_audit_ordered_bounded_next_gain_2026_09_06.md`. On each fixed
variance set `v_i>=epsilon`, (3) removes the regression correction
`rho_i/sqrt(v_i)` in average absolute value. The omitted average square
root is at most sqrt(epsilon). Take n first and then epsilon down to
zero. The scalar next-gain theorem and (9) give exactly (1).

This final step is not an assertion of a cutoff-free state evolution.
Both the retained coherent return and the current feasible u are actual
Boolean-seeded fields. The existing degree approximation and variance
cutoff orders remain in force.

For every seed realization the two cube points `u+J sign(Bu)` and
`u-J sign(Bu)` satisfy
`e(u+J sign(Bu))-e(u-J sign(Bu))=2 sum_i J_i |(Bu)_i|/n`,
when `e(x)=x^TBx/(2n)`.
Their maximum absolute half-energy therefore dominates the gain in (1).
Hollowness permits exact independent Boolean rounding. Thus this is a
certificate for `Q(A)/(n sqrt(n-1))`, with its exact normalization.

It is not a new lower bound beyond .4333221116640807. Its new content is
that second-return innovation cannot disappear for this entire genuine
scalar-source/feasible-mask class. The analogous assertion for arbitrary
high-value marked sources f(G,Y), or a mechanism forcing convergence of
the optimal normalized caps, remains open.

## 6. Positive next gain is an actual energy-improvement step

The following finite algebraic statement applies to ANY feasible random
u, not only the source class above. Write

```math
 \Lambda=\max_{x\in\{\pm1\}^n}|e(x)|,\qquad
 v=(1-|u|)\operatorname{sign}(Bu),\qquad
 j=u^{\mathsf T}Bv/n\ge0.
```

The cube points u+v and u-v give `j<=Lambda<=L/2`. Set `a=j/L<=1/2`
and `epsilon=sign(e(u))`, with epsilon=1 at zero. Then u+epsilon a v
is feasible, and exactly

```math
 \begin{aligned}
 |e(u+\varepsilon a v)|
 &\ge |e(u)|+a j-a^2\Lambda\\
 &\ge |e(u)|+j^2/(2L).
 \end{aligned}
```

Taking expectations and applying Jensen proves the operational inequality

```math
 \boxed{\ \Lambda\ge E|e(u)|+
       \frac{1}{2L}\left(\frac1n E\sum_i(1-|u_i|)|(Bu)_i|\right)^2.\ }
 \tag{10}
```

In particular a verified asymptotic next-gain floor kappa implies
`liminf_n [Lambda_n-E|e(u_n)|]>=kappa^2/(2L)`. Thus (1), and the
monotone variant below, give a strictly positive actual improvement over
the current fractional energy, not merely a separate small lower bound
on Lambda. The two-moment next-gain theorem also plugs directly into
(10) for general marked sources whenever its weighted innovation can be
bounded below. That latter lower bound is still open; we do not attach
the .4333221116640807 value to a source for which it was not proved.

## 6a. Independent finite-sum and normalization audit

The feedback researcher freshly read this entire proof, the full scalar
zero-first covariance theorem, and the full ordered next-gain theorem.
The finite-sum extension in Section 2 passes direct reconstruction.

More precisely, for a fixed finite family of bounded even old responses
H_l and bounded odd Lipschitz new responses psi_l, allow bounded
deterministic coefficients c_(l,i) and form
`sum_l c_(l,i) H_l(G_i) psi_l(Z_i)`.
Expand its covariance into the finitely many pairs (l,m). In each mixed
pair, both forests still have even old-branch count and odd new-branch
count. Equality of the two counts was never used: after whole matches,
a lone star still has an odd positive number of transported leaves and
an even number of old leaves. The exact BW bound in the original proof
therefore survives unchanged. Two-star and non-star pairings also keep
their two gains. Left/right multiplication by bounded diagonal coefficient
matrices preserves each operator and nuclear estimate. Whole matches give
the Gaussian cross term for H_l,H_m and psi_l,psi_m directly.

The fixed finite catalog bounded approximation handles all these terms
simultaneously. The clipped linear second summand is in that scope.
Removing its clip after n uses the previously proved uniform integrability
of averaged actual Z squares, not an unsupported raw row moment bound.
Finally mu_n tends to mu as one scalar number, so its weighted linear
coefficient error is uniform and has vanishing averaged L2 cost.
This proves exactly the covariance of the weighted sum; it does not
infer a nonsymmetric cross matrix by scalar polarization.

The positive fourth-power trace, ideal returned row variance cap, and
replacement by ideal square roots before multiplying the averaged mask
error all pass independently. The exact normalization is
`e(u+J sign(Bu))-e(u-J sign(Bu))=2 j` for the per-seed gain j and
`e(x)=x^T Bx/(2n)`. Half this difference is j, and the maximum absolute
endpoint energy dominates it. Thus no factor of two is lost.

## 7. A monotone erf-response corollary

The director proposed the following variant, and the independent auditor
reconstructed its covariance, constants, and cutoff passage before adding
this section. Keep exactly the source and mask hypotheses of Section 1,
but put `s=L tau` and use the bounded smooth odd increasing response

```math
\psi(z)=\operatorname{erf}\!\left(\frac{z}{\sqrt2s}\right),
\qquad C=H(G)\psi(Z),\quad u=f(G)+C,\quad J=1-|u|.
```

Then an explicit positive next-gain bound is

```math
\boxed{\quad
\liminf_n\frac1n E\sum_iJ_i|(Bu)_i|
 \ge\frac{\mu^3}{h}\,
       \frac{\sqrt3}{48\pi\sqrt\pi}\,L^{-9}>0.
\quad}                                                    (11)
```

This response also obeys `z psi(z)>=0`; no oscillating final response
is required. Combining (11) with (10) gives a strictly positive improvement
over the actual current fractional energy. It does not assign a previously
proved numerical certificate to this different scalar-source construction.

To verify (11), let `d_i^2=T_ii+s^2` and
`P=D_d^{-1}T D_d^{-1}`. The ideal response covariance is exactly

```math
K_\psi=\frac2\pi\arcsin^{\circ}(P).
```

Indeed psi(z) is the conditional expectation of sign(z+sN). In a product
of two conditional expectations use independent auxiliary noises, also
when the two root indices coincide. The Gaussian angle formula then
gives the displayed covariance, INCLUDING its diagonal; its diagonal
is not replaced by one. Equivalently the same identity follows from
the positive Hermite expansion of erf.

The exact first-history linear coefficient is
`a_i=mu_n sqrt(2/pi)/d_i`, with the same mu_n from Section 2.
The finite-sum covariance proof gives the corresponding returned
residual covariance and vanishing averaged absolute cross with its
coherent partner. Replacing mu_n by mu in the ideal matrix is harmless.
Since the cubic arcsin coefficient is 1/6, the ideal residual satisfies

```math
K_{\rm res}\succeq\frac{\mu^2}{3\pi}
                       D_d^{-3}T^{\circ3}D_d^{-3}.
```

Now `d_i^2<=2L^2 tau^2`, so the same trace argument as in Section 3 gives

```math
\frac1n\operatorname{Tr}(QK_{\rm res})
 \ge\frac{\mu^2}{3\pi nL^2\tau^2}
           \sum_i\frac{T_{ii}^4}{d_i^6}
 \ge\frac{\mu^2}{24\pi L^8}=:a_{\rm erf}.
```

The normalized matrix P has diagonal at most 1/2 and operator norm at
most L^2. The positive arcsin Schur series therefore obeys

```math
\left\|\frac2\pi\arcsin^{\circ}(P)\right\|_{op}
 \le\frac2\pi L^2\frac{\arcsin(1/2)}{1/2}
 =\frac23L^2.
```

Schur multiplication by K_H has norm at most h^2, and
`0<=K_res<=K_C` in positive order. Thus each ideal returned variance
is at most `M_erf=(2/3)h^2L^2`, using the unit row norm of B.

Finally for every sigma>=0,

```math
E\left|\operatorname{erf}
       \left(\frac{\sigma N}{\sqrt2s}\right)\right|
 =\frac2\pi\arctan(\sigma/s).
```

This follows by writing the left side as
`E sign(N) sign(sigma N+sN')` and integrating the rotationally invariant
two-dimensional Gaussian angle. Since sigma_i=sqrt(T_ii)<=s, the result
is at most 1/2. Local old/noise separation therefore gives
`E J_i>=mu/2-epsilon_i` with averaged absolute error tending to zero.
Use the capped ideal square roots before multiplying these errors, as
in Section 4. The ordered scalar regression cutoff then gives
`sqrt(2/pi)(mu/2) a_erf/sqrt(M_erf)`, which simplifies exactly to (11).

All approximation and variance-cutoff limits remain those of Sections
2--5. This is still a scalar zero-first source theorem, not a marked
f(G,Y) closure or a pointwise innovation assertion.
