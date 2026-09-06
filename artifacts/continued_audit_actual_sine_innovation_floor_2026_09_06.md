# A positive actual innovation floor for one explicit first-stage choice

Date: 2026-09-06. This is a SPECIAL bounded zero-first sine-response
construction on actual hollow signings. It does not assert a floor for
arbitrary first-stage choices, improve the preserved universal constant,
or establish convergence.

## 1. Explicit admissible source and exact fields

Let B=A/sqrt(n-1) be symmetric hollow with off-diagonal signs divided by
sqrt(n-1), and suppose ||B||op<=L with fixed L>=1. Set Q=B^2. Choose

```math
\alpha=e^{3/2}/2,\qquad
f(g)=\frac{\sin g-\alpha\sin(2g)}{4(1+\alpha)}.
```

Then f is smooth bounded odd, ||f||infinity<=1/4, and
`E[N f(N)]=0`, since `E[N sin(cN)]=c exp(-c^2/2)`.
Its Gaussian variance tau^2 is strictly positive. More generally the
proof below permits any fixed nonzero bounded odd zero-first scalar f,
with norm at most F<1, and a constant h>0 with `s=1-F-h>0`.

Write the normalized Hermite expansion and its exact ideal matrices as

```math
f=\sum_{p\ge3,\ p\ odd} f_p h_p,\qquad
R_f=\sum_p f_p^2 Q^{\circ p},\qquad T=B R_f B,
\qquad \tau^2=\sum_p f_p^2.
```

The ACTUAL first stage is

```math
G=BS,\quad Z=Bf(G),\quad C=h\sin(tZ),\quad u=f(G)+C,
\quad J=1-|u|\ge s.
```

This is inside the previously proved first marked-history theorem with
no dependence on old Y, zero coherent source coefficients, V=0, and
constant mask H=h. Its exact coherent conditional mean is c0=0. Put

```math
e_i=\exp(-t^2T_{ii}/2),\qquad a_i=ht e_i,
\quad L_0=B D_aZ,\qquad \eta=B(C-D_a Z).
```

The symbol L_0 denotes the coherent field, to distinguish it from the
operator cap L. Thus exactly `Bu=Z+L_0+eta`. Define the actual statistics
`v_i=E eta_i^2` and `rho_i=E eta_i (L_0)_i`.

## 2. Required full covariance theorem and the linear-channel extension

The input is the already independently audited theorem
`continued_feedback_zero_first_masked_covariance_2026_09_06.md`, Section
1, equations (1)--(2), with constant H. It states normalized NUCLEAR-norm
comparison of the full covariance of any bounded odd Lipschitz function
of the actual Z with its Gaussian vector of covariance T. It explicitly
licenses bounded-op transport by B. This is not merely a local CLT.

Its Section 6 explicitly permits bounded deterministic row-dependent
polynomial coefficients and a fixed finite catalog of approximating
responses. Uniformly quantize `a_i in [0,ht]` into finitely many bins.
Apply that same proof to the row responses
`h sin(tz)-a_i clip(z,-K,K)` and to the corresponding sums, holding
the bins and K fixed before n. Refine the coefficient bins afterwards.
The Lipschitz and L2 estimates are uniform over this finite catalog.

To include the linear response, let K tend to infinity AFTER n.
The actual Z has uniformly integrable
averaged squares: approximate f at a fixed Hermite-polynomial stage;
its transported squarefree polynomial has uniformly bounded row moments;
bounded-op transport makes the approximation error arbitrarily small
in limiting averaged L2. The inequality splitting Z into this polynomial
and its L2 error proves uniform averaged square tails. The ideal Gaussian
rows have uniformly bounded variances, as shown below, and hence uniform
square tails as well. Finally

```math
\|E[XY^T]\|_*\le
 \sqrt{E\|X\|^2\,E\|Y\|^2}
```

converts all these L2 errors into normalized nuclear errors. This proves
the covariance comparison for C, D_aZ, and C-D_aZ. Polarization gives
the SYMMETRIC cross matrix with the actual row coefficient D_a inside
the compared response. Polarizing two scalar responses alone would not
identify a nonsymmetric cross covariance, particularly when a_i varies.
No such full nonsymmetric claim is needed or made here.

For a Gaussian vector N_T of covariance T, exact Gaussian integration
gives `Cov(h sin(tN_T),N_T)=D_a T`. Consequently, with entrywise sinh,

```math
K_{\rm res}=h^2D_e[\sinh(t^2T)-t^2T]D_e
 =h^2\sum_{p\ge3,\ p\ odd}\frac{t^{2p}}{p!}
                   D_e T^{\circ p}D_e\succeq0.
```

The actual residual C-D_aZ has covariance K_res plus a nuclear o(n)
error. Its symmetric cross covariance with D_aZ is nuclear o(n), by
expanding the three compared covariance matrices and using Gaussian
orthogonality to the linear channel. Therefore

```math
\frac1n\|\operatorname{Cov}(\eta)-B K_{\rm res}B\|_*\to0,
\qquad \frac1n\sum_i|\rho_i|\to0.                    (1)
```

The second conclusion uses bounded-op transport and duality against
arbitrary bounded diagonal sign matrices. If M is the nonsymmetric
cross covariance, each real quadratic form `b_i M b_i^T` equals
`b_i(M+M^T)b_i^T/2`; the symmetric estimate therefore suffices. A signed
average alone would not suffice for the regression quotient below.

## 3. Positive variance from actual B identities

The correlation Schur multiplier inequality gives

```math
0\preceq R_f\preceq L^2\tau^2 I,
\qquad T\preceq L^2\tau^2 Q,
\qquad T_{ii}\le L^2\tau^2.                           (2)
```

The last bound uses that EVERY row of B has Euclidean norm one; it is
stronger than the crude operator bound `||T||<=L^4 tau^2`.
Furthermore, odd source degrees give the exact trace inequality

```math
\operatorname{Tr}T=\operatorname{Tr}(QR_f)
 =\sum_p f_p^2\sum_{i,j}Q_{ij}^{p+1}\ge n\tau^2.      (3)
```

All powers p+1 are even and Q_ii=1. No variance-floor assumption or
leverage theorem is required for (3).

Set `w=t^2 L^2 tau^2`. Keeping only the positive cubic term in K_res,
then using `Q>=T/(L^2 tau^2)`, yields

```math
\begin{aligned}
\frac1n\operatorname{Tr}(QK_{\rm res})
&\ge\frac{h^2t^6}{6nL^2\tau^2}
       \operatorname{Tr}(T D_e T^{\circ3}D_e)\\
&=\frac{h^2t^6}{6nL^2\tau^2}
       \sum_{i,j}e_i e_j T_{ij}^4\\
&\ge\frac{h^2t^6 e^{-w}}{6nL^2\tau^2}\sum_iT_{ii}^4\\
&\ge \frac{h^2t^6\tau^6e^{-w}}{6L^2}
       =:a_*>0.                                      (4)
\end{aligned}
```

The last line is Jensen applied to (3). The common sign of the cubic
sine coefficients is essential: it makes every summand in the second
line nonnegative. This argument must not be applied to arbitrary
root-dependent response coefficients.

An explicit ideal ROW variance upper bound is also available. For every
p, the positive Schur multiplier bound gives
`||T^{circ p}||<=max_i Tii^{p-1} ||T||`. Summing the sine series and using
(2) gives

```math
\|K_{\rm res}\|_{op}
 \le h^2 L^2[\sinh(w)-w]=:M_*.
```

Because rows of B have norm one, every diagonal
`v_i^*=(B K_res B)_ii` is at most M_*. Combining this with (1) and (4),

```math
\liminf_n\frac1n\sum_i\sqrt{v_i}
 \ge\frac{a_*}{\sqrt{M_*}}.                           (5)
```

Indeed `sqrt(v_i^*)>=v_i^*/sqrt(M_*)`, and the average difference of the
actual and ideal square roots is bounded by the square root of their
average absolute diagonal error. This avoids inferring a square-root
floor from an averaged variance floor alone.

## 4. Actual second-query gain and explicit specialization

Apply the audited scalar next-gain corollary in
`continued_audit_ordered_bounded_next_gain_2026_09_06.md`, Section 7.
At each fixed cutoff v_i>=epsilon, (1) implies

```math
\frac1n\sum_{v_i\ge\varepsilon}
 \left|\left|\sqrt{v_i}+\frac{\rho_i}{\sqrt{v_i}}\right|
                                      -\sqrt{v_i}\right|\to0.
```

The omitted average square root is at most sqrt(epsilon). Since J>=s,
first taking n to infinity and then epsilon down to zero gives the
positive actual lower bound

```math
\liminf_n\frac1n E\sum_iJ_i|(Bu)_i|
 \ge s\sqrt{\frac2\pi}\frac{a_*}{\sqrt{M_*}}
 =s h\sqrt{\frac2\pi}
   \frac{w^3e^{-w}}{6L^9\sqrt{\sinh(w)-w}}.             (6)
```

For the explicit f in Section 1 take h=1/4, s=1/2, and
`t=1/(L tau)`, so w=1. Then the bound is

```math
\boxed{\quad
\frac{\sqrt{2/\pi}\,e^{-1}}
 {48\sqrt{\sinh(1)-1}}\,L^{-9}>0.
\quad}                                                (7)
```

Every parameter is fixed before n. The endpoint sign is applied only
after the actual u is formed, and exact multilinear rounding makes this
a lower certificate for `max_x |x^T Bx|/(2n)`. The constant is deliberately
small and does not improve 0.4333221116640807. The useful conclusion is
structural: cancellation `v_i+rho_i=0` cannot eliminate the entire next
innovation for this explicit admissible choice. It may still occur for
other choices or on some roots; no general pointwise floor is asserted.

## 5. Reproducible finite stress tests, separate from the proof

The companion script
`computations/continued_audit_sine_innovation_replay_2026_09_06.py`
constructs the integer-verified Steiner signings, evaluates the exact
Gaussian sine-source covariance formula, checks the deterministic matrix
inequalities numerically, and samples the actual Boolean source and both
returns. Its saved result is
`computations/results/continued_audit_sine_innovation_replay_2026_09_06.json`.

For cap L=2 and 4,000 independent seeds at each order 28,120,496, the
estimated actual innovation variances are respectively
`0.00010668, 0.00012321, 0.00012619`; the ideal variances are
`0.00014513, 0.00012864, 0.00012737`. The signed mean cross covariances
are `0.00032024, 0.00003224, 0.00000610`, with Monte Carlo standard
errors `0.00002288, 0.00000878, 0.00000395`. These finite values are
consistent with, but do not establish a rate for, the proved limits.
The saved average absolute empirical cross covariances include sampling
noise; the script also records the rowwise standard-error scale rather
than presenting their positive bias as a limiting covariance.

A separate 70-digit outward interval check of the coefficient in (7)
gives `0.01460951544402038782864...`. It uses interval exponentials to
form sinh(1) and makes no optimization claim. All matrix and sampling
checks are explicitly labeled finite floating diagnostics; the proof
of (4)--(7) is the analytic argument above.

## 6. Independent audit of the feasible-mask extension without strict slack

The director's complete
`continued_director_masked_innovation_floor_2026_09_06.md` was freshly
read after the constant-mask proof. Its broader theorem passes at the
stated scalar zero-first scope, including the actual remaining-slack
mask `H=1-|f|` when its Gaussian mean is positive.

The needed covariance extension is a FINITE SUM of even-old times
odd-noise responses. Expanding both roots' forests preserves exactly
the even old count and odd noise count used in the original lone-star
classification. Bounded deterministic row coefficients and fixed finite
catalogs preserve every matrix estimate. Thus this extension follows
directly from the mixed forest proof; it is not an inference of an
unrestricted cross covariance from scalar polarization.

The exact coefficient uses `mu_n=E H(G_i)`, independent of i because
all rows are the same normalized signed-sum law. Bounded Gaussian-a.e.
continuity gives mu_n tending to mu. Its uniform deterministic coefficient
error is harmless in the weighted L2 and nuclear estimates. The ideal
mask covariance decomposes as `mu^2 J` plus positive even Hermite terms.
Subtracting only its constant-mask linear noise term leaves a positive
residual, bounded below by the constant-mask cubic term with h^2 replaced
by mu^2. Consequently the director's constants are exactly

```math
a_*=\frac{\mu^2 w^3e^{-w}}{6L^8},\qquad
M_*=h^2L^2\sinh(w),\qquad h=\|H\|_\infty.
```

There is no need for deterministic positive slack. Actual feasibility
gives `J>=H(G)(1-|sin(tZ)|)`. Local old/noise separation and
`|sin z|<=|z|` yield the row expectation lower bound
`mu(1-sqrt(2w/pi))` with averaged absolute error tending to zero.
Before multiplying this error by the innovation standard deviation,
replace the latter by its ideal counterpart, capped by sqrt(M_*).
Normalized nuclear diagonal control makes the replacement's average
absolute error vanish. This is sufficient for the weighted square-root
bound and does not multiply two uncontrolled averaged errors.

The ordered scalar cutoff then removes rho/sqrt(v), using the average
ABSOLUTE symmetric-cross estimate, and produces exactly

```math
\sqrt{\frac2\pi}\left(1-\sqrt{\frac{2w}{\pi}}\right)
 \frac{\mu^3}{h}
 \frac{w^3e^{-w}}{6L^9\sqrt{\sinh w}},\qquad 0<w\le1.
```

All factors are strictly positive under the declared assumptions. The
energy difference of the two final cube endpoints is twice the next
cross gain under `e(x)=x^T Bx/(2n)`, giving the exact original-signing
normalization. This proves a positive next innovation for that whole
scalar feasible-mask class, not an improvement over the preserved
universal lower bound, an additive improvement over the first-stage
certificate, or an extension to arbitrary marked f(G,Y).

The director subsequently sharpened the sine result to every fixed w>0.
This also passes: the absolutely convergent cosine series of |sin x|
has constant term 2/pi and all nonconstant coefficients negative.
Gaussian averaging makes all cosine expectations nonnegative. Hence
the mean slack factor is uniformly at least `mu(1-2/pi)` for all w>0;
the earlier small-w bound remains optionally sharper near zero.

The added finite energy-improvement statement also passes. For any cube
point u let `v=(1-|u|)sign(Bu)` and `j=u^T Bv/n`. The two endpoints
give `j<=Lambda<=L/2`. Choose `a=j/L` and the sign of e(u); the single
feasible point `u+sign(e(u)) a v` has absolute energy at least
`|e(u)|+aj-a^2Lambda>=|e(u)|+j^2/(2L)`. Averaging and Jensen give the
director's operational gap above the current fractional energy. This
is distinct from raising the preserved universal constant.

The monotone erf-response variant was independently derived and, at the
director's request, appended as Section 7 of that artifact. Its cubic
trace floor is `mu^2/(24 pi L^8)` and ideal row cap is
`(2/3)h^2 L^2`. Its actual mean slack is at least mu/2 in the averaged
comparison. The resulting next-gain constant is
`(mu^3/h) sqrt(3)/(48 pi sqrt(pi)) L^-9`. Independent auxiliary Gaussian
noises are used in the covariance product even on its diagonal; this
is why the normalized arcsin matrix has diagonal at most 1/2 rather
than one. The exact first coefficient retains mu_n before its scalar
limit. No oscillating final response, pointwise variance floor, or
unrestricted nonsymmetric covariance claim is needed.
