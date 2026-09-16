# Fresh adversarial reconstruction of the Boolean BH refinement

Date: 2026-09-16 UTC. Scope: the proposed `r^{-5}` weighted inequality and its `m^{11/2}` ordinary BH consequence only. This is secondary analysis, **not a new bound for the quadratic signing quantity `M_n`**, and no novelty claim is made.

## Verdict and evidence

I reconstructed the argument from the [primary paper](https://arxiv.org/html/2609.12427v1) and the candidate derivation in Sections 2–3 of `artifacts/bh_2026_09_16_independent_paper_audit.md`, checking the inequalities rather than relying on its audit verdict. I found no gap in the refinement. In particular, its use of complex scalars costs no hidden factor in hypercontractivity, and its contraction holds uniformly over every retained degree, not merely those checked numerically.

The provided replay `computations/bh_2026_09_16_independent_audit_checks.py` was also run independently; all 31 exact assertions passed. The analytic steps, which that replay cannot certify, are reconstructed below. An optional stronger constant follows by closing the *squared* inequality; neither exponent changes.

## 1. Low levels: independently checked Joukowski extraction

Let `p` be a complex polynomial of degree at most `m`, with `|p(t)| <= 1` on the real interval `[-1,1]`. Then

\[
G(w)=w^m p((w+w^{-1})/2)
\]

is an ordinary polynomial: each term of degree `k <= m` has lowest power `m-k >= 0`. On `|w|=1`, its argument is real and in `[-1,1]`, so `|G(w)| <= 1`; the maximum principle gives this bound throughout the disk.

For arbitrary complex `z`, the two roots of `w^2-2zw+1=0` have product 1. Choose a nonzero root with `|w| <= 1`. Its reciprocal is the other root; hence

\[
|p(z)|\le |w|^{-m}
\le (|z|+\sqrt{|z|^2+1})^m
=\exp(m\operatorname{arsinh}|z|)
\le e^{m|z|}.
\]

Cauchy's coefficient estimate on radius `R=r/m` therefore gives

\[
|[t^r]p|\le (em/r)^r,\qquad 1\le r\le m.
\]

For a cube polynomial `f` normalized by `||f||_infinity <= 1`, the multilinear extension satisfies the same bound on `[-1,1]^n`: the absolute value of each complex affine coordinate slice is convex and maximized at an endpoint. Apply the preceding extraction to `p_x(t)=f(tx)` at each vertex. Thus, writing

\[
q_r=\frac{2r}{r+1},\quad L_r(f)=\|\widehat f^{=r}\|_{q_r},\quad
W_s(f)^2=\sum_{r=1}^m r^{-2s}L_r(f)^2,
\]

and taking any finite homogeneous Boolean BH constants `C_1,C_2,C_3`,

\[
H(f)^2:=\sum_{r=1}^{\min(3,m)}r^{-5}L_r(f)^2
\le D_0^2m^6,\qquad
D_0^2=\sum_{r=1}^3 C_r^2(e/r)^{2r}/r^5.
\]

Only these three old fixed-degree BH constants are inputs. If their chosen source states a real-valued theorem, applying it separately to the real and imaginary parts gives complex constants with at most a factor 2 by coefficient-norm triangle inequality. No version of the new theorem is needed for these initial cases.

## 2. Finite-dimensional closure and the exact coefficient criterion

Fix `m,n` and let `M` be the supremum of `W_{5/2}(g)` over complex cube polynomials of degree at most `m` and supremum norm at most 1. This `M` is finite before proving any dimension-free estimate, since there are finitely many coefficients and each has modulus at most 1. Freezing any coordinates leaves the supremum norm at most 1 and degree at most `m`; re-embedding in the original cube by adding unused coordinates preserves every nonempty level norm. Thus every restriction has weighted norm at most this same `M`.

Randomly partition coordinates into two parts. On level `r >= 4`, keep `ceil(r/3) <= d <= floor(2r/3)`, with `e=r-d`. Put

\[
p_r=2^{-r}\sum_{d=\lceil r/3\rceil}^{\lfloor2r/3\rfloor}\binom rd,
\quad N_r=\lfloor2r/3\rfloor-\lceil r/3\rceil+1.
\]

The captured coefficient norm `C_r` obeys exactly `E C_r^{q_r}=p_r L_r(f)^{q_r}`. Probability-measure norm monotonicity then gives `L_r(f)^2 <= p_r^{-2/q_r} E C_r^2`. If `b_{d,e}` is the coefficient `ell_{q_r}` norm of one bidegree, counting-measure norm comparison gives

\[
C_r^2\le N_r^{1/r}\sum_{d+e=r\ \mathrm{retained}}b_{d,e}^2.
\]

For the matrix of coefficients of a retained bidegree define

\[
X_{d,e}=\|c\|_{\ell_A^{q_d}(\ell_B^2)},\qquad
Y_{d,e}=\|c\|_{\ell_B^{q_e}(\ell_A^2)},\qquad\theta=d/r.
\]

The identities `1/q_r=theta/q_d+(1-theta)/2=theta/2+(1-theta)/q_e` permit first row interpolation and then outer Hölder. The remaining norm is `ell_A^2(ell_B^{q_e})`, which is **at most** `ell_B^{q_e}(ell_A^2)` by Minkowski, as `q_e <= 2`. Therefore `b_{d,e} <= X_{d,e}^theta Y_{d,e}^{1-theta}`. In particular the claimed mixed-norm direction is correct.

Write `rho_d=sqrt((d-1)/(d+1))`, `Xhat=rho_d^e X`, `Yhat=rho_e^d Y`, `U=Xhat/d^{5/2}`, and `V=Yhat/e^{5/2}`. All retained `d,e` are at least 2. Direct substitution gives

\[
\frac{b_{d,e}}{r^{5/2}}\le D(d,e)h(\theta)^{5/2}U^\theta V^{1-\theta},
\]

where

\[
D(d,e)=\left[\frac{d+1}{d-1}\frac{e+1}{e-1}\right]^{de/(2r)},\qquad
h(\theta)=\theta^\theta(1-\theta)^{1-\theta}.
\]

Define the exact prefactor `K_{r,d}=p_r^{-1/q_r}N_r^{1/(2r)}D(d,e)h(theta)^{5/2}`. Squaring and applying weighted AM–GM yields

\[
T(f)^2\le\mathbb E\sum_{(d,e)\ \mathrm{retained}}
K_{r,d}^2\bigl(\theta U_{d,e}^2+(1-\theta)V_{d,e}^2\bigr).
\]

The sufficient pairwise criterion is precisely

\[
K_{r,d}\sqrt{2\max(\theta,1-\theta)}\le c<1.
\]

Indeed it gives `K_{r,d}^2 theta <= c^2/2` and `K_{r,d}^2(1-theta) <= c^2/2` **simultaneously for each pair**. There is no illicit averaging of the AM–GM weights or interchange of separate suprema.

## 3. Full-row bound, including complex scalars

For each fixed `d >= 2`, define the full complementary polynomial

\[
R_A(y)=\sum_{B\subseteq J}c_{A,B}y^B,\qquad |A|=d.
\]

For `z_{A,e}=(sum_{|B|=e}rho_d^{2e}|c_{A,B}|^2)^{1/2}`, the relevant Minkowski inequality is

\[
\|z\|_{\ell_e^2(\ell_A^{q_d})}\le
\|z\|_{\ell_A^{q_d}(\ell_e^2)}.
\]

Raising to power `q_d` identifies this with the triangle inequality in `ell_{2/q_d}`, so its direction is independently established. Adding complementary degrees 0 and 1 is allowed on the right side. Parseval followed by `q_d -> 2` real-noise hypercontractivity gives

\[
\left(\sum_{e\ge2}\widehat X_{d,e}^2\right)^{1/2}
\le\left(\sum_{|A|=d}\|T_{\rho_d}R_A\|_2^{q_d}\right)^{1/q_d}
\le\left(\mathbb E_y L_d(f(\cdot,y))^{q_d}\right)^{1/q_d}
\le\left(\mathbb E_y L_d(f(\cdot,y))^2\right)^{1/2}.
\]

For complex `R_A`, no extra hypercontractive constant is needed: the real noise operator has a positive probability kernel, so `|T_rho R_A| <= T_rho |R_A|`; apply real hypercontractivity to the real nonnegative function `|R_A|`. Its possibly high degree is irrelevant because the inequality holds for every cube function. Restoring the `d^{-5}` weight and summing gives `sum_{d,e>=2}U_{d,e}^2 <= M^2`; exchanging the parts gives the same bound for `V`.

The row must be kept in full. Replacing it by a homogeneous complementary level and then claiming its supremum is bounded by 1 would be invalid. The actual proof does not make that replacement.

## 4. Uniform contraction: independent arithmetic reconstruction

The function `k log((k+1)/(k-1))` decreases for `k>1`, as follows by expanding it as `2 sum_{j>=0}[(2j+1)k^{2j}]^{-1}`. Consequently `D <= 3` for `d,e >= 2`, and `D <= 2^{3/2}` for `d,e >= 3`.

For `r=4`, the only split is `(2,2)`. Thus `p_4=3/8`, `N_4=1`, `D=3`, `h=1/2`, and the sufficient factor is exactly

\[
c_4=3(8/3)^{5/8}2^{-5/2}=(27/32)^{1/8}<1.
\]

For `r=5,6,7`, use `D <= 3`, exact binomial probabilities and block counts, and the extreme retained ratio. Independently factoring these expressions gives respectively

\[
c_5^5=2^{17}3^{15}/5^{18},\qquad
c_6^{12}=2^{67}/(3^{23}5^{14}),\qquad
c_7^{14}=2^{110}3^{29}/(5^8 7^{50}).
\]

Each is less than `(9/10)` to its corresponding power, and `(9/10)^8 < 27/32`; these are exact rational comparisons.

For `r >= 8`, both split degrees are at least 3. Convexity and symmetry give `h <= (4/27)^{1/3}` and `max(theta,1-theta) <= 2/3`. Also `N_r <= r+1` and `log(r+1)/r` decreases, yielding `N_r^{1/(2r)} <= 9^{1/16}`. Finally `1/q_r <= 9/16`.

The finite rational checks `8 <= r <= 31` give `p_r >= 21/32`, with equality at `r=10`. For all `r >= 32`, Hoeffding gives `p_r >= 1-2e^{-r/18} >= 1-2e^{-16/9} > 21/32`; the strict inequality is certified by the positive Taylor truncation

\[
e^{16/9}>\sum_{j=0}^{5}(16/9)^j/j!=5189207/885735>64/11.
\]

Hence the tail factor is bounded by

\[
c_{\rm tail}=(32/21)^{9/16}9^{1/16}2^{3/2}(4/27)^{5/6}\sqrt{4/3}.
\]

Factoring gives `c_tail^{48}=2^{335}/(3^{165}7^{27})`; the exact integer inequality `2^{365}<3^{183}7^{27}` is equivalent to `c_tail<c_4`. Thus all retained bidegrees satisfy the criterion with `c=c_4=0.978986545550759...`. Combining Sections 2–3 now proves `T(f)^2 <= c_4^2 M^2`.

## 5. Conclusion and an optional constant improvement

The advertised closure `W_{5/2} <= H+T` proves `M <= D_0 m^3/(1-c_4)`, so the candidate theorem is valid as stated. A tighter closure is immediate from the *exact* disjoint-level identity

\[
W_{5/2}(f)^2=H(f)^2+T(f)^2\le D_0^2m^6+c_4^2M^2.
\]

Taking the finite-dimensional supremum and rearranging gives

\[
\boxed{W_{5/2}(f)\le C_{\rm sq}m^3\|f\|_\infty,\qquad
C_{\rm sq}=\frac{D_0}{\sqrt{1-(27/32)^{1/4}}}.}
\]

The multiplier `C_sq/D_0` is approximately `4.9037672937964345`, compared with approximately `47.588558198060895` for the stated linear closure. This is only a constant improvement and does not challenge any claim in the candidate audit.

The ordinary BH consequence follows from monotonicity of counting-measure coefficient norms, followed by comparison across the `m` levels:

\[
\|\widehat f\|_{q_m}
\le\|f\|_\infty+m^{1/(2m)}\Big(\sum_{r=1}^mL_r(f)^2\Big)^{1/2}
\le (1+2C_{\rm sq})m^{11/2}\|f\|_\infty.
\]

Here `m >= 1` and `m^{1/(2m)} <= 2`; the constant coefficient is included only in the ordinary BH norm. When `m < 4`, the low-level estimate already covers the whole weighted norm. All constants are independent of dimension, though not numerically specified without selecting the three old fixed-degree BH constants. No assertion about a competitive quadratic signing constant follows from this refinement.
