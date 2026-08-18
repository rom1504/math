# An exact square-polynomial carrier from `L^4` tail approximation

Status: **proved**.  A Hahn--Banach argument, the scalar Hamming-cube Riesz
transform inequality, and the exact scalar tail-space heat estimate give a
dimension-free `L^4` approximation of the square root of every weak-coordinate
row density by a bounded-degree Walsh polynomial.  Squaring and normalizing
then gives an exactly nonnegative polynomial density with a fixed `L^2` bound
and `L^2` error tending to zero.  Consequently both the actual pressure and
the row entropy are recovered uniformly.

This closes the entropy blockade in
[`actual_child_square_root_polynomial_carrier_audit.md`](actual_child_square_root_polynomial_carrier_audit.md).
It also removes the positive-part Fourier leakage of the earlier normed
carrier: the new density has literal Walsh degree at most `2d`.  It does not
compress the exponentially large cross-row response tensor or give a
child-only method for finding the approximating coefficients.

The two imported results are used with their exact hypotheses:

1. Eskenazis--Ivanisvili, *Polynomial inequalities on the Hamming cube*,
   Theorem 11, gives the scalar `L^p` heat-semigroup bound on Walsh tail
   spaces for every `1<p<infinity`.
   [Primary source](https://arxiv.org/abs/1902.02406)
2. Their Theorem 42 records Lust-Piquard's dimension-free scalar Riesz
   transform inequality for `p>=2`; at `p=4` it gives the direction used
   below.  [Original Lust-Piquard paper](https://doi.org/10.1006/jfan.1997.3217)

## 1. A fractional inverse estimate on the `L^(4/3)` tail space

Let `Delta` be the cube number operator,

```math
\Delta\chi_S=|S|\chi_S,                            \tag{L4.1}
```

and let `V_d` be the span of Walsh characters of degree at most `d`.

**Lemma L4.1 (tail fractional inverse).**  There is a universal constant
`C_tail` such that, on every Boolean cube and for every real function `h`
whose Walsh coefficients vanish on levels `0,1,...,d`,

```math
\boxed{
\|\Delta^{-1/2}h\|_{4/3}
\le C_{\rm tail}(d+1)^{-1/3}\|h\|_{4/3}.}         \tag{L4.2}
```

Here `Delta^(-1/2)` is zero on constants; `h` has no constant part, so
there is no ambiguity.

*Proof.*  For `p=4/3`, the angle in Eskenazis--Ivanisvili Theorem 11 is

```math
\theta_p=2\arcsin {2\sqrt{p-1}\over p}={2\pi\over3},
\qquad {\pi\over\theta_p}={3\over2}.              \tag{L4.3}
```

Their tail estimate therefore says that a function supported on levels at
least `D=d+1` obeys

```math
\|e^{-t\Delta}h\|_{4/3}
\le\phi(t)^D\|h\|_{4/3},                           \tag{L4.4}
```

where

```math
\phi(t)=
{(1+e^{-t})^{3/2}-(1-e^{-t})^{3/2}
 \over
 (1+e^{-t})^{3/2}+(1-e^{-t})^{3/2}}
={1-y(t)\over1+y(t)},
\qquad y(t)=\tanh(t/2)^{3/2}.                     \tag{L4.5}
```

For `0<=t<=1`, `tanh(t/2)>=t/4`, and hence

```math
\phi(t)\le e^{-2y(t)}\le e^{-t^{3/2}/4}.          \tag{L4.6}
```

For `t>=1`, apply (L4.4) successively on unit time intervals and use
ordinary `L^(4/3)` contraction on the final interval.  With
`c_1=-log(phi(1))/2>0`, this gives

```math
\|e^{-t\Delta}h\|_{4/3}
\le e^{-c_1Dt}\|h\|_{4/3}.                        \tag{L4.7}
```

(Increasing the harmless constant covers `1<=t<=2`.)  The spectral
identity

```math
\Delta^{-1/2}h={1\over\sqrt\pi}
 \int_0^\infty t^{-1/2}e^{-t\Delta}h\,dt          \tag{L4.8}
```

and Minkowski's inequality reduce the claim to a scalar integral.  On
`[0,1]`, the substitution `s=Dt^(3/2)/4` gives

```math
\int_0^1t^{-1/2}e^{-Dt^{3/2}/4}\,dt
\le C D^{-1/3}\int_0^\infty s^{-2/3}e^{-s}\,ds.
                                                               \tag{L4.9}
```

The integral on `[1,infinity)` from (L4.7) is exponentially small in `D`.
This proves (L4.2). `square`

The exponent is `1/3`, not `1/2`: it is exactly the result of integrating
the small-time tail decay `exp(-cDt^(3/2))` against `t^(-1/2)`.

## 2. Weak-coordinate square roots have bounded fractional gradient

Let `f>0` be a density relative to the fair cube law, and put `r=sqrt(f)`.
Assume

```math
\mathbb E_Uf=1,
\qquad \|f\|_2\le K,
\qquad
|\log f(x)-\log f(x^{(j)})|\le {A\over\sqrt n}    \tag{L4.10}
```

for every bit flip.

**Lemma L4.2 (dimension-free `L^4` fractional gradient).**  There is a
universal `C_R` such that

```math
\boxed{
\|\Delta^{1/2}r\|_4
\le C_R A e^{A/2}\|r\|_4
\le C_R A e^{A/2}K^{1/2}.}                        \tag{L4.11}
```

*Proof.*  With the convention
`partial_jr(x)=(r(x)-r(x^(j)))/2`, (L4.10) gives pointwise

```math
|\partial_jr(x)|
\le {e^{A/(2\sqrt n)}-1\over2}r(x)
\le {Ae^{A/2}\over4\sqrt n}r(x).                 \tag{L4.12}
```

Therefore

```math
\left(\sum_j|\partial_jr|^2\right)^{1/2}
\le {Ae^{A/2}\over4}r.                            \tag{L4.13}
```

The lower side of Lust-Piquard's `p=4` Riesz-transform equivalence is

```math
c_4\|\Delta^{1/2}r\|_4
\le\left\|\left(\sum_j|\partial_jr|^2\right)^{1/2}\right\|_4,
                                                               \tag{L4.14}
```

with a dimension-free `c_4>0`.  Combine (L4.13)--(L4.14) and use
`||r||_4=(E f^2)^(1/4)<=K^(1/2)`. `square`

The direction of (L4.14) is important.  The reverse Riesz inequality fails
below `p=2`, but the proof uses it only at `p=4`.

## 3. `L^4` approximation by low Walsh degree

**Theorem L4.3 (dimension-free square-root approximation).**  Under
(L4.10), for every `d>=0` there is a real Walsh polynomial `g_d in V_d`
such that

```math
\boxed{
\|r-g_d\|_4
\le C(A,K)(d+1)^{-1/3}.}                           \tag{L4.15}
```

*Proof.*  Hahn--Banach duality for distance to a closed subspace gives

```math
\operatorname {dist}_{L^4}(r,V_d)
=\sup\left\{
 \mathbb E_U(rh):h\perp V_d,\ \|h\|_{4/3}=1
 \right\}.                                       \tag{L4.16}
```

The annihilator consists exactly of functions whose Walsh spectrum is on
levels at least `d+1`.  Fourier expansion, followed by Holder, gives

```math
\begin{aligned}
|\mathbb E_U(rh)|
&=|\mathbb E_U{(\Delta^{1/2}r)
                 (\Delta^{-1/2}h)\}|\\
&\le\|\Delta^{1/2}r\|_4
      \|\Delta^{-1/2}h\|_{4/3}.                  \tag{L4.17}
\end{aligned}
```

Lemmas L4.1--L4.2 prove (L4.15).  All spaces and functions can be taken
real; alternatively, taking the real part of a complex approximant cannot
increase its distance from the real function `r`. `square`

There is no illicit identification of the `L^4` best approximant with the
orthogonal Walsh projection.  The theorem asserts existence of a best (or
arbitrarily near-best) polynomial through Banach-space duality.

## 4. Squaring gives a fixed-`L^2`, exact polynomial density

Put

```math
\varepsilon_d=\|r-g_d\|_4,
\qquad z_d=\mathbb E_Ug_d^2,
\qquad q_d={g_d^2\over z_d}.                       \tag{L4.18}
```

For all sufficiently large `d`, `z_d>0`.  The density `q_d` is nonnegative,
has mean one, and has exact Walsh degree at most `2d` (multiplication of
Walsh characters takes symmetric differences).

**Theorem L4.4 (`L^2` density and entropy recovery).**  There are constants
`K_1=K_1(K)<infinity` and `C_1=C_1(K)<infinity` such that, whenever
`epsilon_d<=1/2`,

```math
\boxed{
\|q_d\|_2\le K_1,
\qquad
\|q_d-f\|_2\le C_1\varepsilon_d.}                \tag{L4.19}
```

Consequently the dimension-free entropy modulus from Lemma LDC.2 gives

```math
\boxed{
|D(q_dU\Vert U)-D(fU\Vert U)|
\le\omega_{\max\{K,K_1\}}(C_1\varepsilon_d)
\longrightarrow0.}                               \tag{L4.20}
```

*Proof.*  Since `||r||_2=1` and `||r-g_d||_2<=epsilon_d`,

```math
(1-\varepsilon_d)^2\le z_d\le(1+\varepsilon_d)^2. \tag{L4.21}
```

Also

```math
\|g_d\|_4\le\|r\|_4+\varepsilon_d
\le K^{1/2}+\varepsilon_d.                        \tag{L4.22}
```

It follows that

```math
\|q_d\|_2={\|g_d\|_4^2\over z_d}
\le4(K^{1/2}+1)^2=:K_1.                           \tag{L4.23}
```

Next,

```math
\|g_d^2-r^2\|_2
\le\|g_d-r\|_4\|g_d+r\|_4
\le\varepsilon_d(2K^{1/2}+\varepsilon_d),        \tag{L4.24}
```

and

```math
|z_d-1|
\le\|g_d^2-r^2\|_1
\le\varepsilon_d(2+\varepsilon_d).               \tag{L4.25}
```

Combining (L4.21)--(L4.25) bounds

```math
\left\|{g_d^2\over z_d}-r^2\right\|_2
\le\|g_d^2-r^2\|_2
 +{|1-z_d|\over z_d}\|g_d^2\|_2
\le C_1(K)\varepsilon_d.                         \tag{L4.26}
```

Both densities now lie in one fixed `L^2` ball, so Lemma LDC.2 proves
(L4.20). `square`

This is exactly what the ordinary `L^2` Fourier projection of `sqrt(f)` did
not provide.  The stronger `L^4` approximation controls the `L^2` norm of
the squared polynomial and makes entropy uniformly continuous.
The spike from SQ.16 does not contradict this statement: if its exceptional
atom has fair mass `2^(-d)`, its square-root perturbation has `L^4` size of
order `2^(d/4)/sqrt(d)`, rather than tending to zero.  The old obstruction
separated Hellinger control from entropy precisely by escaping the `L^4`
ball now controlled here.

## 5. Uniform recovery of the actual row-product value

Define the child-independent exact square-polynomial carrier

```math
\mathcal S_{n,d,K_1}
=\left\{
 q=g^2:\ \deg g\le d,\ \mathbb E_Ug^2=1,
 \ \|g\|_4^2\le K_1
 \right\}.                                       \tag{L4.27}
```

Every `q` in this class is a probability density of literal Walsh degree at
most `2d` and satisfies `||q||_2<=K_1`.  It has
`sum_(a=0)^d binom(n,a)` real presentation coefficients before the two
norm constraints.

The presentation is also operational at fixed degree.  Parseval checks
`E g^2=1` from the coefficients, and the degree-`d` evaluation norm gives

```math
0\le q(b)=g(b)^2\le D(n,d):=\sum_{a=0}^d{n\choose a}.            \tag{L4.27a}
```

Thus rejection sampling from `U_n` has expected at most `D(n,d)` proposals.
The entropy integrand is bounded by `D(n,d) log D(n,d)` from above (and by
`e^(-1)` from below), so its row expectation is estimable with polynomially
many fair samples for fixed `d`.  As in Proposition LDC.4, the pressure mean
of any declared carrier product is estimable to additive `epsilon N` from a
number of complete bridge samples independent of the orders, apart from the
polynomial row-rejection cost.  None of these facts solves the global
coefficient optimization.

Let `L` be the exact bridge pressure and

```math
\mathcal F(P)=\mathbb E_PL+{1\over\lambda}
 \sum_iD(P_i\Vert U_n).                            \tag{L4.28}
```

Let `V^row` be its infimum over all row products and let `V^(d,sq)` be its
infimum when every row factor lies in (L4.27).

**Theorem L4.5 (exact-polynomial recovery of the actual product shadow).**
For fixed `beta,lambda>0`, there are a constant `K_1(beta,lambda)` and a
sequence `eta_d^sq->0`, independent of the child orders and their exact
minimizing signings, such that, for all sufficiently large `d`,

```math
\boxed{
0\le V^{(d,{\rm sq})}-V^{\rm row}
\le m\eta_d^{\rm sq}.}                            \tag{L4.29}
```

One may take

```math
\eta_d^{\rm sq}
\le C_{\beta,\lambda}(d+1)^{-1/3}
 +{1\over\lambda}
  \omega_{K_*}\left(C_{\beta,\lambda}(d+1)^{-1/3}\right),
                                                               \tag{L4.30}
```

where `K_*` is fixed.  With the explicit conservative modulus LDC.10a,
this is

```math
\eta_d^{\rm sq}=O_{\beta,\lambda}
 ((d+1)^{-1/12}\log(d+2)).                        \tag{L4.31}
```

*Proof.*  Every row factor `f_i` of a globally optimal actual-child product
shadow satisfies, uniformly,

```math
\|f_i\|_2\le e^{\lambda^2\beta^2/2},
\qquad
\operatorname {osc}_{b_j}\log f_i
\le {2\lambda\beta\over\sqrt n}.                \tag{L4.32}
```

Apply Theorems L4.3--L4.4 to each factor.  This gives
`q_(i,d) in S_(n,d,K_1)` and

```math
\|q_{i,d}-f_i\|_2
\le C_{\beta,\lambda}(d+1)^{-1/3}.                \tag{L4.33}
```

Replace the product factors sequentially.  Averaging `L` over every other
current factor gives a fair-row function `F_i` with one-bit oscillation at
most `2u`, where `u=beta/sqrt(N)`.  Cube Poincare and cancellation of the
constant give

```math
|\mathbb E_{q_{i,d}}F_i-\mathbb E_{f_i}F_i|
\le u\sqrt n\,\|q_{i,d}-f_i\|_2
\le\beta\|q_{i,d}-f_i\|_2.                       \tag{L4.34}
```

Equation (L4.20) controls the entropy change.  Sum over the `m` rows to get
the upper bound in (L4.29); the lower bound follows by restriction of the
competitor class. `square`

The same Gibbs identity as before yields the restricted projection estimate

```math
0\le\mathcal I_{d,{\rm sq}}^{\leftarrow}
      -\mathcal I^{\leftarrow}
\le\lambda m\eta_d^{\rm sq}.                     \tag{L4.35}
```

Thus every fixed extensive reverse-product or coherent-retuning rate is
visible at one fixed exact row Walsh degree.

## 6. Adversarial audit and SML effect

The candidate survives the following possible failure points.

1. **Wrong dual space:** the annihilator of `V_d` in `L^(4/3)` is exactly
   the Walsh tail above `d`; (L4.16) is standard Banach-space distance
   duality, not `L^2` orthogonal projection.
2. **Wrong Riesz direction:** the proof uses the valid `p=4` inequality
   `||Delta^(1/2)r||_4 <= C||grad r||_4`, never the false sub-`2`
   reverse inequality.
3. **Wrong tail exponent:** at `p=4/3`, `theta_p=2pi/3`; the heat exponent is
   `t^(pi/theta_p)=t^(3/2)`, and fractional integration gives `d^(-1/3)`.
4. **Normalization loss:** `z_d` stays in
   `[(1-epsilon_d)^2,(1+epsilon_d)^2]`; normalization cannot amplify the
   error once `d` is fixed sufficiently large.
5. **Entropy loss:** `L^4` approximation of the square root gives `L^2`
   approximation of the density and a uniform `L^2` bound on the new
   density.  This is precisely the missing hypothesis in the prior entropy
   audit.
6. **Hidden positivity leakage:** none.  Squaring gives an exact
   nonnegative polynomial of row Walsh degree at most `2d`.

The result is stronger than the positive-part carrier on the response side:
at fixed macroscopic accuracy, only pressure Walsh coefficients having row
degree at most `2d` can be queried.  It still does **not** bound the number of
rows participating in a coefficient.  The formal response tensor has up to

```math
\left(\sum_{a=0}^{2d}{n\choose a}\right)^m        \tag{L4.36}
```

entries.  Nor does the proof find `g_d` from a polynomial child statistic;
it begins with the unknown optimal factor and uses an existence theorem.

Therefore this is a real additional factor/response reduction and removes
the exact affine positive-part obstruction, but it is not Level 6.  The
strictly narrower missing lemma becomes:

> **Bounded-row-degree cross-row closure.**  At one fixed row degree,
> evaluate or bound the jointly optimized square-polynomial response from an
> actual-child state that does not store the full cross-row tensor, to
> macroscopic accuracy; or prove a scalable actual-minimizer obstruction to
> every such closure.

The unresolved information is now cross-row order, not row Walsh degree,
positivity, density collision norm, or row entropy.
