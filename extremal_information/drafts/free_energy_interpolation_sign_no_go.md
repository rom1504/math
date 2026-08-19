# Free-energy composition: a sublinear universality theorem and the Gaussian sign obstruction

**Status.**  Proof-ready task artifact.  The positive statement below is a
uniform `O(sqrt(N))` Rademacher--Gaussian replacement for the *actual*
bridge soft minimum, with an immediate arrow to `P_N`.  The negative
statement is a **sector-resolved** endpoint comparison theorem: on the natural
variance-preserving interpolation which transfers the radial variance
surplus of two children to the cross block, every generalized disorder
moment of order `s<1` moves in the wrong direction (up to one explicit
`O_beta(1)` diagonal correction).  This includes negative disorder
temperature `s=-lambda`, the quenched case `s=0`, and positive fractional
moments `0<s<1`.

The result does **not** prove an almost-subadditive recurrence for `P_N`.
It rules out the covariance-only Gaussian/Guerra implementation only after
the ordinary exponential sector has been resolved.  The physical augmented
`cosh` interpolation retains a joint sector sign and is not ordered by the
theorem.  A successful interpolation may therefore use a directional
property of that actual joint-sector Gibbs law, or a different endpoint
already supplied with an upper calibration to the child optima.

## 1. Normalizations and the required arrows

For a hollow symmetric signing `A=(a_ij)` of order `k`, put

```math
H_A(x)=\sum_{1\le i<j\le k}a_{ij}x_ix_j,
\qquad
\overline Z_A(t)=\mathbb E_x\cosh(tH_A(x)),              \tag{FI.1}
```

and

```math
F_k(t)=\min_A\log\overline Z_A(t),
\qquad
P_k(\beta)=F_k(\beta/\sqrt k).                          \tag{FI.2}
```

Let

```math
M_k=\min_A\max_x|H_A(x)|.
```

For every `beta>0`, the elementary soft/ground comparison is

```math
\boxed{
 {P_k(\beta)\over\beta k}
 \le {M_k\over k^{3/2}}
 \le {P_k(\beta)\over\beta k}+{\log2\over\beta}.}       \tag{FI.3}
```

Indeed, if `Q(A)=max_x|H_A(x)|` and `t=beta/sqrt(k)`, then
`log Zbar_A(t)<=tQ(A)`.  The two antipodal maximizing spin words give
`log Zbar_A(t)>=tQ(A)-k log2`.  Minimize over `A`.

Consequently, an all-comparable-order recurrence

```math
P_{m+n}(\beta)\le P_m(\beta)+P_n(\beta)+D_{m+n}(\beta), \tag{FI.4}
```

whose defect density is summable on dyadic shells (in particular,
`D_N=O_beta(N^(1-delta))` for any fixed `delta>0`) gives convergence of
`P_k(beta)/k`.  Equation (FI.3), first at fixed `beta` and then with
`beta->infinity`, makes the oscillation of `M_k/k^(3/2)` at most
`log2/beta`; hence it gives convergence of the extremal normalization.
This is the quantitative `P => M` obligation imposed on every auxiliary
quantity below.

Fix a split `m+n=N`, an orientation `epsilon in {+-1}`, and exact children
`A,D` minimizing at the contracted raw temperature `t=beta/sqrt(N)` (the
replacement theorem itself is valid for arbitrary children).  Put

```math
t={\beta\over\sqrt N},
\quad
L_\epsilon(B)=\log\mathbb E_{x,y}\cosh\{t(H_A(x)+
                 \epsilon H_D(y)+x^{\mathsf T}By)\}.     \tag{FI.5}
```

For `s!=0`, define the generalized bridge moment

```math
\Phi_s^{\rm Rad}
 ={1\over s}\log\mathbb E_{\epsilon,B}
                      e^{sL_\epsilon(B)},               \tag{FI.6}
```

where the bridge entries and `epsilon` are fair signs.  At `s=0`, use the
continuous extension `Phi_0^Rad=E L`.  Negative disorder temperature is
`s=-lambda<0`.

Every power mean lies above the minimum.  Parent minimization therefore
gives the direct, assumption-free arrow

```math
\boxed{P_N(\beta)\le\min_{\epsilon,B}L_\epsilon(B)
                    \le\Phi_s^{\rm Rad}\qquad(s\in\mathbb R).} \tag{FI.7}
```

Thus an upper estimate on (FI.6) at the child target, with a
power-saving error, would immediately be the desired recurrence.  No
separate approximation of the exact bridge minimum is required.

## 2. A uniform `O(sqrt(N))` Rademacher--Gaussian replacement

The replacement is stated in a form which includes the augmented cosh
partition in (FI.5).

### Lemma FI.1 (generalized-moment invariance)

Let `Omega` be finite, let `mu` be a positive measure on `Omega`, and let
`q_e:Omega->[ -1,1 ]`, `1<=e<=d`.  For `z in R^d`, set

```math
L(z)=\log\sum_{\omega\in\Omega}\mu(\omega)
 \exp\left\{h(\omega)+t\sum_{e=1}^dz_eq_e(\omega)\right\}. \tag{FI.8}
```

Let `R_e` be iid fair signs and `G_e` iid standard Gaussians.  For fixed
`s in R`, define

```math
\Phi_s(X)=
\begin{cases}
s^{-1}\log\mathbb E e^{sL(X)},&s\ne0,\\
\mathbb E L(X),&s=0.
\end{cases}                                             \tag{FI.9}
```

For every fixed `b,T<infinity`, there is an explicit finite constant
`C_(b,T)` such that, whenever `|s|<=b` and `0<=t<=T`,

```math
\boxed{|\Phi_s(R)-\Phi_s(G)|\le C_{b,T}\,d t^3.}         \tag{FI.10}
```

The constant is independent of `Omega,mu,h,d` and the feature geometry.

#### Proof

For one coordinate `u`, Gibbs differentiation gives

```math
|\partial_uL|\le t,
\qquad |\partial_u^2L|\le t^2,
\qquad |\partial_u^3L|\le8t^3.                          \tag{FI.11}
```

The last derivative is the third centered moment of a random variable in
`[-1,1]`.  If `f(u)=exp(sL(u))`, then

```math
|f'''(u)|
\le |s|t^3(8+3|s|+s^2)f(u).                            \tag{FI.12}
```

Moreover `f(u)/f(0)` lies between
`exp(-|s|t|u|)` and `exp(|s|t|u|)`.  Taylor expansion at zero, with the
matching first two moments of a fair sign and a standard Gaussian, gives

```math
|\mathbb E f(R)-\mathbb E f(G)|
\le c_{b,T}t^3 f(0),                                   \tag{FI.13}
```

where, for `|s|<=b`, one may take the right side constant to be

```math
{|s|(8+3b+b^2)\over6}
\left(e^{bT}+\mathbb E[|G|^3e^{bT|G|}]\right).          \tag{FI.14}
```

Both expectations in (FI.13) are at least

```math
k_{b,T}f(0),
\qquad
k_{b,T}=\min\{e^{-bT},\mathbb E e^{-bT|G|}\}>0.         \tag{FI.15}
```

Hence replacing one coordinate changes the logarithm of the generalized
moment by at most `|s| c'_(b,T)t^3`, for a finite `c'_(b,T)`.
Condition on all other coordinates and telescope through the `d`
coordinates.  Division by `|s|` cancels the factor `|s|` in (FI.14).
The limit `s=0` follows directly by applying the same matched-moment Taylor
argument to `L`, using (FI.11).  This proves (FI.10). `square`

Apply the lemma to (FI.5) after writing

```math
\cosh(tU)=\mathbb E_{\tau\in\{+-1\}}e^{t\tau U}.
```

The bridge features are `q_(ij)=tau x_i y_j`, so they have absolute value
one.  With `d=mn` and `t=beta/sqrt(N)`, Lemma FI.1 gives, for fixed `s,beta`,

```math
\boxed{
|\Phi_s^{\rm Rad}-\Phi_s^{\rm Gau}|
 \le C_{s,\beta}{\beta^3mn\over N^{3/2}}
 \le {C_{s,\beta}\beta^3\over4}\sqrt N.}              \tag{FI.16}
```

Combining (FI.7) and (FI.16), for every fixed negative disorder
temperature `lambda>0`,

```math
\boxed{P_N(\beta)\le\Phi_{-\lambda}^{\rm Gau}
                         +O_{\beta,\lambda}(\sqrt N).}   \tag{FI.17}
```

This is a genuine power-saving reduction for the actual cosh bridge law.
The missing issue is the *direction* of the Gaussian endpoint comparison.

## 3. Exact generalized Gaussian interpolation formula

Let `X_u(omega)`, `0<=u<=1`, be a differentiable centered Gaussian process
on a finite state space with covariance `C_u(omega,omega')`.  Put

```math
Z_u=\sum_\omega\mu(\omega)e^{h(\omega)+X_u(\omega)},
\qquad
\Psi_s(u)=
\begin{cases}
s^{-1}\log\mathbb E Z_u^s,&s\ne0,\\
\mathbb E\log Z_u,&s=0.
\end{cases}                                             \tag{FI.18}
```

Let `E_s<.>` denote Gaussian expectation tilted by `Z_u^s`, followed by
one or two independent Gibbs replicas.  Gaussian integration by parts
gives

```math
\boxed{
\Psi_s'(u)={1\over2}\mathbb E_s\left\langle
 \dot C_u(\omega^1,\omega^1)
 +(s-1)\dot C_u(\omega^1,\omega^2)
\right\rangle.}                                        \tag{FI.19}
```

The formula at `s=0` is its continuous limit.  Notice that negative
disorder temperature changes the replica coefficient from `-1` to
`-(1+lambda)`; it does not reverse its sign.

## 4. The variance-preserving child-to-cross sign theorem

The next comparison is independent of the deterministic Hamiltonian.  It
therefore applies with either ordinary sector of any pair of actual
optimizing children in the base term `h`.

Write

```math
q_x(x,x')={1\over m}\sum_{i=1}^m x_ix_i',
\qquad
q_y(y,y')={1\over n}\sum_{j=1}^n y_jy_j',               \tag{FI.20}
```

and define

```math
t^2={\beta^2\over N},
\qquad
\sigma_m^2=\beta^2\left({1\over m}-{1\over N}\right)
 ={\beta^2n\over mN},
\qquad
\sigma_n^2={\beta^2m\over nN}.                         \tag{FI.21}
```

These are exactly the variance surpluses between child-own scale and
parent scale.  Let

```math
X_{\rm int}(x,y)
=\sigma_m\sum_{i<k}g^m_{ik}x_ix_k
+\sigma_n\sum_{j<l}g^n_{jl}y_jy_l,                      \tag{FI.22}
```

```math
X_{\rm cross}(x,y)
=t\sum_{i,j}g^B_{ij}x_iy_j.                             \tag{FI.23}
```

All Gaussian variables are independent.  Finally set

```math
\widetilde X_{\rm int}=X_{\rm int}+{\beta\over\sqrt2}g_0, \tag{FI.24}
```

where `g_0` is a common scalar Gaussian.  It repairs only the missing
diagonal variance and shifts every generalized pressure by the explicit
constant `s beta^2/4`.

### Theorem FI.2 (all `s<1` move the wrong way)

For every finite base measure and every deterministic `h(x,y)`, let
`Psi_s^int` and `Psi_s^cross` be (FI.18) with Gaussian processes
`X_int` and `X_cross`, respectively.  Then, for every `s<1`,

```math
\boxed{
\Psi_s^{\rm cross}
 \ge \Psi_s^{\rm int}+{s\beta^2\over4}.}                \tag{FI.25}
```

Equivalently, after the harmless common-Gaussian normalization,

```math
\boxed{
\Psi_s^{\rm cross}\ge\widetilde\Psi_s^{\rm int}
\qquad(s<1).}                                           \tag{FI.26}
```

More precisely, along the linear covariance interpolation from
`Xtilde_int` to `X_cross`,

```math
{d\over du}\Psi_s(u)
={1-s\over2}{\beta^2mn\over2N}
 \mathbb E_s\langle(q_x-q_y)^2\rangle\ge0.             \tag{FI.27}
```

Thus:

- for `s=-lambda<0`, the cross soft minimum is no smaller than the
  internal endpoint apart from the explicit constant `lambda beta^2/4`;
- at `s=0`, the quenched cross pressure is no smaller, with no correction;
- for every positive fractional moment `0<s<1`, the inequality is strict
  in the same wrong direction whenever the two overlaps do not coincide;
- only `s>1` reverses the comparison, but that positive-moment endpoint
  has no upper calibration by child minimality.

#### Proof

Let `C_cross,C_int` denote the two covariance kernels.  Direct calculation
using

```math
\sum_{i<k}x_ix_i'x_kx_k'={m^2q_x^2-m\over2}
```

gives

```math
C_{\rm cross}={\beta^2mn\over N}q_xq_y,                 \tag{FI.28}
```

```math
C_{\rm int}
={\beta^2mn\over2N}(q_x^2+q_y^2)-{\beta^2\over2}.      \tag{FI.29}
```

The common scalar in (FI.24) adds `beta^2/2` to every covariance entry.
Consequently

```math
\boxed{
C_{\rm cross}-\widetilde C_{\rm int}
=-{\beta^2mn\over2N}(q_x-q_y)^2.}                      \tag{FI.30}
```

In particular, the self-covariance difference vanishes.  Insert (FI.30)
into (FI.19).  This gives (FI.27), and integration over `u in [0,1]`
gives (FI.26).  Since a common Gaussian `a g_0` changes (FI.18) by
`s a^2/2=s beta^2/4`, (FI.25) follows. `square`

The theorem holds for every split, not merely balanced splits.  On a
comparable window `m/N in [alpha,1-alpha]`, its square coefficient is
`Theta_(alpha,beta)(N)`, so an overlap mismatch of fixed mean size creates
a linear wrong-way increment.

### Actual-child strict witness

Take `m=n=2` and let both children be the one-edge positive signing.  They
are exact pressure minimizers at every temperature because all order-two
signings are switching equivalent.  At the internal endpoint, each
ordinary-sector two-spin Gibbs law has global spin symmetry.  If `q` is
the normalized overlap of two replicas, then, conditional on any finite
realized coupling `J`,

```math
\mathbb E\langle q\rangle=0,
\qquad
\mathbb E\langle q^2\rangle={1+\tanh^2J\over2}\ge{1\over2}. \tag{FI.31}
```

The two child systems are conditionally independent and each conditional
first overlap moment vanishes.  Hence, under every outer `Z^s` tilt,

```math
\mathbb E_s\langle(q_x-q_y)^2\rangle
=\mathbb E_s\langle q_x^2\rangle
 +\mathbb E_s\langle q_y^2\rangle\ge1.                 \tag{FI.31a}
```

The value is generally strictly larger than one when the realized
couplings are nonzero.  Since here
`beta^2mn/(2N)=beta^2/2`, (FI.27) yields at `u=0`

```math
\Psi_s'(0)\ge{(1-s)\beta^2\over4}>0
\qquad(s<1).                                            \tag{FI.32}
```

Thus the wrong sign is strict already for genuine optimizing children; it
is not an artifact of a surrogate or of a limiting argument.

## 5. Why biased-sign radial Gaussianization does not repair the endpoint

For the left child, let `xi_e` be independent signs with

```math
\mathbb E\xi_e=\sqrt{m/N}.
```

Then the full child-scale coefficient

```math
{\beta\over\sqrt m}a_e\xi_e                            \tag{FI.33}
```

has mean `t a_e` and variance `sigma_m^2` from (FI.21).  The analogous
identity holds on the right.  Thus (FI.22) is precisely the Gaussian
moment-matching surrogate for constant-density sign noise around the two
children.  On every comparable split, the same proof as Lemma FI.1 (now
using the uniformly bounded third centered moments of the biased signs)
replaces this biased-sign cloud by (FI.22) with `O_(alpha,beta,s)(sqrt N)`
error.

But every realization `A*xi` is another complete signing.  Exact child
minimality supplies a **lower** bound on its cosh pressure,

```math
\log\overline Z_{A\xi}(\beta/\sqrt m)\ge P_m(\beta),    \tag{FI.34}
```

not the upper bound needed on the right side of a parent recurrence.
Theorem FI.2 says that, for all min-favoring exponents `s<1`, transferring
this cloud variance to the cross block can only raise the generalized
pressure (modulo `O_beta(1)`).  Increasing negative disorder temperature
multiplies the same square by `1+lambda`; it strengthens rather than
removes the sign obstruction.

For `s>1`, (FI.27) reverses, but (FI.34) still points the wrong way and
positive moments emphasize high-pressure noisy children.  No bound of the
form `Psi_s^int<=P_m+P_n+o(N)` follows from actual minimality or from the
moment interpolation.

## 6. Exact radial-scaling audit

The same direction conflict can be seen without Gaussian disorder.  Put

```math
E_k={k\choose2},
\quad \rho=\tanh t,
\quad
W_A(\rho)={\overline Z_A(t)\over(\cosh t)^{E_k}},
\quad
r_k(t)=\log\min_AW_A(\tanh t).                          \tag{FI.35}
```

If `0<=rho_0<=rho_1` and `eta_e` are independent signs of mean
`rho_0/rho_1`, then direct expansion gives

```math
\boxed{
W_A(\rho_0)=\mathbb E_\eta W_{A\eta}(\rho_1).}          \tag{FI.36}
```

Taking minima proves

```math
r_k(t_0)\ge r_k(t_1)\qquad(0\le t_0\le t_1).           \tag{FI.37}
```

At a fixed common raw temperature, uniform bridge and relative-orientation
averaging gives

```math
r_{m+n}(t)\le r_m(t)+r_n(t).                            \tag{FI.38}
```

The two directions oppose at physical temperature.  If
`t_k=beta/sqrt(k)` and `t_N=beta/sqrt(N)`, then

```math
\boxed{
P_N(\beta)\le P_m(\beta)+P_n(\beta)
 +\kappa_{m,n}(\beta)+G_{m,n}(\beta),}                 \tag{FI.39}
```

where

```math
\begin{aligned}
G_{m,n}(\beta)
&=[r_m(t_N)-r_m(t_m)]+[r_n(t_N)-r_n(t_n)]\ge0,\\
\kappa_{m,n}(\beta)
&=E_N\log\cosh t_N-E_m\log\cosh t_m-E_n\log\cosh t_n.
\end{aligned}                                           \tag{FI.40}
```

The purely scalar normalization error is uniformly bounded:

```math
\left|\kappa_{m,n}(\beta)-{\beta^2\over4}\right|
\le{\beta^4\over8}.                                    \tag{FI.41}
```

Indeed `0<=u^2/2-log cosh(u)<=u^4/12`, and each of the three edge
normalization remainders has magnitude at most `beta^4/24`.

Thus the only possibly linear term in the exact radial/random-union
calculation is the nonnegative normalized-pressure gap `G`.  Convexity and
discrete stationarity do not reverse it.  In fact, if `A` minimizes at raw
temperature `t`, flipping one edge gives the exact first-order inequality

```math
a_e\langle\tau x_ix_j\rangle_A\le\tanh t,              \tag{FI.42}
```

so the derivative of
`log Zbar_A(t)-E_k log cosh(t)` is nonpositive at the optimizer.  This is
the infinitesimal version of (FI.37), again in the wrong direction for
bounding `G` from above by `o(N)`.

## 7. Global disorder softening and the incompatible parameter ranges

One could instead soften the minimum over all order-`k` signings:

```math
\mathcal F_{k,\lambda}
=-{1\over\lambda}\log\mathbb E_A
                         e^{-\lambda\log\overline Z_A}. \tag{FI.43}
```

The exact counting sandwich is

```math
P_k\le\mathcal F_{k,\lambda}
\le P_k+{E_k\log2\over\lambda}.                        \tag{FI.44}
```

Using only this uniform information, `F-P=o(k)` requires
`lambda/k->infinity`.  The bound supplied by Lemma FI.1 grows as
`1+lambda+lambda^2` after division by `lambda`, so its all-edge
replacement error is then far larger than `k`.  This is not an
impossibility theorem under a separately proved large basin around the
minimizers; it proves that the generic entropy sandwich and fixed-tilt
universality cannot be combined to manufacture such a basin.

Positive fractional moments have the complementary defect: they remain
above the exact minimum, but actual child minimality supplies no upper
calibration of their noisy-child endpoint.  The sign reversal at `s>1`
therefore has no immediate `P` implication.

## 8. The `s=1+epsilon_N` loophole

The sign in (FI.27) reverses for `s>1`, so one must separately audit a
sequence

```math
s_N=1+\epsilon_N,
\qquad \epsilon_N\downarrow0.                           \tag{FI.45}
```

This does not fall under the `s<1` sign no-go.  It also does not produce a
new power saving from covariance interpolation alone.  The endpoint term
which remains is exactly the normalized radial gap from Section 6.

To state this precisely, take the actual children minimizing at the
contracted raw temperature `t_N=beta/sqrt(N)`.  This choice is intentional:
these are the children inserted unchanged into the parent bridge, whereas
the comparison target `P_m(beta)+P_n(beta)` uses the hotter own-child raw
temperatures `t_m,t_n`.  For `k=m,n`, put
`t_k=beta/sqrt(k)` and give every internal child edge an independent
centered Gaussian of variance

```math
v_k=t_k^2-t_N^2.                                        \tag{FI.46}
```

Let `I_1` be the annealed (`s=1`) generalized pressure of the two
independent augmented-cosh child systems.  Gaussian averaging is exact,
and therefore

```math
I_1=F_m(t_N)+F_n(t_N)+{v_mE_m+v_nE_n\over2}.            \tag{FI.47}
```

Subtracting the physical child target gives

```math
\boxed{
I_1-P_m(\beta)-P_n(\beta)
=G_{m,n}(\beta)+\delta_{m,n}(\beta),}                  \tag{FI.48}
```

where `G_(m,n)` is the nonnegative radial gap in (FI.40), and

```math
\delta_{m,n}
=\sum_{k\in\{m,n\}}E_k\left[
 {t_k^2-t_N^2\over2}
 -\{\log\cosh t_k-\log\cosh t_N\}\right].             \tag{FI.49}
```

Thus `delta_(m,n)>=0`; moreover it is `O_beta(1)` uniformly over all
splits with both children nontrivial.  For example,
`0<=u^2/2-log cosh u<=u^4/12` gives a bound depending only on `beta`.
In particular, the leading endpoint cost in (FI.48) is independent of
`epsilon_N`: it is `G_(m,n)`, not an `epsilon_N N` term.

The change from `s=1` to `s=1+epsilon` is only `O(epsilon N)`.  Indeed, if
`Y` is any of the Gaussian log partitions here, its Euclidean Gaussian
Lipschitz variance proxy is the state-independent self variance
`V=O_beta(N)`.  For

```math
\phi(s)=s^{-1}\log\mathbb E e^{sY},
```

the entropy identity and Gaussian log-Sobolev inequality give

```math
0\le\phi'(s)={D(q_s\Vert\gamma)\over s^2}\le {V\over2}, \tag{FI.50}
```

where `dq_s/dgamma` is proportional to `e^(sY)`.  Hence

```math
|\phi(1+\epsilon)-\phi(1)|\le C_\beta\epsilon N.        \tag{FI.51}
```

The reversed covariance-square gain is of the same maximal order and no
larger.  Since `|q_x-q_y|<=2`, (FI.27) gives

```math
0\le \widetilde\Psi_{1+\epsilon}^{\rm int}
       -\Psi_{1+\epsilon}^{\rm cross}
\le {\epsilon\beta^2mn\over N}
\le {\epsilon\beta^2N\over4}.                          \tag{FI.52}
```

There is no hidden joint-sector loss at this near-one exponent.  For
`a,b>0` and `s=1+epsilon`, convexity and
`a^s+b^s<=(a+b)^s` give

```math
\left({a+b\over2}\right)^s
\le {a^s+b^s\over2}
\le2^{s-1}\left({a+b\over2}\right)^s.                 \tag{FI.53}
```

Thus resolving the two cosh sectors changes the generalized pressure by
at most `(s-1)log2/s=O(epsilon)`.  The same estimate applied to both
children shows that (FI.48)--(FI.52) consistently audit the physical cosh
endpoint as well as the sector-resolved covariance interpolation.

Consequently, if `G_(m,n)>=cN` on a comparable subsequence, every
`epsilon_N=o(1)` leaves the same linear endpoint excess, up to `o(N)`.
If instead `G_(m,n)=o(N)` with a summable uniform rate, then the scalar
radial/random-union inequality (FI.39) has already supplied the desired
recurrence.  The near-one positive-moment route is therefore conditional
on exactly that unresolved radial estimate; it neither proves nor
universally disproves it.  This is the precise residual loophole, and it is
not included in the `s<1` covariance sign no-go.

## 9. Scope and verdict

The proved statements have the following precise reach.

1. **Actual bridge replacement:** (FI.16)--(FI.17) applies directly to the
   cosh pressure of any actual optimizing children.  Its error is genuinely
   sublinear and uniform over every split.
2. **Entire covariance-only comparison for the natural endpoints:**
   (FI.25) is an inequality between the endpoint Gaussian generalized
   pressures, not merely the sign of one chosen parametrization.  Hence a
   different interpolation schedule using only the same two covariance
   kernels cannot reverse it.
3. **Joint-sector qualification:** resolving `cosh` into an ordinary
   sector gives Theorem FI.2.  If one instead keeps the auxiliary sector
   spin inside the Gaussian process, the covariance difference contains
   `tau^1 tau^2(q_x-q_y)^2` and has no pointwise sign.  A proof exploiting a
   special sign after the *actual joint-sector Gibbs average* would use
   structural information beyond covariance ordering and is not ruled out.
   The gap is real rather than notational.  For the one-spin-per-side base
   measure proportional to

   ```math
   \exp\{J\tau xy\},
   ```

   one has `E tau=0`, `E(tau xy)=tanh J`, and hence for two replicas

   ```math
   \mathbb E\bigl[\tau^1\tau^2(q_x-q_y)^2\bigr]
   =-2\tanh^2J<0\qquad(J\ne0).                       \tag{FI.54}
   ```

   Thus the physical joint-sector derivative can reverse the
   sector-resolved sign.  Moreover, for `s<=0` replacing a sum of two
   sectors by separately paid sector moments has no dimension-free loss:
   an exponentially imbalanced pair can lose `Theta(N)`.  Theorem FI.2 is
   therefore not a no-go for the full physical covariance architecture.
4. **No recurrence claimed:** the `O(sqrt N)` universality error is at the
   desired exponent, but the Gaussian endpoint is not upper-bounded by
   `P_m+P_n`; Theorem FI.2 explains why the standard Guerra square cannot
   provide that upper bound.  The exact radial identity leaves the same
   issue as the nonnegative gap `G_(m,n)` in (FI.40).

Therefore the sector-resolved negative-temperature, subunit positive-
fractional-moment, Gaussian/Rademacher, convexity/stationarity, and scalar
radial-scaling package does not by itself yield a genuinely sublinear
composition defect.  The obstruction for `s<1` is a proved sign theorem,
not an unnamed remainder: after the only `O_beta(1)` variance
normalization, the sector-resolved Gaussian covariance transfer is the
negative square (FI.30), and every exponent `s<1` turns it into the
nonnegative derivative (FI.27).  Equation (FI.54) prevents promoting this
statement to the joint `cosh` sector without a new sign theorem.  Exponents
`s=1+o(1)` are not claimed impossible; by (FI.48)--(FI.52) they reduce to
the explicit radial gap `G_(m,n)`, with no new leading cancellation from
the vanishing superunit exponent.
