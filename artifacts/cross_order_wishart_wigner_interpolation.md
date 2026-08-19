# Covariance-matched Wishart--Wigner interpolation

Status: **proved all-order identity up to an explicit `O_beta(sqrt(N))`
universality error for actual pressure-minimizing children**.  It retains the
same-switch cancellation which is lost by separate scalar channels.  The
result gives a direct power-saving defect implication, but minimizer
optimality controls its new endpoint regrets in the wrong direction.

Write

```math
\phi_A(u)=\log\mathbb E_x\cosh(uH_A(x)),
\quad N=m+n,
\quad t={\beta\over\sqrt N},
\quad s_m={\beta\over\sqrt m},
\quad s_n={\beta\over\sqrt n}.
```

Let `A,D` be exact minimizers defining `P_m(beta),P_n(beta)`.  Introduce the
contracted child measures

```math
\Delta_A=\phi_A(s_m)-\phi_A(t),
\qquad
\Delta_D=\phi_D(s_n)-\phi_D(t),                       \tag{0}
```

and

```math
\mu_A(\tau,x)=
 {2^{-m-1}e^{t\tau H_A(x)}\over e^{\phi_A(t)}},
\qquad
\mu_D(\sigma,y)=
 {2^{-n-1}e^{t\sigma H_D(y)}\over e^{\phi_D(t)}}.      \tag{1}
```

## 1. Exact bridge endpoint

For `epsilon in {+-1}`, put

```math
\mathcal F_B^\epsilon=
\mathbb E_{B\sim U_B}\log\!\left[
 2\mathbb E_{\mu_A\otimes\mu_D}
 \mathbf1_{\{\tau\sigma=\epsilon\}}
 e^{t\tau x^{\mathsf T}By}\right],
\qquad
\mathcal F_B={\mathcal F_B^++\mathcal F_B^-\over2}.  \tag{2}
```

The exact likelihood ratio of the joint orientation/bridge output law is

```math
{d\Pi\over dU}(\epsilon,B)
=2(\cosh t)^{-mn}
 \mathbb E_{\mu_A\mu_D}
 \mathbf1_{\{\tau\sigma=\epsilon\}}
 e^{t\tau x^{\mathsf T}By}.                           \tag{3}
```

Consequently

```math
D(U\Vert\Pi)=mn\ell(t)-\mathcal F_B
```

and the exact selected-child construction defect is

```math
\boxed{
G_{m,n}=\mathcal F_B-\Delta_A-\Delta_D,
\qquad E_{m,n}(\beta)\le G_{m,n}.}                    \tag{4}
```

## 2. The random-competitor child endpoint

Let the independent signs `xi^A_ij` have mean
`t/s_m=sqrt(m/N)`, set `C^A=A odot xi^A`, and put

```math
\zeta^A_{ij}=s_mC^A_{ij}-tA_{ij},
\qquad
v_m=\mathbb E(\zeta^A_{ij})^2
=s_m^2-t^2={\beta^2n\over mN}.                        \tag{5}
```

Define `C^D,zeta^D,v_n` analogously and the exact competitor regrets

```math
\mathcal R_A=
\mathbb E_{C^A}[\phi_{C^A}(s_m)-\phi_A(s_m)]\ge0,
\qquad
\mathcal R_D\ge0.                                     \tag{6}
```

The inequalities are pointwise consequences of exact child minimality.
For a realization of the perturbations, define

```math
\begin{aligned}
Z_\epsilon={}&2\mathbb E_{\mu_A\mu_D}
 \mathbf1_{\{\tau\sigma=\epsilon\}}
 \exp\!\left\{
  \tau\sum_{i<j}\zeta^A_{ij}x_ix_j
 +\sigma\sum_{a<b}\zeta^D_{ab}y_ay_b\right\},\\
Z_A={}&\mathbb E_{\mu_A}
 e^{\tau\sum_{i<j}\zeta^A_{ij}x_ix_j},
\qquad
Z_D=\mathbb E_{\mu_D}
 e^{\sigma\sum_{a<b}\zeta^D_{ab}y_ay_b}.
\end{aligned}                                         \tag{7}
```

Since `(Z_++Z_-)/2=Z_AZ_D`, the orientation arithmetic/geometric gap

```math
\mathcal O=\mathbb E\left[
 \log Z_A+\log Z_D-{\log Z_++\log Z_-\over2}\right]
\ge0                                                   \tag{8}
```

is exact, and the discrete child endpoint equals

```math
{1\over2}\sum_\epsilon\mathbb E\log Z_\epsilon
=\Delta_A+\Delta_D+\mathcal R_A+\mathcal R_D-
 \mathcal O.                                          \tag{9}
```

## 3. A signed same-switch interpolation identity

For each orientation sector use the finite base measure
`2 1_(tau sigma=epsilon) mu_A mu_D`.  With all Gaussians independent, put

```math
\begin{aligned}
\mathcal H_B={}&t\tau\sum_{ia}g_{ia}x_iy_a,\\
\mathcal H_C={}&\tau\left[
 \sqrt{v_m}\sum_{i<j}g^A_{ij}x_ix_j
 +\epsilon\sqrt{v_n}\sum_{a<b}g^D_{ab}y_ay_b
 +{\beta\over\sqrt2}g_0\right],\\
\mathcal H_u={}&\sqrt u\,\mathcal H_B+
 \sqrt{1-u}\,\mathcal H_C.
\end{aligned}                                         \tag{10}
```

For two replicas write `X=<x^1,x^2>` and `Y=<y^1,y^2>`.  Because
`sigma=epsilon tau` inside one sector and
`sum_(i<j)x_i^1x_j^1x_i^2x_j^2=(X^2-m)/2`, the covariance gap is exactly

```math
\boxed{
C_C(1,2)-C_B(1,2)
=\tau^1\tau^2{\beta^2\over2N}
 \left(\sqrt{n/m}\,X-\sqrt{m/n}\,Y\right)^2.}        \tag{11}
```

The one-dimensional Gaussian in (10) repairs the otherwise missing
diagonal constant.  Gaussian integration by parts, whose diagonal terms
cancel, therefore gives

```math
F_\epsilon(1)-F_\epsilon(0)
={\beta^2\over4N}\int_0^1
 \mathbb E\left\langle\tau^1\tau^2
 (\sqrt{n/m}X-\sqrt{m/n}Y)^2\right\rangle_{\epsilon,u},du. \tag{12}
```

Define

```math
\begin{aligned}
\mathcal S_{\rm sgn}={}&{\beta^2\over8N}
 \sum_{\epsilon=\pm1}\int_0^1
 \mathbb E\left\langle\tau^1\tau^2
 (\sqrt{n/m}X-\sqrt{m/n}Y)^2\right\rangle du,\\
\mathcal S_+={}&{\beta^2\over8N}
 \sum_{\epsilon=\pm1}\int_0^1
 \mathbb E\left\langle\mathbf1_{\{\tau^1=\tau^2\}}
 (\sqrt{n/m}X-\sqrt{m/n}Y)^2\right\rangle du.
\end{aligned}                                         \tag{13}
```

Then `S_sgn<=S_+` and `S_+<=beta^2 mn/N<=beta^2N/4`.  In particular,
opposite-orientation replicas remain as favorable cancellation in the
signed term; the channels have not been bounded separately.

## 4. Uniform discrete-to-Gaussian error

For

```math
f(z)=\log\int\exp\{\sum_ez_eh_e\}\,d\nu,
\qquad h_e\in\{\pm1\},
```

the third partial derivative is a third centered moment and has absolute
value at most eight.  Coordinatewise matched-mean/variance Lindeberg
replacement therefore costs at most

```math
{4\over3}\sum_e
 (\mathbb E|X_e|^3+\mathbb E|Y_e|^3).                 \tag{14}
```

For the bridge this is `O(beta^3 mn/N^(3/2))`.  For the two child
perturbations, use `|zeta^A|<=2s_m` and
`E(zeta^A)^2=v_m` (and analogously on the right) to obtain
`O(beta^3(sqrt(m)+sqrt(n)))`.  The scalar Gaussian changes a log partition
by at most `beta E|g_0|/sqrt(2)=beta/sqrt(pi)`.  Thus, uniformly in the
children and every split,

```math
|\operatorname{Err}_{m,n}|
\le C\beta^3\sqrt N+{\beta\over\sqrt\pi}              \tag{15}
```

for a universal finite `C`.

Combining (4), (9), (12), and (15) proves

```math
\boxed{
G_{m,n}=\mathcal R_A+\mathcal R_D-\mathcal O
 +\mathcal S_{\rm sgn}+\operatorname{Err}_{m,n}.}    \tag{16}
```

## Direct cross-order implication and boundary

For any `delta>0`, (13), (15), and (16) give the requested direct arrow

```math
\boxed{
\mathcal R_A+\mathcal R_D+\mathcal S_+
\le K_\beta N^{1-\delta}
\quad\Longrightarrow\quad
E_{m,n}(\beta)
\le K_\beta N^{1-\delta}+C\beta^3\sqrt N
 +{\beta\over\sqrt\pi}.}                             \tag{17}
```

Equivalently, the conclusion is
`O_beta(N^(1-min(delta,1/2)))`.  The sharper premise with
`R_A+R_D-O+S_sgn` in place of the nonnegative upper envelope yields the same
conclusion.  The theorem covers all orders and every split at every fixed
`beta`, however large.  If `beta=beta_N` grows, this Lindeberg error is
`o(N)` whenever `beta_N=o(N^(1/6))`.

Exact minimizer optimality supplies only `R_A,R_D>=0`, the wrong direction
for (17), and `S_sgn` has no fixed sign.  Thus the covariance-matched
Wishart--Wigner interpolation by itself does not prove a sublinear defect.

There is also no planted-factor shortcut.  Under the output law,

```math
\sum_aB_{ia}B_{ja}
=x_ix_j\sum_a\eta_{ia}\eta_{ja},                     \tag{18}
```

but the cross-order reverse divergence is averaged at `B~U`, not at the
planted output.  Moreover, `x_ix_j` is a switching gauge:
`phi_(gamma diag(x) A diag(x))=phi_A` exactly.  Recovering the raw factor
therefore contributes no child-pressure increment.  Only cycle/syndrome
information and the joint signed term can help.  This does not exclude a
more complete optimizer-specific use of the full Wishart law; it excludes
closing this covariance-matched interpolation from minimizer stationarity
alone.
