# Low conference pressure does not regularize the bridge gradient

**Status.**  Task-local theorem and audit draft.  This file does not modify
canonical sources.  It answers the question whether the lower-pressure
sublevel itself can replace the operator-norm conditioning in
`conference_regular_conditioned_all_tilts.md`.

The answer is negative in the strongest form relevant to the proposed
Talagrand argument.  A strict lower-pressure sublevel contains two genuine
sign bridges whose pressure difference is `Theta(r)` but whose Frobenius
distance is only `Theta(r^(3/4))`.  At one of them the pressure gradient is
therefore `Omega(r^(1/4))`.  Consequently neither the raw pressure nor any
extension agreeing with it on that sign sublevel can have a
dimension-free Frobenius Lipschitz constant.

There is also a positive localization result.  The operator-regular theorem
is valid after two small repairs, and its speed-`r^2` tail extends to every
low-pressure bridge within `Theta(r^(3/2))` Hamming distance of the regular
sign set.  Thus any hypothetical speed-`r` favorable basin must live deep
inside the operator-irregular sector.

## 1. Exact differential, covariance, and entropy identities

Let `A=A_r` be a symmetric conference signing, fix an orientation
`epsilon in {+-1}`, and put

```math
S_{epsilon,B}=
\begin{pmatrix}A&B\\B^T&\epsilon A\end{pmatrix},
\qquad
t={\beta\over\sqrt{2r}}.
\tag{SG.1}
```

For a real `r` by `r` bridge `B`, write

```math
f_\epsilon(B)
=\log\left[2^{-2r}\sum_{x,y}
 \cosh\left(tG_{\epsilon,B}(x,y)\right)\right],
\tag{SG.2}
```

where

```math
G_{\epsilon,B}(x,y)
=H_A(x)+\epsilon H_A(y)+x^TBy,
\qquad H_A(x)={1\over2}x^TAx.
\tag{SG.3}
```

Introduce an auxiliary `sigma in {+-1}` and let `mu_B` be the Gibbs law on
`(sigma,x,y)` with density proportional to
`exp(t sigma G_(epsilon,B)(x,y))` relative to the uniform law `nu`.  Define

```math
M_B=\mathbb E_{\mu_B}[\sigma xy^T].
\tag{SG.4}
```

Direct differentiation gives the exact identities

```math
\boxed{\nabla_B f_\epsilon(B)=tM_B,}
\tag{SG.5}
```

```math
\boxed{
\nabla^2 f_\epsilon(B)_{ij,k\ell}
=t^2\operatorname {Cov}_{\mu_B}
 (\sigma x_i y_j,\sigma x_k y_\ell),}
\tag{SG.6}
```

and hence

```math
\operatorname {Tr}\nabla^2 f_\epsilon(B)
=t^2\{r^2-\|M_B\|_F^2\}.
\tag{SG.7}
```

In replica form,

```math
\|\nabla_Bf_\epsilon(B)\|_F^2
=t^2\mathbb E_{\mu_B^{\otimes2}}
 [\sigma\sigma'(x^Tx')(y^Ty')].
\tag{SG.8}
```

The Gibbs entropy is

```math
D(\mu_B\Vert\nu)
=t\mathbb E_{\mu_B}[\sigma G_{\epsilon,B}]-f_\epsilon(B),
\tag{SG.9}
```

while the bridge part of its energy is exactly

```math
B:\nabla f_\epsilon(B)
=t\mathbb E_{\mu_B}[\sigma x^TBy].
\tag{SG.10}
```

These formulas distinguish pressure, response, and covariance.  In
particular, a small value of `f` does not by itself upper-bound the entropy
or the bridge response.

There is a useful universal entropy inequality.  For every real matrix `C`,

```math
\log\mathbb E_\nu e^{\sigma x^TCy}
\le {r\over2}\|C\|_{op}^2.
\tag{SG.11}
```

Indeed, condition on `(sigma,x)`, average the independent `y_j`, use
`log cosh u<=u^2/2`, and then
`||C^Tx||_2^2<=r||C||_(op)^2`.  Entropy duality with `C=M_B/r` therefore
gives

```math
\boxed{
D(\mu_B\Vert\nu)\ge {\|M_B\|_F^2\over2r},
\qquad
\|\nabla f_\epsilon(B)\|_F^2
\le\beta^2D(\mu_B\Vert\nu).}
\tag{SG.12}
```

Since `D(mu_B||nu)<=(2r+1)log 2`, this recovers only the global
`O_beta(sqrt(r))` gradient scale.  Crucially, (SG.9) supplies no better
bound from the scalar condition `f(B)<=ar`.

## 2. A sign-bridge counterexample inside a strict low-pressure sublevel

Use the Paley conference sequence in the range

```math
0<\beta<{\sqrt2\over6}.
\tag{SG.13}
```

The archived conference calculation gives

```math
h_\beta=2\psi(\beta/\sqrt2)+{\beta^2\over4},
\qquad
\tau_\beta=2\psi(\beta),
\qquad
\gamma(\beta)=h_\beta-\tau_\beta>0.
\tag{SG.14}
```

For orientation `epsilon=-1`, the universal-double bridge

```math
B_r^0=A_r+I
\tag{SG.15}
```

is a sign matrix and the archived universal-double pressure identity says

```math
{1\over r}f_-(B_r^0)\longrightarrow\tau_\beta.
\tag{SG.16}
```

For completeness, the input is independently visible from

```math
(S_r^0)^2
=\operatorname {diag}
 \big((2r-1)I+2A_r,(2r-1)I+2A_r\big).
\tag{SG.16a}
```

Thus the empirical law of `tS_r^0` tends to the symmetric Bernoulli law on
`{+-beta}`.  Every fixed power reduces, blockwise, to a linear combination
of `I` and `A_r`, so the required diagonal/off-diagonal delocalization also
holds.  The archived strict-high-temperature pressure theorem then gives
`f_-(B_r^0)/(2r)->psi(beta)`, which is precisely (SG.16).

Also, for

```math
S_r^0=\begin{pmatrix}A_r&B_r^0\\(B_r^0)^T&-A_r\end{pmatrix},
\tag{SG.17}
```

the conference identity implies the elementary bound

```math
\|S_r^0\|_{op}
\le2\sqrt{r-1}+1.
\tag{SG.18}
```

### Theorem SG.1 (macroscopic planted increment with a mesoscopic edit)

Fix `delta` with `0<delta<gamma(beta)` and set

```math
c=\sqrt{{\sqrt2\,\delta\over\beta}},
\qquad
k=\lfloor c r^{3/4}\rfloor.
\tag{SG.19}
```

Choose any `k` rows `I` and `k` columns `J`.  Form the sign bridge `B_r^1`
by overwriting the `I` by `J` block of `B_r^0` with all `+1` entries and
leaving every other entry unchanged.  Then

```math
\boxed{
f_-(B_r^1)-f_-(B_r^0)=\delta r+o(r),}
\tag{SG.20}
```

```math
\boxed{
\|B_r^1-B_r^0\|_F=(\sqrt2+o(1))c r^{3/4},}
\tag{SG.21}
```

and

```math
\boxed{
\|\nabla f_-(B_r^1)\|_F
\ge {\delta\over\sqrt2c}r^{1/4}+o(r^{1/4}).}
\tag{SG.22}
```

Consequently, for every fixed

```math
\tau_\beta+\delta<a<h_\beta,
\tag{SG.23}
```

both `B_r^0` and `B_r^1` lie in the sign-bridge sublevel
`{B:f_-(B)<=ar}` for all large `r`, but every Frobenius-Lipschitz function
agreeing with `f_-` on that sublevel has Lipschitz constant
`Omega(r^(1/4))`.

**Proof.**  Put `D_r=B_r^1-B_r^0` and let `1_I,1_J` denote the all-one
vectors on the overwritten block.  Since
`||B_r^0||_(op)<=sqrt(r-1)+1`,

```math
\left|1_I^TB_r^0[I,J]1_J\right|
\le k(\sqrt{r-1}+1)=o(k^2).
\tag{SG.24}
```

Every nonzero entry of `D_r` equals `2`, and therefore

```math
d_r:=1_I^TD_r1_J
=\sum_{ij}|(D_r)_{ij}|
=k^2-o(k^2),
\qquad
\|D_r\|_F^2=2d_r.
\tag{SG.25}
```

In particular `sup_(x,y)|x^TD_ry|=d_r`, so pointwise comparison of the
Boltzmann weights gives

```math
f_-(B_r^1)\le f_-(B_r^0)+td_r.
\tag{SG.26}
```

For the reverse inequality use the auxiliary-`sigma` partition function.
Map every `(sigma,x,y)` to the state with the same outside coordinates and
with

```math
x_I=1_I,
\qquad y_J=\sigma1_J.
\tag{SG.27}
```

Each pinned state has exactly `2^(2k)` preimages.  If `z=(x,y)` and `z'`
is its pinned image, then (SG.18) gives

```math
t|H_{S_r^0}(z')-H_{S_r^0}(z)|
\le4t\|S_r^0\|_{op}\sqrt{kr}
=O_\beta(\sqrt{rk})=o(r).
\tag{SG.28}
```

On every pinned state,
`sigma x^TD_ry=d_r`.  Summing the pinned weights and comparing them with
all baseline weights consequently yields

```math
f_-(B_r^1)
\ge f_-(B_r^0)+td_r
 -4t\|S_r^0\|_{op}\sqrt{kr}-2k\log2.
\tag{SG.29}
```

Both errors are `o(r)`, while

```math
{td_r\over r}\longrightarrow{\beta c^2\over\sqrt2}=\delta.
\tag{SG.30}
```

This proves (SG.20)--(SG.21).  Convexity at the endpoint gives

```math
\nabla f_-(B_r^1):D_r
\ge f_-(B_r^1)-f_-(B_r^0).
\tag{SG.31}
```

Cauchy--Schwarz and (SG.20)--(SG.21) prove (SG.22).  Equations
(SG.14), (SG.16), and (SG.20) prove the sublevel assertion.  Finally, any
extension agreeing at the two displayed bridges has Lipschitz constant at
least their pressure difference divided by their Frobenius distance, the
same `Omega(r^(1/4))` lower bound. `square`

### Consequence for the proposed proof method

The real sublevel of the convex function `f_-` is convex and contains these
two sign points.  Theorem SG.1 therefore rules out both of the following:

1. a dimension-free bound on `||nabla f||_F` throughout a strict
   conference lower-pressure sublevel;
2. a dimension-free convex Frobenius-Lipschitz extension agreeing with the
   pressure on all sign points of that sublevel.

Thus the operator conditioning in the regular-sector proof cannot simply
be replaced by the scalar condition `f<=ar`.  The theorem does **not**
disprove a speed-faster-than-`r` lower tail by another method.  Its
`r^(1/4)` obstruction would turn the naive convex-Lipschitz exponent at a
linear deviation into order `r^(3/2)`, which would still be
superexponential at speed `r` if a matching global upper bound were proved.
No such upper bound follows from the scalar sublevel alone here.

## 3. Audit and repair of the operator-regular theorem

Consider the draft `conference_regular_conditioned_all_tilts.md`.  Its
fixed-tilt theorem is valid.  The claimed extension to every positive
sequence `lambda_r=o(sqrt(r))` is also valid, but the written proof needs a
separate argument when `lambda_r` is exponentially small.

Fix `delta_0>0` and `kappa<1/2` with

```math
{\beta(3+\delta_0)\over\sqrt2}<\kappa,
\tag{SG.32}
```

and set

```math
\mathcal K_{\epsilon,r}
=\{B:\|tS_{\epsilon,B}\|_{op}\le\kappa\}.
\tag{SG.33}
```

The formatting error in (RC.5) should read

```math
\|tS_{\epsilon,B}\|_{op}
\le{\beta\over\sqrt2}
\left(\sqrt{1-1/r}+{\|B\|_{op}\over\sqrt r}\right).
\tag{SG.34}
```

On the convex real set (SG.33), Theorem 1.3 of
`high_temperature_frobenius_pressure_stability.md` gives, for a bridge
increment `E`,

```math
|df_B[E]|
\le {K_\kappa\over2}
 \left\|t\begin{pmatrix}0&E\\E^T&0\end{pmatrix}\right\|_*
=K_\kappa t\|E\|_*
\le {K_\kappa\beta\over\sqrt2}\|E\|_F.
\tag{SG.35}
```

Here the block matrix has nuclear norm `2||E||_*`; this checks the factor
in the dimension-free Frobenius constant.  At a boundary point of `K`, the
same estimate for the full gradient follows directly from the archived
covariance bound: each of the two Ising components has covariance operator
norm at most `K_kappa`, its cross block has Frobenius norm at most
`K_kappa sqrt(r)`, and the `cosh` Gibbs law is their convex mixture.  Thus
there is no hidden one-sided-tangent assumption.  The supporting-plane supremum
over `K_(epsilon,r)` is consequently a convex
`K_kappa beta/sqrt(2)`-Lipschitz extension `g_epsilon` which agrees with
`f_epsilon` on that set.

The rectangular Rademacher norm tail gives

```math
q_r:=\Pr\{B\notin\mathcal K_{\epsilon,r}\}\le e^{-c r}.
\tag{SG.36}
```

Convex concentration yields

```math
1\le\mathbb E e^{-\lambda(g-\mathbb Eg)}
\le e^{C\lambda^2}.
\tag{SG.37}
```

By Cauchy--Schwarz the omitted centered moment is at most

```math
\mathbb E[1_{\mathcal K^c}e^{-\lambda(g-\mathbb Eg)}]
\le\exp\{-cr/2+2C\lambda^2\}.
\tag{SG.38}
```

For `lambda_r=o(sqrt(r))`, this is `e^(-Omega(r))`.  If
`lambda_r>=e^{-c_1r}`, (SG.37)--(SG.38), divided by `lambda_r r`, give the
claimed rate after choosing `c_1` smaller than the exponent in (SG.38).
If `0<lambda_r<e^{-c_1r}`, use instead `0<=f<=kappa r` on `K` and the
bounded-range exponential lemma under the conditioned law:

```math
0\le\mathbb E[f\mid\mathcal K]
 +{1\over\lambda_r}\log
   \mathbb E[e^{-\lambda_r f}\mid\mathcal K]
\le {\lambda_r\kappa^2r^2\over8}=o(r).
\tag{SG.39}
```

Thus, including the two orientations whose means have the same rate,

```math
\boxed{
{\mathcal R^K_{\lambda_r,r}\over r}\longrightarrow h_\beta
\quad\hbox{for every positive }\lambda_r=o(\sqrt r).}
\tag{SG.40}
```

This validates the draft's range.  The proof should use the split above,
rather than claim that an unspecified `o(1)` normalization error remains
harmless after division by an arbitrarily small `lambda_r`.

## 4. The speed-`r^2` tail extends through a mesoscopic edit collar

The regular-sector result has a useful strengthening which does not require
gradient control outside `K`.

### Theorem SG.2 (Hamming-collar localization)

Fix `eta>0` and one orientation.  Let

```math
\mathcal L_{\eta,r}
=\{B\in\{+-1\}^{r\times r}:
 f_\epsilon(B)\le(h_\beta-\eta)r\},
\tag{SG.41}
```

and let `K^sign_(epsilon,r)` be the sign bridges in (SG.33).  Set

```math
s_r=\left\lfloor{\eta r\over4t}\right\rfloor
=\Theta_{\beta,\eta}(r^{3/2}).
\tag{SG.42}
```

Then there is `c_(beta,kappa,eta)>0` such that

```math
\boxed{
\Pr\{B\in\mathcal L_{\eta,r},
 d_H(B,K^{sign}_{\epsilon,r})\le s_r\}
\le e^{-c_{\beta,\kappa,\eta}r^2}.}
\tag{SG.43}
```

**Proof.**  Suppose `B` is counted on the left and choose a regular sign
bridge `B_0` with `d_H(B,B_0)<=s_r`.  A bridge-bit flip changes `f` by at
most `2t`, so

```math
f_\epsilon(B_0)
\le f_\epsilon(B)+2ts_r
\le(h_\beta-\eta/2)r.
\tag{SG.44}
```

The regular-sector speed-`r^2` theorem, with `eta/2`, shows that the number
of such centers is at most

```math
2^{r^2}e^{-c_0r^2}.
\tag{SG.45}
```

The Hamming ball around each center has volume at most

```math
\sum_{j\le s_r}{r^2\choose j}
\le\left({e r^2\over s_r}\right)^{s_r}
=\exp\{O_{\beta,\eta}(r^{3/2}\log r)\}
=e^{o(r^2)}.
\tag{SG.46}
```

The union bound proves (SG.43), after decreasing `c_0`. `square`

In particular, if the full low-pressure set has probability `e^(-O(r))`,
then a subset of the same speed-`r` mass lies at Hamming distance greater
than `Theta(r^(3/2))` from every operator-regular sign bridge.  Theorem
SG.1 operates exactly at this mesoscopic edit exponent; it is consistent
with, and helps explain, the collar scale.

## 5. Research judgment

The scalar-sublevel replacement is closed:

```text
low pressure
  does not imply dimension-free Frobenius response,
  does not admit a dimension-free agreeing extension,
  and therefore cannot directly trigger the r^2 Talagrand proof.
```

The operator-regular theorem survives audit, including all positive
`lambda_r=o(sqrt(r))`, and the Hamming-collar theorem sharply localizes the
remaining lower-LDP problem.  A finite-tilt phase, or any speed-`r`
target-reaching basin, must be supported on bridges which are deeply
operator-irregular in edit distance.  The next legitimate theorem is
therefore a joint operator-irregularity/pressure estimate; another scalar
sublevel-gradient argument should not be pursued.
