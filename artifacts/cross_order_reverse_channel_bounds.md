# Reverse-channel bounds for the actual-child cross-order defect

Status: **proved finite-order inequalities for actual pressure-minimizing
children**, plus two sharply scoped channel-simulation obstructions.  The
detector improves the coefficient of the linear defect at sufficiently high
temperature but does not change its exponent.

Fix `N=m+n`, `L=mn`, `t=beta/sqrt(N)`, and

```math
s_m={\beta\over\sqrt m},\qquad s_n={\beta\over\sqrt n},
\qquad \rho=\tanh t.
```

Let `A,D` be exact minimizers for `P_m(beta),P_n(beta)`, and use the
notation of
[`cross_order_centered_channel_identity.md`](cross_order_centered_channel_identity.md).
In particular,

```math
E_{m,n}(\beta)\le G_{m,n}
=L\ell(t)-\Delta_A-\Delta_D-D_{\rm KL}(U\Vert\Pi).     \tag{1}
```

Here `Pi` is the actual joint orientation/bridge output law, not a surrogate.

## 1. A universal cut-code detector

The bridge marginal of `Pi` is a mixture of product binary channels

```math
W_t(B\mid q)={\exp\{t\langle B,q\rangle\}
 \over 2^L(\cosh t)^L},
\qquad q=uv^{\mathsf T},                               \tag{2}
```

over at most `2^(N-1)` distinct rank-one sign words.  Define

```math
W(B)=\max_{u,v}\langle B,uv^{\mathsf T}\rangle
    =\max_{u,v}|u^{\mathsf T}Bv|.                      \tag{3}
```

For `0<kappa<1`, put

```math
a_\kappa=2^{N-1}\exp\{-\kappa^2\rho^2L/2\},
\qquad
b_\kappa=\exp\{-(1-\kappa)^2\rho^2L/2\}.             \tag{4}
```

If `a_kappa<1-b_kappa`, then

```math
\boxed{
D_{\rm KL}(U\Vert\Pi)
\ge d_{\rm bin}(1-a_\kappa\Vert b_\kappa).}           \tag{5}
```

Indeed, for `F_kappa={W(B)<kappa rho L}`, the Rademacher subgaussian
bound and a union bound over (3) give `U(F_kappa^c)<=a_kappa`.  Under every
component (2), the planted score is a sum of `L` independent signs of mean
`rho`, so Hoeffding gives `Pi(F_kappa)<=b_kappa`.  Binary data processing
proves (5); passing from the joint law to its bridge marginal only weakens
relative entropy.

It is useful to retain the exact finite-order bound

```math
\mathfrak D_{m,n}(\beta)=
\sup_{\substack{0<\kappa<1\\a_\kappa<1-b_\kappa}}
d_{\rm bin}(1-a_\kappa\Vert b_\kappa).                 \tag{6}
```

For the explicit choice

```math
\kappa_N^*=\sqrt{
 {2((N-1)\log2+\log N)\over\rho^2L}},                 \tag{7}
```

whenever `kappa_N^*<1`, one has `a_(kappa_N^*)=1/N` and

```math
\mathfrak D_{m,n}(\beta)
\ge(1-N^{-1}){(1-\kappa_N^*)^2\rho^2L\over2}
   -h(N^{-1}).                                        \tag{8}
```

## 2. An optimizer-specific radial payment

Let `M(C)=max_x|H_C(x)|` and `M_k=min_C M(C)`.  The two configurations
`x,-x` which attain `M(A)` contribute the same quadratic energy.  Hence

```math
\phi_A(s_m)\ge(1-m)\log2+\log\cosh(s_mM(A)),
\qquad
\phi_A(t)\le\log\cosh(tM(A)).                         \tag{9}
```

The difference of the two log-cosines is increasing in `M(A)>=M_m`, so

```math
\boxed{
\Delta_A\ge\underline\Delta_m:=
\left[(1-m)\log2+
 \log{\cosh(s_mM_m)\over\cosh(tM_m)}\right]_+}        \tag{10}
```

and, more simply,

```math
\Delta_A\ge[(s_m-t)M_m-m\log2]_+.                    \tag{11}
```

The same statements hold for `D`.  Combining (1), (5), and (10) gives the
fully explicit actual-child recurrence

```math
\boxed{
E_{m,n}(\beta)
\le L\ell(t)-\underline\Delta_m-\underline\Delta_n
       -\mathfrak D_{m,n}(\beta).}                    \tag{12}
```

This improves the annealed coefficient above a calculable temperature but
remains linear at comparable splits.  More importantly, if

```math
\mathcal B_{m,n}
=K_m[\ell(s_m)-\ell(t)]+K_n[\ell(s_n)-\ell(t)],
```

then the following is an exact quantitative sufficient lemma:

```math
\boxed{
\underline\Delta_m+\underline\Delta_n
\ge\mathcal B_{m,n}-\mathfrak D_{m,n}(\beta)-\omega_N
\quad\Longrightarrow\quad
E_{m,n}(\beta)
\le\omega_N+{\beta^2\over4}+{\beta^4\over12}.}       \tag{13}
```

Thus a power-saving remainder in the displayed, finite-order inequality
would give the same power-saving cross-order defect without estimating the
three centered terms separately.

### Comparable-split coefficient

If `m/N->theta in (0,1)` and `a=theta(1-theta)`, the detector activates
when

```math
\beta^2a>2\log2.                                      \tag{14}
```

Optimizing `kappa` in (5) gives

```math
{\mathfrak D_{m,n}(\beta)\over N}
\ge {\beta^2a\over2}
 \left(1-\sqrt{2\log2\over\beta^2a}\right)^2-o(1).  \tag{15}
```

The annealed bridge coefficient left after this detector is therefore

```math
R_{\rm det}(\beta,\theta)
=\beta\sqrt{2\theta(1-\theta)\log2}-\log2.            \tag{16}
```

Using the rigorous lower bound
`M_k>=(c_*-o(1))k^(3/2)`, `c_*=0.336493364431...`, in
(11)--(12) yields

```math
\begin{aligned}
\limsup {E_{m,n}(\beta)\over N}\le{}&
R_{\rm det}(\beta,\theta)\\
&-\theta[\beta c_*(1-\sqrt\theta)-\log2]_+\\
&-(1-\theta)
 [\beta c_*(1-\sqrt{1-\theta})-\log2]_+ .             \tag{17}
\end{aligned}
```

At the equal split this becomes

```math
\limsup {E_{N/2,N/2}(\beta)\over N}
\le0.5887050113\,\beta-\log2
 -[0.0985566246\,\beta-\log2]_+.                     \tag{18}
```

The detector changes the high-temperature annealed coefficient from
quadratic in `beta` to linear in `beta`, but not from exponent one to a
sublinear power of `N`.  Even replacing `c_*` by the conjectural optimal
constant `1/2` leaves a positive linear coefficient, so this detector plus
the elementary cap payment cannot close (13).

## 3. The exact reverse-posterior target

Let `pi(epsilon)` be the output orientation law, `mu_epsilon` the
conditional prior on the latent rank-one bridge word, and
`Pi_epsilon=mu_epsilon W_t`.  If `mu_(epsilon,B)` is the posterior of the
word after observing a **uniform** bridge `B`, the reverse-KL chain rule
gives, separately for each orientation,

```math
\boxed{
L\ell(t)=D(U_B\Vert\Pi_\epsilon)
 +\mathbb E_{B\sim U_B}
   D(\mu_\epsilon\Vert\mu_{\epsilon,B}).}             \tag{19}
```

Set

```math
\Lambda_{\rm rev}={1\over2}\sum_{\epsilon=\pm1}
 \mathbb E_{B\sim U_B}
 D(\mu_\epsilon\Vert\mu_{\epsilon,B}).               \tag{20}
```

Then the joint orientation chain rule and (1) give the exact implication

```math
\boxed{
E_{m,n}(\beta)
\le\Lambda_{\rm rev}-\Delta_A-\Delta_D
 -D(U_\epsilon\Vert\pi),}                             \tag{21}
```

and therefore

```math
\boxed{
\Lambda_{\rm rev}
\le\Delta_A+\Delta_D+D(U_\epsilon\Vert\pi)+\omega_N
\quad\Longrightarrow\quad E_{m,n}(\beta)\le\omega_N.} \tag{22}
```

This identity fixes a misleading planted-recovery intuition.  Standard
factor recovery studies a posterior averaged under `B~Pi_epsilon`, usually
in the forward KL direction.  The defect contains the reverse divergence
averaged under `B~U_B`.  Neither its measure nor its KL direction is the
standard planted one.

## 4. Exact reverse Blackwell simulation still leaves optimizer excess

Let `p_(m,u)` be the child cut-channel output law and let `T` be the BSC
degradation from `s_m` to `t`.  At the selected minimizing output `A`, set

```math
\kappa_A(z)=T(A\mid z),
\qquad
\widehat\kappa_A(z)=
 {\kappa_A(z)p_{m,s_m}(z)\over p_{m,t}(A)}.            \tag{23}
```

Direct expansion gives

```math
\boxed{
\Gamma_A=D(\kappa_A\Vert\widehat\kappa_A)+H_A,
\qquad
H_A=\mathbb E_{\kappa_A}
 \log{p_{m,s_m}(z)\over p_{m,s_m}(A)}\ge0.}           \tag{24}
```

Suppose a kernel `K` from bridge/orientation outputs satisfies

```math
K_\#U=\kappa_A\otimes\kappa_D,
\qquad
K_\#\Pi\le e^{\varepsilon_N}
 (\widehat\kappa_A\otimes\widehat\kappa_D)
```

pointwise.  Data processing and (24) yield

```math
D(U\Vert\Pi)
\ge\Gamma_A+\Gamma_D-H_A-H_D-\varepsilon_N.
```

Consequently the promised direct defect implication is

```math
\boxed{
H_A+H_D+\varepsilon_N\le\omega_N
\quad\Longrightarrow\quad
E_{m,n}(\beta)
\le\omega_N+{\beta^2\over4}+{\beta^4\over12}.}       \tag{25}
```

Even an exact reverse Blackwell simulation leaves the optimizer-excess
terms.  Ordinary forward Le Cam comparison is weaker and cannot control
the selected least-likely output.

## 5. A literal edgewise triangle simulator has a linear deficiency

For a fixed left-child edge `ij`, every disjoint bridge triangle gives

```math
B_{ia}B_{ja}=x_ix_j\eta_{ia}\eta_{ja}.                \tag{26}
```

Even granting the unknown common orientation as an oracle, this is a BSC
observation of the desired augmented child edge with correlation `rho^2`.
Together with the direct weak child bit of correlation `rho`, `r` such
independent observations have half-log-likelihood ratio

```math
tw+\operatorname{arctanh}(\rho^2)\sum_{\ell=1}^r c_\ell.
```

If

```math
r\operatorname{arctanh}(\rho^2)<t,                    \tag{27}
```

its sign is always the original weak bit `w`.  The binary experiment's
total variation is still exactly `rho`; it cannot simulate any BSC with
larger correlation.  Hence each target edge needs at least

```math
r_*=left\lceil {t\over\operatorname{arctanh}(\rho^2)}
       \right\rceil
={\sqrt N\over\beta}+O_\beta(1)                       \tag{28}
```

source-disjoint triangle samples merely to improve the weak bit.

At the balanced split there are only `N^2/4` source bridge bits, while
each triangle uses two.  Thus only `O_beta(N^(3/2))` of the
`N^2/4+O(N)` child target edges can be amplified.  Comparing the remaining
binary experiments by data processing and Pinsker gives the minimax
conditional-KL deficiency lower bound

```math
\boxed{
{\beta^2(\sqrt2-1)^2\over8}N-O_\beta(\sqrt N).}       \tag{29}
```

This is deliberately narrow: it excludes source-disjoint, edgewise
triangle tensorization, even with an orientation oracle.  It does not
exclude global reuse of bridge bits.  Such reuse still has to confront the
reverse-uniform posterior in (19)--(22).

## Quantitative conclusion

Equations (12) and (17) are genuine improved all-order bounds for the
actual selected children, but their exponent remains one.  Equations
(13), (22), and (25) state three exact ways a future estimate would imply a
sublinear `E_(m,n)`.  The detector, elementary optimizer payment, ordinary
factor recovery, and source-disjoint triangle amplification cannot provide
that estimate.  A successful channel argument must use globally reused
complete-graph cut constraints in the reverse-uniform direction.
