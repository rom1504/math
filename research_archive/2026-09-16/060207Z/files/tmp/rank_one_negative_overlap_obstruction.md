# Rank-one entropy forces a negative-tilt cavity-overlap phase at large beta

Date: 2026-08-18.

Status: rigorous theorem for the **actual child channel** (indeed for every
prior supported on rank-one sign matrices).  It shows that the current decay
target for `rhohat_N^-(lambda)` is false throughout an explicit large-`beta`
regime.  No conference or Paley surrogate is used.

## 1. Setup

Let `m+n=N`, `d=mn`, `t=beta/sqrt(N)`, and `rho=tanh t`.  Condition on one
relative child orientation and let `mu` be the actual zero-bridge augmented
Gibbs prior on

```math
Q=(Q_{ij})=(\tau X_iZ_j)\in\{-1,1\}^{m\times n}.
```

Every `Q` is rank one, so the number `K` of distinct matrices in the support
of `mu` satisfies

```math
K\le2^{m+n-1}=2^{N-1}.                              \tag{RO.1}
```

For fair bridge disorder `U`, the forward likelihood is

```math
p(B)=E_\mu\prod_e(1+\rho B_eQ_e),
```

and the parent pressure has the form `L(B)=C+log p(B)`.  Delete bridge edge
`e` and write

```math
r_e(B_{-e})=E_\mu[Q_e\mid B_{-e}]
```

for the exact cavity response.  For `s<=0`, put

```math
{d\widehat\Pi_s\over dU}={e^{sL}\over E_Ue^{sL}}
={p^s\over E_Up^s},
\qquad
S_s=E_{\widehat\Pi_s}\sum_e r_e^2.                 \tag{RO.2}
```

## 2. A tilted cavity-energy identity

Let `m_e(B)=E_\mu[Q_e\mid B]` be the full posterior mean.  Inserting the
deleted edge gives

```math
m_e(B)={r_e+\rho B_e\over1+\rho B_er_e}.            \tag{RO.3}
```

Conditionally on `B_(-e)`, under `widehat Pi_s`,

```math
\Pr_s(B_e=b\mid B_{-e})
={ (1+b\rho r_e)^s
  \over(1+\rho r_e)^s+(1-\rho r_e)^s}.             \tag{RO.4}
```

Put `x=r_e^2` and `z=atanh(rho r_e)`.  A direct two-point calculation gives

```math
\begin{aligned}
g_s(r_e)
&:=E_s[B_em_e(B)\mid B_{-e}]\\
&={\rho(1-x)+(1-\rho^2)r_e\tanh(sz)
   \over1-\rho^2x}.                                 \tag{RO.5}
\end{aligned}
```

For `s<=0`, the second numerator term is nonpositive.  Convexity of
`atanh` on `[0,1)` gives `|z|<=t|r_e|`, and hence

```math
|\tanh(sz)|\le |s|t|r_e|.
```

and `1-rho^2x>=1-rho^2` imply

```math
\boxed{
g_s(r_e)\ge\rho\left\{1-
 \left(1+{|s|t\over\rho}\right)r_e^2\right\}.}    \tag{RO.6}
```

This distinction from the planted/output law is exact: (RO.4)--(RO.6) are
computed under the inverse pressure tilt of the **fair** bridge law.

## 3. Support entropy bounds the tilted bridge energy

For each bridge `B`,

```math
\sum_eB_em_e(B)
\le f(B):=\max_{Q\in\operatorname {supp}\mu}\langle B,Q\rangle.
                                                               \tag{RO.7}
```

The log-sum-exp bound for fair Rademacher `B` is

```math
\log E_Ue^{af(B)}\le\log K+{a^2d\over2}.            \tag{RO.8}
```

Also, the fair-cube logarithmic-Sobolev inequality and the exact cavity
gradient give, with the repository's half-flip normalization, the
**self-consistent** entropy bound

```math
D(\widehat\Pi_s\Vert U)
\le s^2E_{\widehat\Pi_s}\Gamma_L
\le s^2t^2S_s.                                      \tag{RO.9}
```

Entropy duality applied to (RO.8), optimized over `a>0`, therefore yields

```math
E_{\widehat\Pi_s}f
\le\sqrt{2d\{\log K+D(\widehat\Pi_s\Vert U)\}}
\le\sqrt{2d\{(N-1)\log2+s^2t^2S_s\}}.             \tag{RO.10}
```

Summing (RO.6), comparing with (RO.7)--(RO.10), and rearranging proves the
main finite-order inequality.

### Theorem RO.1 (rank-one negative-tilt overlap floor)

For every actual child pair, every orientation, and every `s<=0`, put

```math
a_N={2(N-1)\log2\over \rho^2d},
\qquad b_{N,s}={2s^2t^2\over\rho^2},
\qquad c_{N,s}=1+{|s|t\over\rho}.                  \tag{RO.11}
```

If `a_N<1`, then

```math
\boxed{
{S_s\over d}\ge
\chi(a_N,b_{N,s},c_{N,s})
:={2c_{N,s}+b_{N,s}
 -\sqrt{b_{N,s}^2+4b_{N,s}c_{N,s}+4c_{N,s}^2a_N}
 \over2c_{N,s}^2}>0.}                               \tag{RO.12}
```

Indeed, with `x_s=S_s/d`, (RO.6)--(RO.10) give

```math
1-c_{N,s}x_s\le\sqrt{a_N+b_{N,s}x_s}.              \tag{RO.13}
```

The difference between the two sides is strictly decreasing in `x_s`, and
its unique zero is (RO.12).  No minimality identity is needed beyond
selecting the actual children; the mechanism is the `O(N)` logarithmic
support entropy of their exact rank-one latent word.

At the fair endpoint this specializes to

```math
\boxed{
{1\over mn}E_U\sum_e r_e^2
\ge
1-{\sqrt{2(N-1)\log2/(mn)}\over\tanh(\beta/\sqrt N)}.}
                                                               \tag{RO.14}
```

Thus, if `m/N->theta in (0,1)`, the fair normalized cavity overlap has
liminf at least

```math
\left[1-{1\over\beta}
 \sqrt{{2\log2\over\theta(1-\theta)}}\right]_+.    \tag{RO.13}
 \sqrt{{2\log2\over\theta(1-\theta)}}\right]_+.    \tag{RO.15}
```

## 4. A fixed negative interval, hence a rhohat obstruction

Assume a comparable split window

```math
\vartheta N\le m,n\le(1-\vartheta)N,
\qquad 0<\vartheta\le1/2,
```

and define

```math
\beta_c(\vartheta)
=\sqrt{{2\log2\over\vartheta(1-\vartheta)}}.        \tag{RO.16}
```

Fix `beta>beta_c(vartheta)` and `lambda>0`.  Put

```math
a={\beta_c(\vartheta)^2\over\beta^2}<1,
b=2\lambda^2,
\qquad c=1+\lambda.                                  \tag{RO.17}
```

Uniformly for every `s in [-lambda,0]`, the finite coefficients in (RO.13)
are no worse asymptotically than `a,b,c`.  Consequently

```math
\boxed{
\liminf_{N\to\infty}\widehat\rho_N^-(\lambda)
\ge\chi(a,b,c)
={2c+b-\sqrt{b^2+4bc+4c^2a}\over2c^2}>0.}          \tag{RO.18}
```

The liminf is uniform over all comparable splits, both orientations, and
all actual contracted-temperature minimizing children.  For asymptotically
balanced splits the explicit threshold is

```math
\beta_c(1/2)=\sqrt{8\log2}=2.35482\ldots;
```

for the window `[N/4,3N/4]` it is

```math
\beta_c(1/4)=\sqrt{32\log2/3}=2.71911\ldots.
```

Hence `rhohat_N^-(lambda)=o(1)`, and therefore every proposed power decay,
is false in this large-fixed-temperature regime.

## 5. What rectangle parity alone sees

The deterministic rectangle relation

```math
Q_{ij}Q_{i\ell}Q_{kj}Q_{k\ell}=1
```

places, for every edge `e`, `(m-1)(n-1)` degree-three path monomials with
coefficient `rho^3` in the edge-deleted numerator.  Their orthogonal fair
`L^2` mass is

```math
(m-1)(n-1)\rho^6=\Theta_\beta(N^{-1}).              \tag{RO.19}
```

Thus any fixed even-cycle level sees only the familiar `1/N` scale per
edge.  The constant floor in RO.11 is not a fixed-cycle estimate: it is a
nonperturbative resummation forced by the complete rank-one support having
only `exp(O(N))` states.  Above `beta_c`, a fair random bridge cannot sustain
the direct independent-edge response on all `Theta(N^2)` edges without its
Gibbs cross energy exceeding the maximum of those `exp(O(N))` rank-one
words.

## 6. Frontier consequence

RO.18 is a scalable actual-minimizer obstruction with an explicit
macroscopic cause.  It rules out `L_raw-negative-overlap` as a decay lemma
for all large fixed `beta`.  It does **not** by itself identify whether the
positive overlap creates extensive reverse-product dependence or coherent
retuning, and therefore does not supply a Level-6 recurrence.  The only
nonredundant continuation of this implementation is a directional theorem
extracting one of those branches from the low-entropy rank-one overlap
phase.
