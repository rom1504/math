# Sector-contrast quotient for canonical actual-child rows

Status: **rigorous exact quotient theorem and actual-minimizer scope
falsifier**.  This note isolates the smallest pointwise statistic needed to
compare the two orientation-conditioned canonical inverse row escorts.  The
quotient is exact, but it controls only the erased-row part of the bridge
law.  An existing pair of exact pressure-minimizing children shows that it
does not determine the joint interaction or its canonical cumulant.

## 1. Sector output likelihoods

Fix a right child signing `D` of order `n`, internal temperature `t`, and
row-channel amplitude `u`.  As in DER.3, let

```math
R_D^a(b;t,u)
=E_{Y\sim\mu_{D,a,t}}\cosh(u\langle b,Y\rangle),
\qquad a\in\{\pm1\},
```

and normalize each sector output relative to the fair row law:

```math
z_a(b)={R_D^a(b;t,u)\over(\cosh u)^n}.
                                                               \tag{SQ.1}
```

Both `z_+` and `z_-` are strictly positive and have uniform mean one.
Define their equal-sector scale and contrast by

```math
s(b)={z_+(b)+z_-(b)\over2},
\qquad
c(b)={z_+(b)-z_-(b)\over z_+(b)+z_-(b)}.            \tag{SQ.2}
```

Thus `s>0`, `-1<c<1`, and `z_a=s(1+ac)`.  Equivalently,

```math
\operatorname {atanh}c(b)={1\over2}\log{z_+(b)\over z_-(b)}.  \tag{SQ.3}
```

The function `c` is the sector Bayes factor in bounded coordinates.  It is
not assumed to have few values: as a table on the row cube it can still be
exponentially large.

For the left and right child sector biases, put

```math
\gamma_A={1\over2}\log{Z_A^+\over Z_A^-},
\qquad
\gamma_D={1\over2}\log{Z_D^+\over Z_D^-},
\qquad
\theta_\epsilon=\tanh(\gamma_D+\epsilon\gamma_A).  \tag{SQ.4}
```

## 2. Exact factorization and neutral bracketing

**Theorem SQ.1 (sector-contrast factorization).**  For the canonical erased
row of the actual two-child channel in orientation `epsilon`, its forward
output likelihood is exactly

```math
\boxed{
z_{A\to D}^{\epsilon}(b)
=s(b)\{1+\theta_\epsilon c(b)\}.}                  \tag{SQ.5}
```

The neutral augmented extension response of `D` is

```math
\boxed{
z_{D,t,u}^{0}(b)
=s(b)\{1+\tanh(\gamma_D)c(b)\}.}                   \tag{SQ.6}
```

If `Pi(epsilon)` is the forward zero-bridge orientation law, then

```math
\boxed{
z_{D,t,u}^{0}
=\Pi(+)z_{A\to D}^{+}+\Pi(-)z_{A\to D}^{-}.}       \tag{SQ.7}
```

Consequently the neutral response is bracketed pointwise:

```math
\boxed{
\min_{\epsilon}z_{A\to D}^{\epsilon}(b)
\le z_{D,t,u}^{0}(b)
\le\max_{\epsilon}z_{A\to D}^{\epsilon}(b).}      \tag{SQ.8}
```

Here `s` alone is the response with *equal sector weights*.  Unless
`gamma_D=0`, it is not the neutral augmented response in (SQ.6), whose
sector weights are proportional to `Z_D^a`.

*Proof.*  Reindex the sector in DER.19 by the right-child sign `a`.  Its
canonical weight is

```math
\widehat\omega_a^{\epsilon}
={Z_A^{\epsilon a}Z_D^a
  \over\sum_{d=\pm1}Z_A^{\epsilon d}Z_D^d}
={e^{a(\gamma_D+\epsilon\gamma_A)}
  \over2\cosh(\gamma_D+\epsilon\gamma_A)}
={1+a\theta_\epsilon\over2}.                       \tag{SQ.9}
```

Taking the mixture of `z_a=s(1+ac)` gives (SQ.5).  The neutral augmented
measure has sector weights

```math
\omega_D^a={Z_D^a\over Z_D^++Z_D^-}
={1+a\tanh\gamma_D\over2},                         \tag{SQ.10}
```

which proves (SQ.6).

Put `p=tanh gamma_A` and `q=tanh gamma_D`.  EE.22 gives

```math
\Pi(\epsilon)={1+\epsilon pq\over2},
```

while the addition formula gives

```math
\theta_+={p+q\over1+pq},
\qquad
\theta_-={q-p\over1-pq}.                           \tag{SQ.11}
```

Direct substitution yields

```math
\Pi(+)\theta_++\Pi(-)\theta_-=q.                  \tag{SQ.12}
```

Equations (SQ.5)--(SQ.6) prove the barycenter identity (SQ.7).
The weights `Pi(+)` and `Pi(-)` are strictly positive at finite
temperature, so (SQ.8) follows. `square`

The inverse powers obey the analogous *pointwise unordered* bracket

```math
\min_\epsilon (z_{A\to D}^{\epsilon})^{-\lambda}
\le (z_{D,t,u}^{0})^{-\lambda}
\le\max_\epsilon (z_{A\to D}^{\epsilon})^{-\lambda}. \tag{SQ.13}
```

There is no fixed choice of endpoint in (SQ.13): it changes with the sign
of `c(b)`.  After separately normalizing the three inverse escorts, (SQ.13)
does not give stochastic domination, a KL sandwich, or a common mode.  The
all-minimizer counterexample DER.26 already shows that neutral reinsertion
preference can reverse under a biased canonical sector law.

## 3. Minimal quotient for the two inverse row escorts

For `lambda>0`, define the two normalized canonical inverse row escorts

```math
{dr_\epsilon\over dU_n}(b)
={\{z_{A\to D}^{\epsilon}(b)\}^{-\lambda}
  \over
  \mathcal Z_\epsilon},
\qquad
\mathcal Z_\epsilon
=E_{U_n}\{z_{A\to D}^{\epsilon}\}^{-\lambda}.      \tag{SQ.14}
```

**Theorem SQ.2 (minimal binary-experiment statistic).**  If
`\gamma_A\ne0`, then

```math
\boxed{
{dr_+\over dr_-}(b)
={\mathcal Z_-\over\mathcal Z_+}
 \left{{1+\theta_-c(b)\over1+\theta_+c(b)}\right\}^{\lambda}.} \tag{SQ.15}
```

The function of `c` on the right is strictly monotone.  Hence, modulo null
sets and invertible relabelling, `c` is the minimal exact statistic for the
binary experiment `{r_+,r_-}`: it is sufficient, and every statistic which
makes the likelihood ratio measurable must also determine `c`.

If `gamma_A=0`, then `theta_+=theta_-` and the two row escorts coincide, so
the minimal statistic is trivial.

*Proof.*  Substitute (SQ.5) in (SQ.14).  The common factor `s(b)` cancels
from the pointwise likelihood ratio, giving (SQ.15).  At finite temperature
all denominators are positive.  Moreover

```math
{d\over dc}\log
 \left({1+\theta_-c\over1+\theta_+c}\right)
={\theta_--\theta_+
  \over(1+\theta_-c)(1+\theta_+c)}.                 \tag{SQ.16}
```

Since `tanh` is injective, `\gamma_A\ne0` implies
`\theta_+\ne\theta_-`.  The derivative never vanishes, so the likelihood
ratio and `c` generate the same level-set partition.  This is the standard
minimal-sufficiency criterion for a two-law finite experiment. `square`

For `m` erased rows the corresponding canonical products satisfy

```math
\log{d(r_+^{\otimes m})\over d(r_-^{\otimes m})}(B)
=m\log{\mathcal Z_-\over\mathcal Z_+}
 +\lambda\sum_{i=1}^m
  \log{1+\theta_-c(B_i)\over1+\theta_+c(B_i)}.      \tag{SQ.17}
```

Thus `c` captures the entire orientation change of the canonical iid-row
certificate.  This is an exact sector-oriented quotient, but not by itself
a low-information quotient: storing arbitrary values of `c` on
`{+-1}^n` can retain an exponential row table.  Also, the scalar
`mathcal Z_-/mathcal Z_+` depends on the joint distribution of `(s,c)`
under the fair row law; `c` is minimal as an observation statistic once
the two models are fixed, not a claim that the models can be reconstructed
from the range of `c` alone.

**Proposition SQ.2a (exact contrast mass).**  If `P_a` is the
sector-conditioned row-output law with density `z_a` relative to `U_n`,
then

```math
\boxed{
\operatorname {TV}(P_+,P_-)=E_{U_n}s|c|,
\qquad E_{U_n}sc=0.}                              \tag{SQ.17a}
```

Unless `c=0` almost surely, both `{c>0}` and `{c<0}` have positive
`sU_n`-mass.  Hence neither sector output, and therefore neither
orientation response when the two biases differ, can dominate the other
pointwise.

*Proof.*  Since `z_+-z_-=2sc`, the total-variation identity is immediate.
The equality `E_Usc=0` follows from `E_Uz_+=E_Uz_-=1`.  Since `s>0`, a
nonzero function of only one sign could not have zero `sU_n`-mean.
`square`

## 4. What the quotient does not determine

Let `p_\epsilon(B)` be the full forward bridge likelihood and let

```math
h_\epsilon(B)
=\log p_\epsilon(B)-\sum_{i=1}^m
 \log z_{A\to D}^{\epsilon}(B_i)                  \tag{SQ.18}
```

be the exact joint interaction in CR.4.  The orientation difference is

```math
\boxed{
h_+(B)-h_-(B)
=\log{p_+(B)\over p_-(B)}
 -\sum_{i=1}^m
  \log{1+\theta_+c(B_i)\over1+\theta_-c(B_i)}.}    \tag{SQ.19}
```

The second term is completely determined by the sector-contrast quotient.
The first is a genuinely joint sector-response ratio.  It depends on
multi-row correlations of the left child and on their coupling to the
right child; it is not determined by the erased-row channel.  Therefore
`(s,c,gamma_A,gamma_D)` determines both canonical row products but not
`h_epsilon`, the endpoint escort `q_epsilon`, or

```math
\mathcal J_\epsilon
=\log E_{r_\epsilon}
 \exp\{-\lambda(h_\epsilon-E_{r_\epsilon}h_\epsilon)\}. \tag{SQ.20}
```

This failure occurs already on actual pressure-minimizing children.  Let
`A_0,A_1` be the two order-eight minimizing classes in EO.4 and use the
unique order-two minimizer `D` as the common right child.  The two `A_r`
have the same complete signed energy histogram (RP.7), hence the same
`gamma_A(t)` for every `t`; the right child, and therefore its complete
functions `s` and `c`, is identical in the two experiments.  Consequently
SQ.5 gives identical canonical erased-row likelihoods for both
orientations and both `r=0,1`.

Nevertheless EO.4 proves, for all sufficiently large finite internal `t`,
that both left children are exact pressure minimizers while their
low-channel joint row-interaction coefficients are different:

```math
K_\epsilon(A_0,D;t)\ne K_\epsilon(A_1,D;t),
```

with exact zero-temperature limits `20` and `12`.  By EO.2 this gives
different leading `u^4` coefficients for both weighted row total
correlation and `mathcal J`, while the marginal-retuning contribution
starts only at order `u^8`.  Thus even the entire row-scale/contrast pair
`(s,c)` plus both scalar sector biases does not determine the joint
interaction response.

This is a fixed-order, small-channel actual-minimizer falsifier.  It is not
an asymptotic statement at `t=u=beta/sqrt(N)`.  Its exact lesson is narrower:
the sector contrast is the minimal quotient for **orientation comparison of
the canonical row escorts**, but any theorem controlling the joint
`h/J` resource needs an additional cross-row, sector-oriented statistic.

## 5. A constant-dimensional carrier for the response tangent

The missing cross-row statistic has a finite exact form at the first
nonzero order in the channel.  For a child signing `C` of order `d`, let

```math
v_C^a=\bigl(E_{\mu_{C,a,t}}X_iX_j\bigr)_{1\le i<j\le d},
\qquad
G_C(a,b)=\langle v_C^a,v_C^b\rangle,
\qquad a,b\in\{\pm1\}.                            \tag{SQ.21}
```

At a fixed internal temperature, define the sector--Gram state

```math
\mathsf S_2(C;t)
=\bigl(d,\gamma_C(t),G_C(t)\bigr),
\qquad
\gamma_C={1\over2}\log{Z_C^+\over Z_C^-}.         \tag{SQ.22}
```

Apart from the order, this contains four real numbers: one sector bias and
the three entries of a symmetric `2 by 2` Gram matrix, independently of the
number of spins.  It does not
contain the energy landscape, the sector correlation vectors themselves,
or a row table.

**Theorem SQ.3 (finite carrier for the exact interaction tangent).**  For
two children `A,D`, the oriented overlap-curvature coefficient in EO.2 is
determined exactly by their two sector--Gram states:

```math
\boxed{
K_\epsilon(A,D;t)
=\sum_{a,b=\pm1}
  \pi_a^\epsilon\pi_b^\epsilon
  G_A(a,b)\{n+2G_D(\epsilon a,\epsilon b)\},}      \tag{SQ.23}
```

where

```math
\pi_a^\epsilon
={Z_A^a Z_D^{\epsilon a}
  \over\sum_{c=\pm1}Z_A^cZ_D^{\epsilon c}}.       \tag{SQ.24}
```

Consequently the same fixed-dimensional state composes exactly to the
first nonzero joint-cancellation response:

```math
T_u(\epsilon)
={\lambda^2u^4\over2}
  \mathcal K_\epsilon(\mathsf S_2(A;t),
                       \mathsf S_2(D;t))
 +O_{A,D,t,\lambda}(u^6),
```

```math
\mathcal J_\epsilon(u)
={\lambda^2u^4\over2}
  \mathcal K_\epsilon(\mathsf S_2(A;t),
                       \mathsf S_2(D;t))
 +O_{A,D,t,\lambda}(u^6).                          \tag{SQ.25}
```

The canonical marginal-retuning term in EO.7 is `O(u^8)`.

*Proof.*  EO.2 gives

```math
K_\epsilon
=\sum_{1\le i<k\le m}\sum_{j,l=1}^n
  \bigl(\Gamma_{ik,jl}^{\epsilon}\bigr)^2.         \tag{SQ.26}
```

Condition on the shared zero-bridge sector `a`.  The two children are then
independent, and the signed sector convention gives

```math
\Gamma_{ik,jl}^{\epsilon}
=\sum_{a=\pm1}\pi_a^\epsilon
  C_A^a(i,k)C_D^{\epsilon a}(j,l),                 \tag{SQ.27}
```

where `C_C^a(r,s)=E_{\mu_{C,a,t}}X_rX_s` and
`C_C^a(r,r)=1`.  Expanding the square in (SQ.26) and interchanging sums,

```math
K_\epsilon
=\sum_{a,b}\pi_a^\epsilon\pi_b^\epsilon
 \left(\sum_{i<k}C_A^a(i,k)C_A^b(i,k)\right)
 \left(\sum_{j,l}C_D^{\epsilon a}(j,l)
                   C_D^{\epsilon b}(j,l)\right).  \tag{SQ.28}
```

The ordered-pair sum for a child of order `d` is

```math
\sum_{r,s}C_C^a(r,s)C_C^b(r,s)
=d+2G_C(a,b).                                      \tag{SQ.29}
```

The left sum in (SQ.26) is over `i<k`, so it is exactly `G_A(a,b)`;
the right sum is over all ordered `(j,l)`, so (SQ.29) gives its factor.
This proves (SQ.23).  Combining it with EO.2 proves (SQ.25). `square`

This is the first explicitly named, demonstrably sub-landscape child state
in this phase that composes to a nonzero joint resource.  Its limitation is
equally precise: (SQ.25) is an infinitesimal `u\to0` statement at fixed
children.  The remainder is not uniform when the child orders grow and
`u=t=\beta/\sqrt N`; block-parity examples show that fourth-order data do
not control a physical endpoint in general.  Nor is a finite-precision
encoding or closure under repeated composition proved.  A Level-6 route therefore
requires an optimizer-specific synchronization or tail theorem promoting
`\mathsf S_2`, or a controlled finite extension of it, from tangent scale
to the contracted-temperature path.

**Corollary SQ.4 (tangent spectral-or-harmless dichotomy).**  Put

```math
g_C=\max_{a=\pm1}G_C(a,a).                         \tag{SQ.30}
```

Then every child pair satisfies

```math
\boxed{
0\le K_\epsilon(A,D;t)
\le g_A\{|D|+2g_D\}.}                           \tag{SQ.31}
```

For comparable child orders with `m+n=N`, if `K_epsilon=o(N^3)`, the
formal physical-amplitude leading term
`(lambda^2 beta^4/(2N^2))K_epsilon` is `o(N)`.  Conversely, if
`K_epsilon>=eta N^3` along a subsequence, then for all large `N` at least
one child has a sector covariance matrix with

```math
\boxed{\lambda_{\max}(E_{\mu_{C,a,t}}XX^T)
\ge c_{\eta,\kappa}\sqrt N}                    \tag{SQ.32}
```

for some sector `a`; equivalently, one unit-vector aggregate spin
observable has variance at least `c_(eta,kappa)sqrt(N)`.

*Proof.*  Cauchy--Schwarz gives
`|G_C(a,b)|<=sqrt(G_C(a,a)G_C(b,b))<=g_C`.  The
full-correlation inner product `|D|+2G_D(a,b)` is nonnegative and at
most `|D|+2g_D`; applying these bounds termwise in (SQ.23), whose total is
the squared norm `K_epsilon>=0`, gives (SQ.31).
If both `g_A` and `g_D` were `o(N^(3/2))`, its right side would be
`o(N^3)`.  Quantitatively, writing `g=max(g_A,g_D)`,
`K_epsilon<=g(N+2g)`; a lower bound `eta N^3` therefore forces
`g>=c_eta N^(3/2)` for large `N`.  In the witnessing sector, the
correlation matrix `C=E[XX^T]` is positive semidefinite with trace equal
to the child order `d`, and

```math
\operatorname {tr}C^2=d+2G_C(a,a).               \tag{SQ.33}
```

Thus `lambda_max(C)>=tr(C^2)/tr(C)`, which is at least a constant times
`sqrt(N)` when `d` is comparable to `N`. `square`

This corollary turns an extensive **tangent** into a named aggregate
child mode rather than an arbitrary bridge oracle.  It still does not
exclude a physical-scale phase generated first by higher-order parity
structure with small pair covariance, so it is not a physical endpoint
dichotomy.

## 6. Exact contrast-fibre disintegration of the canonical row reference

The contrast quotient separates canonical row-orientation dependence from
a common reference kernel at the full physical amplitude.  It does not
remove orientation dependence from the full escort.  Let

```math
C_i=c(B_i),\qquad \mathbf C=(C_1,\ldots,C_m),      \tag{SQ.34}
```

and let `R_epsilon` and `Q_epsilon` denote the laws of `mathbf C` under
the canonical product `r_epsilon^{otimes m}` and the full inverse escort
`q_epsilon`, respectively.  For a value `z` in the finite range of `c`,
define

```math
\kappa_z(b)
={\boldsymbol 1_{\{c(b)=z\}}s(b)^{-\lambda}U_n(b)
  \over
  E_{U_n}[\boldsymbol 1_{\{c=z\}}s^{-\lambda}]}. \tag{SQ.35}
```

**Theorem SQ.5 (canonical-row orientation fibre dichotomy).**  For either
orientation and every contrast vector in the support,

```math
\boxed{
r_\epsilon^{\otimes m}(dB\mid\mathbf C)
=\bigotimes_{i=1}^m\kappa_{C_i}(dB_i),}          \tag{SQ.36}
```

and the right side is independent of `epsilon`.  The exact canonical
cumulant therefore disintegrates as

```math
\boxed{
\mathcal J_\epsilon
=D(R_\epsilon\Vert Q_\epsilon)
 +E_{R_\epsilon}
   D\!\left(\bigotimes_i\kappa_{C_i}
       \middle\Vert q_\epsilon(\,\cdot\mid\mathbf C)\right).} \tag{SQ.37}
```

Both terms are nonnegative.  Hence if `J_epsilon>=eta N`, then at least
one of the following holds:

1. the scalar row feature `c(B_i)` carries a reverse-product image gap at
   least `eta N/2`; or
2. after the canonical row-orientation statistic has been revealed, a
   common orientation-independent reference product still has average
   conditional KL at least `eta N/2` from the actual bridge escort.

*Proof.*  On a fibre `c(b)=z`, (SQ.5) and (SQ.14) give

```math
r_\epsilon(b)\ \propto\
s(b)^{-\lambda}(1+\theta_\epsilon z)^{-\lambda}U_n(b). \tag{SQ.38}
```

The second factor is constant on the fibre, so normalization leaves
exactly `kappa_z`, independently of orientation.  Rowwise conditioning of
the product law proves (SQ.36).  Apply the ordinary KL chain rule to the
deterministic map `B mapsto C` to obtain (SQ.37). `square`

The first term in (SQ.37) is a concrete high-transport image certificate;
unlike a fixed coordinate or Walsh parity, `c` can use the whole row.
The second term is a more localized conditional object, but no smaller
information footprint is yet proved: only canonical row orientation has
been exhausted, while the conditional full escort can retain joint
orientation through (SQ.19).  This is a classification rather than a
no-gain theorem.  The
table `c` may have exponentially many values, and no theorem currently
makes the conditional residual sublinear or gives it a tight latent-product
representation.  When `gamma_A=0`, the two canonical row escorts already
coincide and their minimal orientation statistic is trivial; conditioning
on a nonconstant `c` is then optional rather than minimal.

**Theorem SQ.6 (polynomial approximate orientation quotient).**  Suppose
the opposite child `A` is an exact order-`m` pressure minimizer at
`t=beta/sqrt(N)`, with `m<=N`.  Let

```math
\ell(b)=\log{dr_+\over dr_-}(b).               \tag{SQ.39}
```

Its oscillation obeys

```math
\boxed{
\operatorname {osc}\ell
\le4\lambda|\gamma_A|
\le4\lambda N\left(\log2+{\beta^2\over4}\right).} \tag{SQ.40}
```

For every `eta>0`, partition the range of `ell` into intervals of length
at most `eta`, and let `T_eta(b)` be the interval label.  The quotient
has at most

```math
1+\left\lceil{4\lambda N(\log2+\beta^2/4)\over\eta}\right\rceil
                                                               \tag{SQ.41}
```

states, and on every nonempty cell

```math
\boxed{
D_\infty(r_+(\,\cdot\mid T_\eta)
 \Vert r_-(\,\cdot\mid T_\eta))\le\eta,
\qquad
D_\infty(r_-(\,\cdot\mid T_\eta)
 \Vert r_+(\,\cdot\mid T_\eta))\le\eta.}   \tag{SQ.42}
```

For `m` canonical rows, conditioning on their label vector gives the
same two-sided bound `m eta` between the two product kernels.  Thus any
choice `eta_N=o(1)` removes canonical sector orientation at `o(N)` KL
cost using only a polynomial-size alphabet per row; for example
`eta_N=N^(-alpha)` uses `O(N^(1+alpha))` labels.

*Proof.*  The nonconstant part of (SQ.15) is

```math
\lambda\{\log(1+\theta_-c)-\log(1+\theta_+c)\}. \tag{SQ.43}
```

It is monotone, and its endpoint difference is `-4lambda gamma_A`,
because `(1+tanh x)/(1-tanh x)=e^(2x)`.  This proves the first
inequality in (SQ.40).  Also `|gamma_A|<=tQ(A)`.  A maximizing spin and
its global negative give

```math
F_A(t)=\log E_x\cosh(tH_A(x))\ge tQ(A)-m\log2. \tag{SQ.44}
```

Exact pressure minimality and annealed averaging give
`F_A(t)<=binom(m,2)log cosh(t)`, while
`log cosh(t)<=t^2/2`; this proves the second inequality in (SQ.40).

On a quotient cell, the likelihood ratio `L=dr_+/dr_-` has
`max L/min L<=e^eta`.  The conditional likelihood ratio is
`L/E_(r_-|cell)L`, so it lies between `e^(-eta)` and `e^eta`.
This proves (SQ.42).  Tensorization over the rowwise cells adds the
max divergences. `square`

For the label vector `mathbf T_eta`, the ordinary chain rule also gives
the full-amplitude split

```math
\mathcal J_\epsilon
=D((r_\epsilon^{\otimes m})_{\mathbf T_\eta}
       \Vert(q_\epsilon)_{\mathbf T_\eta})
 +E D(r_\epsilon^{\otimes m}(\,\cdot\mid\mathbf T_\eta)
       \Vert q_\epsilon(\,\cdot\mid\mathbf T_\eta)). \tag{SQ.45}
```

If the first term is sublinear while `J_epsilon` is linear, the conditional
residual is linear; by (SQ.42), changing the canonical orientation changes
its reference kernel by only `o(N)` max divergence when `eta_N=o(1)`.

The theorem rules out a superpolynomial **per-row canonical orientation
alphabet** at any inverse-polynomial accuracy.  It does not compress the
vector of row labels, compute the map without the one-row response table,
or control the second term in
(SQ.37).  Extensive complexity can therefore survive only in the
contrast-image dynamics across rows or in the orientation-blind
canonical reference versus the within-cell full-escort residual.

## 7. Frontier implication

The quotient identifies the first exact scalar coordinate in which sector
orientation acts:

```text
row orientation response = common scale s + scalar sector contrast c.
```

It cleanly separates the biased-extension issue from the joint-dependence
issue.  What remains is not another scalar correction to `c`; the EO.4
collision shows that multi-row oriented overlap information is already
visible while `(s,c)` is held fixed.  A useful next state must therefore
augment the row quotient with a cross-row statistic, yet remain strictly
smaller than the complete joint bridge likelihood.  The finite tensor in
EO.2 is the exact tangent candidate, but no uniform physical-scale closure
is currently proved.

The independent audit of SQ.3--SQ.6, including the information-footprint
qualification for the fibre quotient, is
[`../audits/actual_child_sector_gram_tangent_carrier_audit.md`](../audits/actual_child_sector_gram_tangent_carrier_audit.md).
