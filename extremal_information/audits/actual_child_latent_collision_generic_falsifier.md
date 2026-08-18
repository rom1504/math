# A maximally diffuse generic channel can have exponential latent collision

**Status.**  Rigorous scalable **generic, non-actual** falsifier for deriving
`log overline K_del=o(N)` from the currently available rank-one, spread,
flip-inequality, and strong-channel inputs alone.  The latent prior is the
uniform projective rank-one law.  It has maximal conditional entropy on both
factors.  Nevertheless, in an explicit nonempty strong-channel endpoint
window, its inverse-escort-averaged deleted posterior collision is
`exp{Omega(N)}`.

The proof has two parts.  First, a deterministic entropy--barycenter
inequality says that a posterior with macroscopic rank-one mean must spend
linear KL relative to the uniform projective prior.  Second, the existing
bipartite-ground-state overlap inequality supplies such a mean at the actual
negative endpoint.  No conference or Paley surrogate is used.

This does **not** construct actual pressure-minimizing sign children.  Uniform
factor laws are zero-interaction Gibbs laws, not Gibbs laws of hollow
`{+-1}` signings.  The conclusion is therefore a ceiling on arguments using
only the structural consequences listed in Section 6, not a falsification of
an optimizer-specific collision theorem.

## 1. Uniform projective rank-one prior

Put

```math
d=mn,qquad N=m+n,qquad t={\beta\over\sqrt N},
\qquad \rho=\tanh t.
```

Let `X` and `Y` be independent fair sign vectors of lengths `m,n`, and let
`mu` be the induced uniform law on the `2^(N-1)` projective matrices

```math
Q=XY^{\mathsf T}.
```

Here "projective" quotients only the factor gauge
`(x,y)sim(-x,-y)`; the two channel atoms `Q` and `-Q` remain distinct.

The normalized binary-channel likelihood and its inverse escort are

```math
P_t(B)=E_\mu\prod_{e=1}^d(1+\rho B_eQ_e),
\qquad
{dq_\lambda\over dU}
={P_t^{-\lambda}\over E_UP_t^{-\lambda}}.        \tag{LC.1}
```

Let `nu_B` be the complete forward posterior, let

```math
M(B)=E_{\nu_B}Q,
```

and let `K_0(B)` and `K_e(B_(-e))` be the complete and deleted collision
factors from (GC.9)--(GC.11).  Thus

```math
K_0(B)=\exp D_2(\nu_B\Vert\mu).                  \tag{LC.2}
```

## 2. A bilinear Rademacher MGF bound

### Lemma LC.1 (exact determinant domination)

For every real `m by n` matrix `A`, independent fair sign vectors `X,Y`,
and every real `theta` satisfying
`theta^2||A||_op^2<1`,

```math
\boxed{
\log E_{X,Y}e^{\theta X^{\mathsf T}AY}
\le {\theta^2\|A\|_F^2
       \over2\{1-\theta^2\|A\|_{\rm op}^2\}}.}   \tag{LC.3}
```

*Proof.*  Conditional on `X`, `log cosh z<=z^2/2` gives

```math
E_Ye^{\theta X^{\mathsf T}AY}
\le e^{\theta^2\|A^{\mathsf T}X\|_2^2/2}.
```

Introduce an independent standard Gaussian vector `g in R^n`.  The
Gaussian linearization identity and a second use of `log cosh z<=z^2/2`
give

```math
\begin{aligned}
E_Xe^{\theta^2\|A^{\mathsf T}X\|^2/2}
&=E_gE_Xe^{\theta X^{\mathsf T}Ag}\\
&\le E_ge^{\theta^2\|Ag\|^2/2}\\
&=\det(I-\theta^2A^{\mathsf T}A)^{-1/2}.
\end{aligned}                                      \tag{LC.4}
```

If `s_j` are the singular values of `A`, then

```math
-{1\over2}\sum_j\log(1-\theta^2s_j^2)
\le {\theta^2\sum_js_j^2
       \over2(1-\theta^2\|A\|_{\rm op}^2)},
```

which proves (LC.3).  `square`

## 3. Posterior barycenter forces projective information

### Theorem LC.2 (entropy--barycenter inequality)

Let `nu` be any law on projective rank-one sign matrices, absolutely
continuous with respect to the uniform projective law `mu`, and put

```math
A=E_\nu Q.
```

Then

```math
\boxed{
D(\nu\Vert\mu)
\ge {\|A\|_F^2\over2\sqrt d}
     -{\|A\|_F^2\over6d}
\ge {\|A\|_F^2\over3\sqrt d}.}                  \tag{LC.5}
```

Consequently, pointwise at every bridge word,

```math
\boxed{
K_0(B)\ge
\exp\left\{{\|M(B)\|_F^2\over3\sqrt d}\right\}.} \tag{LC.6}
```

*Proof.*  Lift `nu` to the full product cube by assigning half of the mass
of each projective atom `Q=xy^T` to each of its two representatives
`(x,y)` and `(-x,-y)`.  The uniform product law has the same two-to-one
disintegration.  Therefore the lift `nu_tilde` satisfies exactly

```math
D(\widetilde\nu\Vert U_X\otimes U_Y)
=D(\nu\Vert\mu),
\qquad E_{\widetilde\nu}XY^{\mathsf T}=A.         \tag{LC.7}
```

There is no extra `log 2`: the posterior/reference density ratio is
constant on each two-point gauge fibre.

Apply the entropy variational inequality with the test
`theta X^TAY`, where

```math
\theta={1\over2\sqrt d}.
```

Because `A` is a convex combination of sign matrices,
`||A||_op<=sqrt(d)`.  Moreover,

```math
E_{\widetilde\nu}X^{\mathsf T}AY
=\langle A,E_{\widetilde\nu}XY^{\mathsf T}\rangle
=\|A\|_F^2.                                      \tag{LC.8}
```

Lemma LC.1 and
`theta^2||A||_op^2<=1/4` therefore give

```math
\log E_{U_X,U_Y}e^{\theta X^{\mathsf T}AY}
\le {2\over3}\theta^2\|A\|_F^2
={\|A\|_F^2\over6d}.                             \tag{LC.9}
```

Equations (LC.8)--(LC.9) prove the first inequality in (LC.5).  Since
`sqrt(d)>=1`, its coefficient is at least `1/(3sqrt(d))`, proving the
second.  Finally `D_2>=D` and (LC.2) give (LC.6).  `square`

The scale in (LC.5) is the important point: a barycenter with squared
Frobenius norm `Theta(d)` costs `Theta(sqrt(d))=Theta(N)` nats for a
comparable split.  The antipodal two-word example does not contradict this;
its reference prior has only two atoms rather than `2^(N-1)` uniform atoms.

## 4. Collision lower bound from an endpoint overlap

Let `r(B)` be the matrix of exact deleted-edge cavity responses evaluated
at the complete bridge word.  The exact insertion comparison gives

```math
\|M(B)-r(B)\|_F^2\le4\rho^2d.                    \tag{LC.10}
```

Hence, if a bridge law `q` obeys

```math
{1\over d}E_q\|r(B)\|_F^2\ge\eta,               \tag{LC.11}
```

then

```math
E_q\|M(B)\|_F^2
\ge(\eta/2-4\rho^2)d.                            \tag{LC.12}
```

Jensen and Theorem LC.2 now imply

```math
\boxed{
\log E_qK_0(B)
\ge {\eta/2-4\rho^2\over3}\sqrt d.}             \tag{LC.13}
```

The one-edge comparison in Lemma GC.3 is pointwise:

```math
K_0(B)\le e^{4t}K_e(B_{-e}).
```

Therefore the deleted average satisfies the fully explicit implication

```math
\boxed{
\log\overline K_{\rm del}
\ge-4t+{\eta/2-4\rho^2\over3}\sqrt d.}           \tag{LC.14}
```

This conversion uses the annealed collision factor itself, not merely the
average of its logarithm.

## 5. An explicit negative-endpoint window

Although the headline of Theorem 37.56 names actual children, direct
inspection of its proof shows that the bipartite-ground-state overlap
estimate (37.184) uses only
the rank-one channel, the one-bit oscillation of its likelihood, and the
bipartite sign ground-state envelope.  It therefore applies to the uniform
projective prior above without any optimizer assumption.

Assume

```math
{mn\over N^2}\ge\gamma_0>0,
```

and define

```math
\beta_{\rm BG}(\gamma_0)
=\sqrt{2\over\pi}
 \sqrt{\gamma_0^{-1}+2\gamma_0^{-1/2}}.           \tag{LC.15}
```

Fix an endpoint parameter satisfying

```math
\boxed{
0<\lambda<1,
\qquad
\beta>{\beta_{\rm BG}(\gamma_0)\over1-\lambda}.} \tag{LC.16}
```

Taking `delta=lambda` in the pointwise estimate (37.184) gives, at the
actual endpoint `q_lambda=q_(-lambda)`,

```math
\liminf_N{1\over d}E_{q_\lambda}\|r(B)\|_F^2
\ge\eta_*,
\qquad
\eta_*={1-\beta_{\rm BG}(\gamma_0)/\beta-\lambda
             \over1+\lambda}>0.                  \tag{LC.17}
```

Combining (LC.14) and (LC.17) proves the scalable collision obstruction

```math
\boxed{
\liminf_{N\to\infty}{1\over N}
 \log\overline K_{\rm del}
\ge {\sqrt{\gamma_0}\over6}
 {1-\beta_{\rm BG}(\gamma_0)/\beta-\lambda
       \over1+\lambda}>0.}                       \tag{LC.18}
```

At balanced splits, `gamma_0=1/4` and
`beta_BG=4/sqrt(pi)`.  Thus (LC.18) holds, for example, for every

```math
0<\lambda<1,
\qquad \beta>{4\over(1-\lambda)\sqrt\pi}.        \tag{LC.19}
```

The endpoint restriction is exact for this proof.  The current overlap
inequality loses the coefficient `lambda`; it gives no positive endpoint
constant at `lambda>=1`.  No claim about `lambda=1` is made here.

## 6. Which actual-child consequences the example passes

The uniform projective prior has stronger diffuse properties than the
currently proved actual-child spread estimates:

1. each factor is uniform, so every factor subset, after conditioning on
   every exterior coordinate, has maximum atom exactly `2^(-|U|)`;
2. every projective Hamming cap has the corresponding uniform-cube
   exponential bound;
3. the row/column switching group is transitive, and every
   likelihood-dependent bridge tilt has averaged latent posterior equal to
   the prior;
4. regarded as a zero-interaction child law, it satisfies all flip-moment
   inequalities `E exp(-2t sum_(e in S)Y_e)>=1`, the optimizer tangent
   inequalities, the annealed pressure bound, and the entropy/cap bounds
   with slack or equality;
5. its rank-one support size is exactly the maximal permitted
   `2^(N-1)`.

Yet (LC.18) says that its deleted doubled-temperature ratio has a fixed
positive exponential rate.  Consequently none of central symmetry,
rank-one support, maximal sector spread, narrow-cap decay, averaged-
posterior nonretuning, or the **inequality directions** of the optimizer
flip identities can imply

```math
\log\overline K_{\rm del}=o(N).                  \tag{LC.20}
```

The scope boundary is essential.  A zero-interaction law is not generated
by a hollow sign matrix, so the example need not obey an unknown equality
or nonradial rigidity specific to exact pressure-minimizing signings.  A
positive actual-child theorem must use precisely such extra information--
for example quantitative nonradial flip **values**, an optimizer-specific
two-temperature cancellation, or another identity absent from the abstract
uniform law.

Thus the scalar SML remains open for actual children, but it cannot be
settled positively by recombining only the structural inequalities already
proved.
