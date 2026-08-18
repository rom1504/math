# A diffuse common phase survives factorwise macroscopic spread

Status: **rigorous scalable generic counterexample**.  This is not an
actual minimizing-child channel.  It strengthens the partial-common-block
examples by removing every frozen macroscopic factor block.  Both factor
laws have exponential conditional min-entropy on every subset, even after
conditioning on all complementary coordinates, while the canonical inverse
row product is still linearly worse than a different row product.

All logarithms are natural.  The construction is deliberately tested at the
physical scale `u=beta/sqrt(N)`.

## 1. A BSC-smoothed common phase

Let `m+n=N`, with `m/N -> theta in (0,1)`, and fix
`beta,lambda>0`.  Choose fixed crossover probabilities

```math
0<\epsilon_L,\epsilon_R<\frac12,
\qquad a=1-2\epsilon_L,
\qquad b=1-2\epsilon_R.                              \tag{DC.1}
```

Let `sigma,tau` be independent fair signs.  Independently, let the
coordinates of `xi in {+-1}^m` and `eta in {+-1}^n` be BSC noises with

```math
\Pr\{\xi_i=-1\}=\epsilon_L,
\qquad
\Pr\{\eta_j=-1\}=\epsilon_R.                        \tag{DC.2}
```

Put

```math
X=\sigma\xi,
\qquad Y=\tau\eta,
\qquad Q=XY^{\mathsf T}.                             \tag{DC.3}
```

This is a genuinely diffuse common phase.  If `U subseteq [m]`, then for
every exterior value having positive probability,

```math
\boxed{
 \left\|\mathcal L(X_U\mid X_{U^c})\right\|_\infty
 \le(1-\epsilon_L)^{|U|}.}                           \tag{DC.4}
```

Indeed, conditioning on the exterior only changes the posterior weight of
`sigma`.  Conditional on either value of `sigma`, the coordinates in `U`
remain independent and every atom has probability at most
`(1-epsilon_L)^|U|`; mixing the two components preserves the same bound.
The analogous statement for `Y` has rate `-log(1-epsilon_R)`.  Thus every
positive-density marginal and conditional marginal of **both** factors has
a fixed exponential min-entropy rate.  In particular, this construction
passes the conclusion of Theorem 37.64 rather than merely its whole-word
version.

For `u=beta/sqrt(N)`, the forward bridge likelihood relative to the fair
bridge law `U_(mn)` is

```math
p(B)=\frac{\mathbb E_{\xi,\eta}
             \cosh\{u\xi^{\mathsf T}B\eta\}}
            {(\cosh u)^{mn}}.                        \tag{DC.5}
```

Let its inverse escort and canonical inverse row product be

```math
\frac{dq}{dU_{mn}}=\frac{p^{-\lambda}}{E_Up^{-\lambda}},
\qquad r=\nu_N^{\otimes m}.                          \tag{DC.6}
```

We next identify `nu_N` exactly.

## 2. Exact erased-row law

Write `rho=tanh u`,

```math
v_N=\operatorname {arctanh}(b\rho),
\qquad S(z)=\sum_{j=1}^nz_j.                         \tag{DC.7}
```

The latent row word is a fair global sign times `eta`.  Consequently its
output likelihood is

```math
\begin{aligned}
p_i(z)
 &=\frac{E_\eta\cosh\{u\langle z,\eta\rangle\}}
          {(\cosh u)^n}\\
 &=(1-b^2\rho^2)^{n/2}\cosh\{v_NS(z)\}.             \tag{DC.8}
\end{aligned}
```

The constant disappears under inverse escort, so

```math
\boxed{
 \frac{d\nu_N}{dU_n}(z)
 =\frac{\cosh\{v_NS(z)\}^{-\lambda}}{Z_{n,N}},
 \qquad
 Z_{n,N}=E_{U_n}\cosh(v_NS)^{-\lambda}.}             \tag{DC.9}
```

Let `d_N=D(nu_N||U_n)`.  If `n/N -> 1-theta`, the triangular-array CLT and
bounded convergence give

```math
d_N\longrightarrow d_0(b)
=D\!\left(
 \frac{\cosh(G_b)^{-\lambda}}
      {E\cosh(G_b)^{-\lambda}}\,\gamma_b
 \middle\Vert\gamma_b\right)>0,                    \tag{DC.10}
```

where `G_b` is Gaussian with variance
`b^2 beta^2(1-theta)` and `gamma_b` denotes its law.  (The passage to the
limit uses the CLT together with uniform sub-Gaussian integrability, not
pointwise boundedness of `log cosh`.)  Positivity follows
because `b,beta,lambda>0` make the tilt nonconstant.  Moreover `d_0(b)` is
continuous at `b=1`.

## 3. A uniform noisy-to-hard comparison

The key estimate is that BSC smoothing changes the expectation of the log
likelihood by only `o_(epsilon)(N)`, uniformly for the two bridge products
that will be compared.

Let

```math
p_0(B)=\frac{\cosh\{u\langle B,\mathbf1\mathbf1^T\rangle\}}
             {(\cosh u)^{mn}},
\qquad \delta(B)=\log p(B)-\log p_0(B).              \tag{DC.11}
```

Set `theta_N=m/N`, `nu_N=n/N` (only in this paragraph `nu_N` is a scalar),

```math
\begin{aligned}
A_N={}&-\theta_N\log(1-\epsilon_L)
       -\nu_N\log(1-\epsilon_R),\\
C={}&2(1+\lambda)\beta^2,\\
R_N={}&\theta_N\log\!\left[1+\epsilon_L
                  (e^{C\nu_N}-1)\right]\\
    &+\nu_N\log\!\left[1+\epsilon_R
                  (e^{C\theta_N}-1)\right].          \tag{DC.12}
\end{aligned}
```

Then, for either `P=U_(mn)` or `P=r`,

```math
\boxed{-A_NN\le E_P\delta\le R_NN+\log2.}            \tag{DC.13}
```

Here is the proof.  The latent word `Q_0=11^T` occurs in (DC.5) with
probability at least
`w_0=(1-epsilon_L)^m(1-epsilon_R)^n`, which gives the lower bound.
For the upper bound, if `Delta_Q=Q-Q_0`, then

```math
\cosh\{u\langle B,Q\rangle\}
\le \cosh\{u\langle B,Q_0\rangle\}
       e^{u|\langle B,\Delta_Q\rangle|}.              \tag{DC.14}
```

Under `r`, every sequential bridge-bit conditional mean has magnitude at
most

```math
\kappa_N=\tanh(\lambda v_N)\le\lambda u.             \tag{DC.15}
```

To see this, flipping one bit changes the logarithm of the row density in
(DC.9) by at most `2lambda v_N`; the same ratio bound survives arbitrary
conditioning and marginalization.  Iterating conditional moment-generating
functions therefore gives, for every fixed array `Delta`,

```math
E_r e^{u\langle B,\Delta\rangle}
\le\exp\left\{
 {u^2\over2}\|\Delta\|_2^2
 +\kappa_Nu\|\Delta\|_1\right\}.                    \tag{DC.16}
```

The same bound holds under `U` with `kappa_N=0`.  If `d(Q,Q_0)` is Hamming
distance, then `Delta_Q` has `d` nonzero entries, all of magnitude two.
Combining Jensen, the two signs in the absolute value, (DC.15), and
(DC.16) yields

```math
E_P\delta
\le\log2+log E_Q\exp\{C d(Q,Q_0)/N\}.              \tag{DC.17}
```

If `F` and `G` count the flipped coordinates of `xi` and `eta`, then

```math
d(Q,Q_0)=nF+mG-2FG\le nF+mG.                         \tag{DC.18}
```

Taking the binomial moment-generating functions proves the upper bound in
(DC.13).  Notice that

```math
\lim_{\epsilon_L,\epsilon_R\downarrow0}
 \limsup_N(A_N+R_N)=0.                               \tag{DC.19}
```

This is the only continuity input; it is quantitative and does not appeal
to a thermodynamic limit for the noisy channel.

## 4. Linear canonical-minus-best-product regret

Let

```math
J=D(r\Vert q),
\qquad
I^\leftarrow=\inf_{P=\otimes_iP_i}D(P\Vert q).       \tag{DC.20}
```

The fair bridge law is an admissible row product.  Cancelling the common
normalizer of `q` gives the exact comparison

```math
J-I^\leftarrow
\ge J-D(U\Vert q)
=m d_N+\lambda(E_r-E_U)\log p.                      \tag{DC.21}
```

For the hard likelihood in (DC.11), positivity of `log cosh` and
`E_U|sum_(ij)B_(ij)|<=sqrt(mn)` give

```math
(E_r-E_U)\log p_0\ge-u\sqrt{mn}.                    \tag{DC.22}
```

Together with (DC.13),

```math
\boxed{
J-I^\leftarrow
\ge m d_N-\lambda\left{
 (A_N+R_N)N+\log2+u\sqrt{mn}\right}.}              \tag{DC.23}
```

Because `d_0(1)>0`, continuity in (DC.10) and (DC.19) imply that there are
fixed **positive** crossover probabilities `epsilon_L,epsilon_R` and a
constant `c=c(beta,lambda,theta)>0` for which

```math
\boxed{J-I^\leftarrow\ge cN}                        \tag{DC.24}
```

for all sufficiently large `N`.  Explicitly, it is enough that

```math
\theta d_0(1-2\epsilon_R)
>\lambda\limsup_N(A_N+R_N).                          \tag{DC.25}
```

The square-root term in (DC.23) is negligible.  Thus a diffuse BSC common
phase simultaneously has factorwise conditional min-entropy on every
macroscopic set and extensive coherent canonical product retuning.

## 5. Posterior retuning and exact scope

This is not the archived weak-common-latent directionality example in
`drafts/actual_child_response_directedness_ceiling.md`.  That example starts
from an already constructed bridge law with a constant fair-product
background and consequently has `I^leftarrow=O(1)`.  Here the common latent
is placed in the two **rank-one factors before the physical bridge channel**;
the exact erased-row inverse escort then incurs (DC.10), and (DC.24) proves
linear canonical-minus-best-product regret.

Unlike the subgroup examples in Theorems 37.58 and 37.63, the BSC-weighted
latent prior is not invariant under the full rank-one switching group.
Consequently its averaged forward posterior is not fixed by every
likelihood-dependent bridge tilt.  This is a real distinction, not merely a
missing symmetry proof.  Already for one row and two columns, the two
projective words have prior masses

```math
p_*=(1+b^2)/2,
\qquad 1-p_*=(1-b^2)/2.                              \tag{DC.26}
```

At `b=1/2`, `u=1`, and `lambda=1`, direct substitution in the two-output
Bayes formula gives averaged posterior mass

```math
\bar p_*=0.5444556798\ldots
\ne p_*=0.625.                                      \tag{DC.27}
```

For reproducibility, put

```math
k={\cosh(2u)\over\cosh(u)^2},\quad
\ell={1\over\cosh(u)^2},\quad
A=p_*k+(1-p_*)\ell,\quad
D=p_*\ell+(1-p_*)k.
```

Then the number in (DC.27) is exactly

```math
{p_*\{kA^{-(\lambda+1)}+\ell D^{-(\lambda+1)}\}
 \over
 A^{-\lambda}+D^{-\lambda}}.
```

Hence this counterexample does **not** show that the energy-shell or another
low-information posterior observable must fail.  It proves the sharper
generic ceiling

```math
\boxed{
\begin{gathered}
\text{uniform exponential conditional spread on every factor subset}
\\
\centernot\Longrightarrow
J-I^\leftarrow=o(N).
\end{gathered}}                                     \tag{DC.28}
```

In fact the retuning in this example has an exact polynomial-state
quotient.  The group `S_m times S_n` acts by row and column permutations.
The prior, fair bridge law, likelihood, inverse escort, and Bayes channel
are all equivariant under this action.  A projective rank-one word has an
orbit determined by

```math
K=\min\{\#\{i:\xi_i=-1\},m-\#\{i:\xi_i=-1\}\},
\qquad
L=\min\{\#\{j:\eta_j=-1\},n-\#\{j:\eta_j=-1\}\}.   \tag{DC.29}
```

Both `mu` and the averaged posterior `bar mu` are uniform conditional on
`(K,L)`.  The KL chain rule therefore gives the exact identity

```math
\boxed{
D(\bar\mu\Vert\mu)
=D\!\left(\mathcal L_{\bar\mu}(K,L)
          \middle\Vert\mathcal L_\mu(K,L)\right).}  \tag{DC.30}
```

There are only `O(mn)` such states (and at most one extra orientation bit if
one keeps signed rather than projective words).  Thus the construction
falsifies **spread-only** closure but does not hide its phase in the full
exponential landscape.  Its coherent variable is exposed by the two factor
magnetizations.  An actual-child analogue would need either a comparably
small optimizer-specific quotient or a proof that no such diffuse phase can
occur.

The common phase is now diffuse rather than a frozen block.  Therefore an
actual-child closure cannot follow from Theorem 37.64 alone.  It must use a
further optimizer-specific identity to control or expose this diffuse
phase--for example, a bound on posterior retuning in a genuinely coarser
quotient--and then separately connect that observable to row-product
regret.  The construction remains generic and makes no claim that BSC
factor laws can be exact quadratic Gibbs sectors of minimizing children.
