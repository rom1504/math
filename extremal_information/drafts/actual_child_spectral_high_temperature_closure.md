# Spectral high-temperature closure for actual-child row ANOVA

Status: **task-local rigorous theorem and scalable ceiling**.  This note asks
whether the infinitesimal row-ANOVA classification can be transported to the
physical bridge amplitude uniformly in order.  It gives a positive answer in
an explicit strict-high-temperature spectral regime, for the actual
contracted-temperature child minimizers.  It then shows why asymptotic
pressure near-minimality alone cannot verify the required spectral premise.
No conference or Paley surrogate is used.

We use the setup of
[`actual_child_row_anova_infinitesimal.md`](actual_child_row_anova_infinitesimal.md).
Let `m+n=N`, `t=beta/sqrt(N)`, fix an orientation `epsilon`, and let

```math
L_\epsilon(B)=\log\overline Z_N(A,\epsilon D,B;t),
\qquad B\in\{\pm1\}^{m\times n}.                     \tag{SH.1}
```

The children `A,D` below are actual minimizers of their raw-temperature
pressures.

## 1. The weaker typical-response criterion

The operator norm is only one way to verify the response property actually
used in the proof.

**Lemma SH.0 (quantitative typical-gradient closure).**  Suppose there are
convex sets `K_N` in the real bridge space, fixed constants `c,h>0`, and
numbers `R_N` such that

```math
\Pr_{U_B}(B\notin\mathcal K_N)\le e^{-cN},
\quad
\sup_{B\in\mathcal K_N}\|\nabla_BL_\epsilon(B)\|_F\le R_N,
\quad
\mathbb E_BL_\epsilon(B)\le hN,                       \tag{SH.0}
```

and `L_epsilon` has polynomial range (as every signing log-partition does).
Also assume `L_epsilon>=0`, again automatic here.
Then

```math
\operatorname{Var}_{U_B}(L_\epsilon)=O(R_N^2)+o(1),
\qquad \sigma_{\rm cross}^2(L_\epsilon)=O(R_N^2)+o(1),          \tag{SH.0a}
```

and every fixed `0<lambda<c/h` satisfies

```math
\boxed{G_{\lambda,\epsilon}
\le C_0\lambda R_N^2+o(N),}                          \tag{SH.0b}
```

for a universal concentration constant `C_0`.  In particular,
`R_N=o(sqrt(N))` excludes a linear negative-disorder gain.

*Proof.*  Take the supremum of the supporting affine functions to the convex
function `L_epsilon` over `K_N`.  This gives a convex `R_N`-Lipschitz extension
`g` agreeing with `L_epsilon` on `K_N`.  Convex-Lipschitz concentration gives
a bound `log E exp[z(g-Eg)]<=C_0z^2R_N^2`.  Polynomial range and
the exceptional probability give `E|L_epsilon-g|^2=o(1)`, proving (SH.0a).
For the negative moment, its exceptional contribution after centering is at
most `exp[-cN+lambda hN]`; the good contribution is bounded by the
subgaussian moment of `g`.  Jensen supplies the lower bound one.  This proves
(SH.0b). `square`

For the actual Gibbs law,

```math
\nabla_BL_\epsilon(B)
=t\,\mathbb E_{\nu_{\epsilon,B}}[\tau XY^{\mathsf T}].          \tag{SH.0c}
```

Thus the intrinsic premise in (SH.0), before introducing spectra, is the
typical parent-response estimate

```math
\left\|\mathbb E_{\nu_{\epsilon,B}}
 [\tau XY^{\mathsf T}]\right\|_F=O(\sqrt N)           \tag{SH.0d}
```

on a convex bridge set with exponentially small complement.  This is
strictly less information than the full bridge pressure landscape.  The
spectral theorem below is one sufficient certificate for (SH.0d).
Equivalently, for two independent replicas from the parent Gibbs law, its
square is the signed overlap product

```math
\left\|\mathbb E[\tau XY^{\mathsf T}]\right\|_F^2
=\mathbb E[\tau\tau'\langle X,X'\rangle
                       \langle Y,Y'\rangle].          \tag{SH.0e}
```

Thus (SH.0d) asks only for an `O(N)` two-replica response, not a full
coset/energy histogram.

More generally, inserting (SH.0c) into (SH.0b) yields the following exact
phase dichotomy.  If, for a fixed `lambda<c/h`,

```math
G_{\lambda,\epsilon}\ge\eta N,                        \tag{SH.0f}
```

then on **every** convex bridge carrier with complement at most `e^(-cN)`,

```math
\boxed{
\sup_{B\in\mathcal K_N}
 {\left\|\mathbb E_{\nu_{\epsilon,B}}
          [\tau XY^{\mathsf T}]\right\|_F^2\over N^2}
\ge {\eta+o(1)\over C_0\lambda\beta^2}.}             \tag{SH.0g}
```

Thus a real fixed-tilt phase forces a fixed mass of macroscopic two-replica
row response on every exponentially typical convex carrier.  Conversely,
the strictly weaker condition

```math
\sup_{B\in\mathcal K_N}
 \|\mathbb E[\tau XY^{\mathsf T}]\|_F=o(N)            \tag{SH.0h}
```

already rules the phase out.  This is much weaker than a uniform covariance
operator bound or the `O(sqrt(N))` estimate in (SH.0d).

## 2. Uniform physical-scale closure under a spectral margin

Let

```math
S_\epsilon(B)=
\begin{pmatrix}A&B\\B^{\mathsf T}&\epsilon D\end{pmatrix}.      \tag{SH.2}
```

**Theorem SH.1 (strict-high-temperature actual-child closure).**  Fix a
balanced window `theta_0<=m/N<=1-theta_0`, constants `s>0` and `kappa<1/2`,
and suppose the actual children obey

```math
\boxed{
t\max\{\|A\|_{\rm op},\|D\|_{\rm op}\}
+\beta\left(\sqrt{m/N}+\sqrt{n/N}+s\right)
\le\kappa.}                                           \tag{SH.3}
```

Then, uniformly in the orientation,

```math
\boxed{
\operatorname{Var}_{U_B}(L_\epsilon)=O_{\beta,\kappa,s}(1),
\qquad
\sigma_{\rm cross}^2(L_\epsilon)=O_{\beta,\kappa,s}(1).}       \tag{SH.4}
```

There is a constant `lambda_0=lambda_0(beta,kappa,s)>0` such that every fixed
`0<lambda<lambda_0` also satisfies

```math
\boxed{
G_{\lambda,\epsilon}
=\mathbb E_{U_B}L_\epsilon
 +{1\over\lambda}\log\mathbb E_{U_B}e^{-\lambda L_\epsilon}
=O_{\beta,\kappa,s,\lambda}(1).}                      \tag{SH.5}
```

Consequently both terms in the exact row-shadow decomposition (AC.24) are
`O(1)`:

```math
0\le\mathbb E_UL_\epsilon-V_\lambda^{\rm row}\le O(1),
\qquad
0\le\lambda^{-1}\mathcal I_\lambda^{\leftarrow}\le O(1).       \tag{SH.6}
```

In this regime there is no linear negative-disorder phase, and all physical
cross-row ANOVA mass is tight rather than merely its zero-amplitude tangent.

*Proof.*  A rectangular Rademacher matrix satisfies

```math
\Pr\{\|B\|_{\rm op}>
 \sqrt m+\sqrt n+s\sqrt N\}\le2e^{-c_sN}.             \tag{SH.7}
```

On the complementary event, (SH.3) and the triangle inequality give
`\|tS_epsilon(B)\|op<=kappa`.

On the convex set

```math
\mathcal K=\{B:\|tS_\epsilon(B)\|_{\rm op}\le\kappa\},         \tag{SH.8}
```

the dimension-free high-temperature covariance bound in
`artifacts/high_temperature_frobenius_pressure_stability.md` applies to both
signs of the interaction.  If `C` is either spin covariance, then
`\|C\|op<=K_kappa`.  Its `m`-by-`n` cross block has Frobenius norm at most
`K_kappa sqrt(min(m,n))`.  Differentiating the symmetrized log-partition
therefore gives

```math
\|\nabla_BL_\epsilon(B)\|_F
\le tK_\kappa\sqrt{\min(m,n)}
\le\beta K_\kappa.                                   \tag{SH.9}
```

The supremum of the supporting affine functions to `L_epsilon` on `K` is a
convex `beta K_kappa`-Lipschitz extension `g` to the whole bridge cube and
agrees with `L_epsilon` on `K`.  Convex-Lipschitz concentration on the
Rademacher cube gives, for fixed real `z`,

```math
\log\mathbb E e^{z(g-\mathbb Eg)}
\le C_{\beta,\kappa}z^2.                              \tag{SH.10}
```

Both `L_epsilon` and `g` have polynomially bounded range, so (SH.7) implies
`E|L_epsilon-g|^2=o(1)`.  Equations (SH.9)--(SH.10) prove
`Var(L_epsilon)=O(1)`, and row-ANOVA orthogonality gives the second part of
(SH.4).

It remains to guard the negative moment from the exceptional set.  Jensen's
inequality in the bridge and
`cosh(a+b)<=2cosh(a)cosh(b)` give

```math
\mathbb E_BL_\epsilon(B)
\le\log2+F_m(t)+F_n(t)+mn\log\cosh t
\le\log2+{\beta^2(N-1)\over4}.                       \tag{SH.11}
```

The last inequality uses actual child minimality and the random-sign
annealed competitor.  Since `L_epsilon>=0`, the bad-set contribution to

```math
\mathbb E e^{-\lambda(L_\epsilon-\mathbb EL_\epsilon)}
```

is at most

```math
2\exp\{-c_sN+\lambda[\log2+\beta^2(N-1)/4]\}.         \tag{SH.12}
```

Choose `lambda_0` so that `lambda_0 beta^2/4<c_s/2`.
On the good set replace `L_epsilon` by `g`, use (SH.10), and use the
`o(1)` mean difference.  The centered negative moment is between one and
`exp(O(lambda^2))+o(1)`.  Dividing its logarithm by `lambda` proves (SH.5).
Finally, both terms of (AC.24) are nonnegative and sum to `G`, proving
(SH.6). `square`

For example, if one could prove

```math
\|A\|_{\rm op}\le K\sqrt m,
\qquad \|D\|_{\rm op}\le K\sqrt n                   \tag{SH.13}
```

uniformly for the selected minimizing children, then (SH.3) holds throughout
the explicit interval

```math
\beta\left[
 K\max\{\sqrt\theta,\sqrt{1-\theta}\}
 +\sqrt\theta+\sqrt{1-\theta}+s\right]<1/2.          \tag{SH.14}
```

The theorem therefore identifies a concrete optimizer-specific route to an
order-uniform physical result.  It is only a high-temperature closure; it
does not approach the large-`beta` zero-temperature objective.

## 3. Pressure near-minimality does not supply the spectral premise

The missing spectral input cannot be obtained from asymptotic pressure
optimality alone.

**Proposition SH.2 (localized spectral spikes are pressure-invisible).**
Fix `beta>0` and let `A_n` be any exact minimizer at raw temperature
`t_n=beta/sqrt(n)`.  For every exponent

```math
{1\over2}<\alpha<{3\over4},                           \tag{SH.15}
```

there is a hollow sign matrix `A_n'` such that

```math
\boxed{
0\le
 \log\overline Z_n(A_n',t_n)-F_n(t_n)
\le 2\beta n^{2\alpha-1/2}=o(n),}                    \tag{SH.16}
```

while

```math
\boxed{\|A_n'\|_{\rm op}\ge n^\alpha-2.}             \tag{SH.17}
```

More precisely, `A_n'` is an `o(n)`-near-minimizer and has an operator norm
larger than `sqrt(n)` by a power.

*Proof.*  Choose `k=floor(n^alpha)` vertices and replace every sign in their
principal block by `+1`.  At most `k(k-1)/2` edge signs change.  Flipping one
edge changes the log-partition by at most `2t_n`, so

```math
|\log\overline Z_n(A_n',t_n)
  -\log\overline Z_n(A_n,t_n)|
\le t_nk(k-1)\le\beta n^{2\alpha-1/2}.                \tag{SH.18}
```

(The factor two in (SH.16) is a harmless convention-safe relaxation.)  The
edited principal submatrix is `J_k-I_k`, whose operator norm is `k-1`.
Operator norm cannot decrease under compression to a principal subspace, so
`\|A_n'\|op>=k-1>=n^alpha-2`. `square`

The proposition is a scalable obstruction, not a claim that the exact
minimizer itself has a spike.  It proves that any argument using only an
`o(n)` pressure-near-minimizer hypothesis, a limiting pressure value, or an
`o(n)` variational defect cannot establish (SH.13).  Exact finite optimizer
inequalities would have to be used in an essentially stronger way.

## 4. Director judgment

The physical-scale theorem is sharp about what high-temperature covariance
can currently do:

```math
\boxed{
\text{uniform spectral margin}
\Longrightarrow
\text{tight physical row ANOVA and no fixed-small-}lambda\text{ gain}.}     \tag{SH.19}
```

It does not reset the campaign's fixed-`lambda` SML, for two reasons.

1. No uniform `O(sqrt n)` operator bound is known for the actual exact
   thermal minimizers, and Proposition SH.2 rules out deriving it from
   asymptotic near-optimality alone.
2. Even an optimizer spectral theorem would close only an explicit small-
   `beta` interval.  The convergence problem ultimately needs every fixed
   `beta` before `beta` tends to infinity.

The exact next lemma on this branch would be:

> **Exact-minimizer typical-response lemma.**  Prove directly from the full
> flip-minimality inequalities that (SH.0d) holds along a convex
> exponentially typical bridge interpolation set, or exhibit a quantified
> localized response mode that forces a competing pressure-decreasing flip.

This is strictly narrower than arbitrary child-law classification, but its
large-`beta` relevance remains doubtful.
