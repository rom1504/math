# Bernoulli extremizes the symmetrized high-temperature spectral pressure

Status: **proved free-probability inequality and exact pressure mapping**.  Among compact
mean-zero spectral laws of fixed variance, the symmetric Bernoulli law uniquely minimizes the
larger of the two replica-symmetric free energies associated with `mu` and `-mu`.  Through
Corollary 2.10(c) of Fan--Misiakiewicz--Wang--Wen, this is exactly the spectral quantity selected
by the normalized `cosh` pressure.

The theorem gives an unconditional all-order upper bound for the minimized pressure by
restricting asymptotically nearby Paley conference matrices.  A matching lower bound still needs
a pressure-preserving regularization theorem for exact or near minimizers.  That remaining
obligation is stated precisely in Section 5.

The three logical scopes in this note are independent and should not be conflated:

1. **Exact spectral inequality (unconditional at the level of measures).**  Theorem 1.1 proves
   `max{J(mu),J(-mu)} >= psi(c)` for every admissible compact mean-zero law of variance `c^2`;
   it neither assumes nor proves that a signing sequence has the replica-symmetric pressure
   predicted by its spectral law.
2. **Pressure limsup (unconditional for signings).**  Near-order Paley restrictions prove
   `limsup P_n(beta)/n <= psi(beta)` for each fixed `0<beta<1/2`, without any hypothesis on
   pressure minimizers.
3. **Pressure liminf (conditional).**  The reverse inequality requires the explicit
   pressure-preserving regularization target in Section 5.2.  Even if that target is proved, the
   present argument applies only at fixed `0<beta<1/2`; it does not by itself reach the
   zero-temperature regime needed for convergence of `M_n/n^(3/2)`.

## 1. The spectral theorem

Let `mu` be a compactly supported probability law, let `X` have law `mu`, and assume

```math
\mathbb EX=0,
\qquad \mathbb EX^2=c^2.
\tag{1.1}
```

Assume that the exterior Cauchy-transform inverse defining `R_mu(t)` and `R_(-mu)(t)` is valid
for `0<=t<=1`.  It is enough, and is the only case used below, that both supports lie strictly
inside `(-1,1)`.  Define

```math
J(\mu)={1\over2}\int_0^1R_\mu(t)\,dt
\tag{1.2}
```

and

```math
\psi(c)={1\over4}\left[
 \sqrt{1+4c^2}-1
 -\log\!\left({1+\sqrt{1+4c^2}\over2}\right)
 \right].
\tag{1.3}
```

**Theorem 1.1 (symmetrized Bernoulli extremizer).**

```math
\boxed{\max\{J(\mu),J(-\mu)\}\ge\psi(c).}
\tag{1.4}
```

If `c>0`, equality holds if and only if

```math
\mu={1\over2}(\delta_{-c}+\delta_c).
\tag{1.5}
```

For `c=0`, equality means the degenerate law `delta_0`, the natural zero-variance version of
(1.5).

## 2. Variational identity, including the domain check

Write

```math
G_\mu(z)=\mathbb E{1\over z-X},
\qquad b=\sup\operatorname{supp}\mu.
```

On the exterior branch, let `K_mu(t)=G_mu^(-1)(t)` and
`R_mu(t)=K_mu(t)-1/t`.  An antiderivative calculation gives

```math
\begin{aligned}
\int_0^1R_\mu(t)\,dt
&=K_\mu(1)-1-\mathbb E\log(K_\mu(1)-X)\\
&=\inf_{z>b}\left\{z-1-\mathbb E\log(z-X)\right\}.
\end{aligned}
\tag{2.1}
```

For completeness, if `z=K_mu(t)`, then

```math
{d\over dt}\{tz-\mathbb E\log(z-X)-\log t\}
=z-{1\over t}=R_\mu(t).
```

As `t` decreases to zero, the expression in braces tends to `1`, proving the first line of
(2.1).  The derivative of the variational functional is `1-G_mu(z)`, so its unique stationary
point is `K_mu(1)` and is its minimum.

There is no hidden boundary issue in the high-temperature application.  If
`supp(mu) subset (-1,1)`, Jensen's inequality for the strictly convex map
`x -> 1/(1-x)` gives

```math
G_\mu(1)=\mathbb E{1\over1-X}>1
\tag{2.2}
```

when `c>0`; meanwhile `G_mu(z)` decreases to zero as `z` tends to infinity.  Hence
`K_mu(1)>1>b` is an interior minimizer.  The same argument applies to `-mu`.  More generally,
(2.1) remains true with an infimum and a limiting exterior point whenever the inverse branch is
valid up to `t=1`.

Thus

```math
2J(\mu)=\inf_{z>b}F_+(z),
\qquad
2J(-\mu)=\inf_{w>-\inf\operatorname{supp}\mu}F_-(w),
\tag{2.3}
```

where

```math
F_+(z)=z-1-\mathbb E\log(z-X),
\qquad
F_-(w)=w-1-\mathbb E\log(w+X).
```

## 3. Product Jensen proof and equality

For any admissible positive `z,w`, put `p=zw`.  Both factors below are pointwise positive, and

```math
\mathbb E[(z-X)(w+X)]
=zw+(z-w)\mathbb EX-\mathbb EX^2
=p-c^2>0.
\tag{3.1}
```

Concavity of `log`, followed by arithmetic--geometric mean, gives

```math
\begin{aligned}
F_+(z)+F_-(w)
&=z+w-2-\mathbb E\log[(z-X)(w+X)]\\
&\ge z+w-2-\log(p-c^2)\\
&\ge 2\sqrt p-2-\log(p-c^2).
\end{aligned}
\tag{3.2}
```

This proof does not require the variational infima to be attained.  Taking both infima in (3.2)
and using `max(a,b)>=(a+b)/2` yields

```math
\max\{J(\mu),J(-\mu)\}
\ge {1\over4}\inf_{p>c^2}
 \{2\sqrt p-2-\log(p-c^2)\}.
\tag{3.3}
```

Put `q=sqrt(p)`.  The unique minimizer satisfies

```math
q^2-q-c^2=0,
\qquad
q={1+\sqrt{1+4c^2}\over2},
\qquad
p-c^2=q.
\tag{3.4}
```

Substitution in (3.3) gives exactly (1.3)--(1.4).

Suppose now that `c>0` and equality holds.  Equality must hold in `max >= average`, in AM--GM,
and in Jensen.  Hence

```math
J(\mu)=J(-\mu),
\qquad z=w=q,
```

and `(q-X)(q+X)=q^2-X^2` is constant almost surely.  Therefore `X^2=c^2` almost surely.
The mean-zero condition forces equal masses at `-c` and `c`, proving (1.5).  Conversely, direct
inversion of

```math
G(z)={z\over z^2-c^2}
```

gives

```math
R(t)={\sqrt{1+4c^2t^2}-1\over2t},
```

whose integral is `psi(c)`, so the symmetric Bernoulli law attains equality.

## 4. Exact mapping to normalized `cosh` pressure

Let `X_n` be a deterministic symmetric sequence satisfying Assumption 2.9 of
Fan--Misiakiewicz--Wang--Wen,
[“Dynamical mean-field limit and replica-symmetric free energy for the
orthogonally-invariant SK model”](https://arxiv.org/abs/2607.10102), with limiting spectral law
`mu` and an eventual operator-norm bound strictly below `1/2`.  Corollary 2.10(c), at zero
external field, gives

```math
{1\over n}\log Z(X_n)\longrightarrow\log2+J(\mu),
\qquad
{1\over n}\log Z(-X_n)\longrightarrow\log2+J(-\mu),
\tag{4.1}
```

where

```math
Z(X)=\sum_{x\in\{-1,1\}^n}\exp\!\left({1\over2}x^{\mathsf T}Xx\right).
```

The normalized `cosh` partition is exactly

```math
\overline Z(X)
=2^{-n}\sum_x\cosh\!\left({1\over2}x^{\mathsf T}Xx\right)
=2^{-n-1}\{Z(X)+Z(-X)\}.
\tag{4.2}
```

Therefore the logarithmic sum selects the larger orientation:

```math
\boxed{
{1\over n}\log\overline Z(X_n)
\longrightarrow\max\{J(\mu),J(-\mu)\}
\ge\psi(c).}
\tag{4.3}
```

For a signing `A_n` at scaled temperature `beta`, take

```math
X_n={\beta\over\sqrt n}A_n.
```

Zero diagonal gives mean spectral value zero, while flat off-diagonal signs give

```math
{1\over n}\operatorname{Tr}X_n^2
=\beta^2\left(1-{1\over n}\right)\longrightarrow\beta^2.
\tag{4.4}
```

Thus every signing sequence to which Corollary 2.10(c) applies obeys

```math
\liminf {1\over n}\log\overline Z(X_n)\ge\psi(\beta),
\tag{4.5}
```

with equality at the spectral level only for the conference law
`(delta_(-beta)+delta_beta)/2`.

## 5. What is proved for the minimized pressure, and what remains

Write

```math
P_n(\beta)=\min_{A_n}
 \log\overline Z\!\left({\beta A_n\over\sqrt n}\right).
\tag{5.1}
```

### 5.1 Unconditional all-order upper bound

For every fixed `0<beta<1/2`,

```math
\boxed{\limsup_{n\to\infty}{P_n(\beta)\over n}\le\psi(\beta).}
\tag{5.2}
```

Indeed, the prime number theorem in the progression `1 mod 4` supplies Paley conference orders
`r=r(n)>=n` with `r/n -> 1`.  Take an `n` by `n` principal restriction `C_n` of the order-`r`
conference matrix.  At a fixed raw temperature `u`, adding a signed vertex can only increase the
normalized `cosh` partition:

```math
\mathbb E_y\cosh(u\{H(x)+yL(x)\})
=\cosh(uH(x))\cosh(uL(x))
\ge\cosh(uH(x)).
\tag{5.3}
```

With `u=beta/sqrt(n)`, this gives

```math
P_n(\beta)
\le\log\overline Z_n(C_n,u)
\le\log\overline Z_r(A_r,u).
\tag{5.4}
```

The scaled full conference matrix `uA_r` has limiting spectral law
`(delta_(-beta)+delta_beta)/2`, operator norm tending to `beta<1/2`, and satisfies Assumption
2.9 exactly through `A_r^2=(r-1)I`.  Equation (4.3) therefore makes the last term in (5.4)
`r psi(beta)+o(r)=n psi(beta)+o(n)`, proving (5.2).

### 5.2 Exact remaining regularization lemma

The missing lower bound is not another spectral optimization: Theorem 1.1 has completed that
part.  It is the following signing regularity statement.

**Pressure-preserving regularization target.**  For each fixed `0<beta<1/2` and
each exact pressure minimizer `A_n^*`, construct a signing `A_tilde_n` such that

```math
\log\overline Z\!\left({\beta\widetilde A_n\over\sqrt n}\right)
\le P_n(\beta)+o_\beta(n),
\tag{5.5}
```

and, for some `eta_beta>0`,

```math
\left\|{\beta\widetilde A_n\over\sqrt n}\right\|_{\rm op}
\le {1\over2}-\eta_\beta,
\tag{5.6}
```

while the sequence satisfies the fixed-power diagonal/off-diagonal delocalization condition
(2.1).  It is enough to establish these properties directly for `A_n^*`, taking
`A_tilde_n=A_n^*`.

To see sufficiency, pass from any subsequence to a further subsequence on which the empirical
spectral laws converge.  Equations (4.4), (5.6), and Theorem 1.1 give

```math
{1\over n}\log\overline Z
 \!\left({\beta\widetilde A_n\over\sqrt n}\right)
\ge\psi(\beta)-o(1).
```

Equation (5.5) then yields `P_n(beta)/n >= psi(beta)-o(1)`.  Together with (5.2), this would prove

```math
{P_n(\beta)\over n}\longrightarrow\psi(\beta)
\tag{5.7}
```

throughout the temperature range where the regularization holds.

This obligation is substantial and currently open.  Neither pressure minimization nor the flat
edge magnitudes currently imply the operator-norm or power-delocalization hypotheses.  The result
also remains intrinsically high-temperature: convergence for `beta<1/2` alone is insufficient for
the zero-temperature squeeze, which requires fixed `beta` arbitrarily large before sending
`beta` to infinity.
