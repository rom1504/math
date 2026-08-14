# Spectral extremality and the exact finite-temperature regularity dichotomy

Status: **proved reduction**.  This note removes the symmetry assumption from
the spectral extremality argument and identifies the exact remaining escape
route for a signing whose fixed-temperature `cosh` pressure is strictly below
the conference value.  It does not prove that pressure minimizers are regular.

## 1. Spectral extremality without a symmetry assumption

Let `mu` be a compactly supported probability law with

```math
\mathbb E X=0,
\qquad \mathbb E X^2=\sigma^2>0.
\tag{1.1}
```

Assume that the real inverse Cauchy-transform branches of `mu` and `-mu`
both extend from zero through `u=1`.  In particular this holds when
`supp(mu)` is contained in `(-1/2,1/2)`, which is the regime used by the
high-temperature theorem.  Define

```math
J(\mu)={1\over2}\int_0^1 R_\mu(u)\,du
\tag{1.2}
```

and

```math
\psi(\sigma)={1\over4}\left[
 \sqrt{1+4\sigma^2}-1
 -\log\left({1+\sqrt{1+4\sigma^2}\over2}\right)
 \right].
\tag{1.3}
```

**Lemma 1.1 (two-orientation spectral extremality).**  Under (1.1)--(1.2),

```math
\boxed{J(\mu)+J(-\mu)\ge2\psi(\sigma),
\qquad \max\{J(\mu),J(-\mu)\}\ge\psi(\sigma).}
\tag{1.4}
```

Equality in the second inequality holds if and only if

```math
\mu={1\over2}(\delta_{-\sigma}+\delta_\sigma).
\tag{1.5}
```

Thus the symmetric Bernoulli law uniquely minimizes the two-orientation
free energy; symmetry of the competing law is not an assumption.

**Proof.**  Let `G_mu` be the Cauchy transform and let

```math
G_\mu(z)=1,
\qquad G_{-\mu}(y)=1,
\qquad
z>\sup\operatorname{supp}\mu,
\quad y>-\inf\operatorname{supp}\mu.
\tag{1.6}
```

These are the values of the two inverse branches at `u=1`.  Integration of
`R_mu(u)=G_mu^{-1}(u)-1/u` gives the exact identities

```math
2J(\mu)=z-1-\mathbb E\log(z-X),
\qquad
2J(-\mu)=y-1-\mathbb E\log(y+X).
\tag{1.7}
```

For example, (1.7) follows by differentiating

```math
uG_\mu^{-1}(u)
-\mathbb E\log(G_\mu^{-1}(u)-X)-\log u
\tag{1.8}
```

and observing that its limit at `u=0` is one.

The product `(z-X)(y+X)` is positive.  Jensen's inequality and (1.1) give

```math
\begin{aligned}
2\{J(\mu)+J(-\mu)\}
&=z+y-2-\mathbb E\log\{(z-X)(y+X)\}\\
&\ge z+y-2-\log(zy-\sigma^2).
\end{aligned}
\tag{1.9}
```

Put `q=sqrt(zy)`.  Positivity in (1.9) implies `q>sigma`, and AM--GM yields

```math
2\{J(\mu)+J(-\mu)\}
\ge 2q-2-\log(q^2-\sigma^2).
\tag{1.10}
```

The right side has its unique minimum at

```math
q_\sigma={1+\sqrt{1+4\sigma^2}\over2},
\qquad q_\sigma^2-\sigma^2=q_\sigma,
\tag{1.11}
```

and the minimum is `4 psi(sigma)`.  This proves (1.4).

If `max{J(mu),J(-mu)}=psi(sigma)`, then their sum is at most
`2 psi(sigma)`, so equality holds throughout (1.9)--(1.10).  Equality in
AM--GM gives `z=y`, while strict equality in Jensen says that

```math
(z-X)(z+X)=z^2-X^2
\tag{1.12}
```

is almost surely constant.  Hence `X^2=sigma^2` almost surely, and the
mean-zero condition forces the two masses to be equal.  Conversely (1.5)
has `J(mu)=J(-mu)=psi(sigma)`.  This proves the equality statement.

There is also an exact nonnegative penalty decomposition.  With

```math
\mathcal E=
 \log(zy-\sigma^2)
 -\mathbb E\log\{(z-X)(y+X)\}\ge0,
\qquad
h(q)=2q-2-\log(q^2-\sigma^2),
\tag{1.13}
```

one has

```math
\max\{J(\mu),J(-\mu)\}-\psi(\sigma)
\ge {1\over4}\left[
 (\sqrt z-\sqrt y)^2+h(q)-h(q_\sigma)+\mathcal E
 \right].
\tag{1.14}
```

Thus spectral asymmetry, displacement of the inverse-transform saddle, and
nonconstant squared spectral magnitude are all paid before any entrywise
regularity question arises.

## 2. Consequence for signing sequences

For a zero-diagonal symmetric signing `A_n`, put

```math
X_n={\beta A_n\over\sqrt n},
\qquad
p_n(A_n,\beta)={1\over n}\log\left[
 2^{-n}\sum_{x\in\{-1,1\}^n}
 \cosh\left({1\over2}x^{\mathsf T}X_nx\right)
 \right].
\tag{2.1}
```

For a fixed positive integer `k`, define the power-delocalization residual

```math
\Delta_{n,k}=
 \max_i\left|(X_n^k)_{ii}-{1\over n}\operatorname{Tr}X_n^k\right|
 +\max_{i\ne j}|(X_n^k)_{ij}|.
\tag{2.2}
```

**Corollary 2.1 (exact regularity dichotomy).**  Fix `beta>0` and a signing
subsequence `A_n`.  Suppose

```math
\limsup_n {\beta\|A_n\|_{\rm op}\over\sqrt n}< {1\over2}
\tag{2.3}
```

and, for every fixed `k` and every `eta>0`,

```math
\Delta_{n,k}<n^{-1/2+\eta}
\quad\hbox{eventually}.
\tag{2.4}
```

Then

```math
\boxed{\liminf_n p_n(A_n,\beta)\ge\psi(\beta).}
\tag{2.5}
```

Consequently, if for some `delta>0`

```math
p_n(A_n,\beta)\le\psi(\beta)-\delta
\tag{2.6}
```

along an infinite subsequence, then on a further infinite subsequence either

```math
\limsup_n {\beta\|A_n\|_{\rm op}\over\sqrt n}\ge {1\over2},
\tag{2.7}
```

or there are a fixed `k` and a fixed `eta>0` for which

```math
\Delta_{n,k}\ge n^{-1/2+\eta}
\tag{2.8}
```

infinitely often.

**Proof.**  Under (2.3), every subsequence has a further subsequence whose
empirical spectral law converges to a compact law `mu` supported strictly
inside `(-1/2,1/2)`.  The signing identities give

```math
{1\over n}\operatorname{Tr}X_n=0,
\qquad
{1\over n}\operatorname{Tr}X_n^2
=\beta^2(1-1/n),
\tag{2.9}
```

so `mu` has mean zero and variance `beta^2`.

Condition (2.4) is exactly the fixed-power hypothesis in Assumption 2.9 of
the deterministic high-temperature theorem of Fan--Misiakiewicz--Wang--Wen.
It applies to both `X_n` and `-X_n`, giving

```math
{1\over n}\log Z(X_n)\longrightarrow\log2+J(\mu),
\qquad
{1\over n}\log Z(-X_n)\longrightarrow\log2+J(-\mu).
\tag{2.10}
```

Since the normalized `cosh` partition function is

```math
{Z(X_n)+Z(-X_n)\over2^{n+1}},
\tag{2.11}
```

its pressure converges along this further subsequence to
`max{J(mu),J(-mu)}`, which is at least `psi(beta)` by Lemma 1.1.  Applying
this argument to a subsequence realizing the liminf proves (2.5).
The final statement is the contrapositive; failure of (2.4) has the explicit
form (2.8).

## 3. What the lemma does and does not settle

Corollary 2.1 is the requested exact bridge for every uniformly
high-temperature, fixed-power-delocalized sequence.  It proves that a
strictly subconference pressure sequence cannot remain in that class.  No
symmetry of its limiting spectrum is needed.

The threshold alternative cannot presently be sharpened from
`limsup >= 1/2` to an eventual pointwise inequality: a sequence may approach
the boundary from below, where the imported theorem has no uniform margin.
Nor does the result prove that exact pressure minimizers satisfy (2.3) or
(2.4).  Those two possibilities are now the only escapes, rather than
unexamined regularity assumptions.

A large operator norm by itself need not create a linear pressure penalty.
Indeed, changing only the edges inside a set of `c sqrt(n)` vertices changes
`log Zbar` by at most `O(beta sqrt(n))=o(n)`, while making that principal
submatrix all positive forces an eigenvalue at least `c sqrt(n)-1`.  Thus a
localized spectral spike can cross the high-temperature threshold without
altering the limiting pressure.  This explains why (2.7) is a genuine escape
case and why a direct norm-only penalty cannot close it.
