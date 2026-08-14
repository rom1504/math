# A finite-rank obstruction to eigenvalue-only finite-temperature lower bounds

Status: **rigorous no-go for the two natural finite-`n` bounds**.  A single
active edge, padded by isolated coordinates, disproves both the rms-only
`n psi(c_n)` lower bound and the spherical/determinant spectral lower bound.
The loss is only `O(1)`, so this does **not** disprove an `o(n)`-error theorem
for dense fixed-modulus signings.

## 1. Setup

For a symmetric zero-diagonal matrix `X`, write

```math
F_n(X)=\log \mathbb E_{\varepsilon\in\{-1,1\}^n}
 \cosh\!\left({1\over2}\varepsilon^{\mathsf T}X\varepsilon\right)
```

and

```math
\psi(c)={1\over4}\left[
 \sqrt{1+4c^2}-1-
 \log\!\left({1+\sqrt{1+4c^2}\over2}\right)
 \right].
```

Fix `0<a<1/2` and, for every `n>=2`, take

```math
X_n=a(e_1e_2^{\mathsf T}+e_2e_1^{\mathsf T}).                 \tag{1.1}
```

Then

```math
\|X_n\|_{\rm op}=a<1/2,
\qquad
\operatorname{spec}(X_n)=\{-a,a,0^{(n-2)}\},
\qquad
\operatorname{Tr}X_n^2=2a^2.                                \tag{1.2}
```

The Boolean Hamiltonian is just `a epsilon_1 epsilon_2`, and hence

```math
\boxed{F_n(X_n)=\log\cosh a.}                                \tag{1.3}
```

## 2. The rms Bernoulli bound is false

Put

```math
c_n=\sqrt{{1\over n}\operatorname{Tr}X_n^2}
    =a\sqrt{2/n}.
```

Since

```math
\psi(c)={c^2\over4}-{c^4\over8}+O(c^6),
```

we have

```math
n\psi(c_n)={a^2\over2}-{a^4\over2n}+O_a(n^{-2})
   \longrightarrow {a^2\over2}.                             \tag{2.1}
```

But

```math
\log\cosh a<{a^2\over2}\qquad(a\ne0),                       \tag{2.2}
```

because the derivative of `a^2/2-log cosh(a)` is
`a-tanh(a)>0`.  Equations (1.3)--(2.2) prove that, for every fixed
`0<a<1/2` and all sufficiently large `n`,

```math
\boxed{F_n(X_n)<n\psi\!\left(
 \sqrt{\operatorname{Tr}X_n^2/n}\right).}                   \tag{2.3}
```

For a completely explicit instance, `a=1/4` already violates (2.3) at
`n=6`.  The limiting deficit is

```math
{a^2\over2}-\log\cosh a>0
```

(approximately `0.00032019638` at `a=1/4`).

Thus variance density plus the operator-norm hypothesis cannot give an exact
finite-`n` Bernoulli-value lower bound.  Padding is the obstruction: the rms
formula Gaussianizes a fixed active core, whereas its exact negative fourth
cumulant remains visible.

## 3. The spherical/determinant relaxation also points the wrong way

Let `U_n` be uniform on the Euclidean sphere of radius `sqrt(n)`.  The
spherical counterpart of (1.3) is

```math
S_n(X_n)=\log\mathbb E
 \cosh\!\left({1\over2}U_n^{\mathsf T}X_nU_n\right)
=\log\mathbb E\cosh(aU_{n,1}U_{n,2}).                        \tag{3.1}
```

The first two coordinates converge jointly to independent standard Gaussians
`G_1,G_2`.  For `a<1/2`, the bound

```math
|aU_{n,1}U_{n,2}|
 \le {a\over2}(U_{n,1}^2+U_{n,2}^2)
```

gives uniform exponential integrability, so convergence also holds for the
expectation in (3.1).  Conditioning on `G_1` gives

```math
\mathbb E\cosh(aG_1G_2)
=\mathbb E\exp(a^2G_1^2/2)
=(1-a^2)^{-1/2}.
```

Consequently

```math
S_n(X_n)\longrightarrow-{1\over2}\log(1-a^2).                \tag{3.2}
```

This is strictly *larger* than the Boolean pressure:

```math
-{1\over2}\log(1-a^2)-\log\cosh a>0,                        \tag{3.3}
```

since the derivative of the left side is
`a/(1-a^2)-tanh(a)>a-tanh(a)>0`.  At `a=1/4`, the gap is
approximately `0.00133945695`.

In particular, the attractive comparison

```math
F_n(X)\stackrel{?}{\ge}S_n(X)
```

is false even under `diag(X)=0` and `||X||_op<1/2`.

The same example disposes of the direct determinant saddle

```math
K_n(X)=\inf_{q>\|X\|_{\rm op}}
 \left\{{n(q-1)\over2}
 -{1\over4}\log\det(q^2I-X^2)\right\}.                       \tag{3.4}
```

Indeed, for (1.1),

```math
K_n(X_n)=\inf_{q>a}\left\{
 {n(q-1)\over2}-{n-2\over2}\log q
 -{1\over2}\log(q^2-a^2)\right\}.                           \tag{3.5}
```

Its minimizer tends to `q=1`, and therefore

```math
K_n(X_n)\longrightarrow-{1\over2}\log(1-a^2)
>F_n(X_n).                                                    \tag{3.6}
```

Thus neither a Hubbard--Stratonovich argument nor a spherical saddle may be
turned into the desired lower bound by simply replacing the Boolean measure by
its rotationally invariant counterpart.  The replacement increases the
pressure on a localized eigenspace.

## 4. Exact false step in the shifted Gaussian route

For any `gamma>||X||_op`, Hubbard--Stratonovich gives the valid identity

```math
2^{-n}\sum_\varepsilon
 \exp\!\left({1\over2}\varepsilon^{\mathsf T}X\varepsilon\right)
=e^{-\gamma n/2}\,
 \mathbb E_{G\sim N(0,\gamma I+X)}\prod_i\cosh G_i.           \tag{4.1}
```

What fails is the proposed spectral lower estimate on the correlated Gaussian
product.  On the padded edge (1.1), such an estimate reduces asymptotically to
the Gaussian product value `(1-a^2)^(-1/2)`, while the left side of (4.1) is
exactly `cosh(a)`.  Inequality (3.3) has the opposite sign.  Equivalently,
there is no useful global positive-quadratic pointwise minorant of `log cosh`:
`log cosh(t)` grows only linearly as `|t|` tends to infinity.

## 5. Scope and the viable weakened target

This counterexample establishes:

1. an exact bound `F_n(X)>=n psi(sqrt(Tr X^2/n))` is false for general hollow
   high-temperature matrices;
2. the finite spherical/free-probability determinant functional is not a
   universal Boolean lower bound;
3. eigenvalues do not by themselves justify Gaussianizing a localized active
   subspace.

It does **not** establish any of the following:

- failure of `F_n(X)>=n psi(c_n)-o(n)`;
- failure for dense fixed-modulus signing matrices;
- failure after adding a localization-sensitive `O(1)` or `o(n)` correction;
- failure of the asymptotic theorem under the existing traffic/delocalization
  assumptions.

Any finite-`n` repair must pay at least

```math
{a^2\over2}-\log\cosh a
```

on the padded-edge family, or explicitly exclude finite-rank localization.
Thus a pressure-preserving regularization/delocalization theorem remains the
right missing bridge; a bare Gaussian-mixture, determinant, or eigenvalue-only
substitution does not remove it.

## 6. The isotropic Gaussian trial for polar involutions has a fixed gap

There is a second exact boundary relevant to the polar target.  Suppose
`U=2P-I` is a symmetric involution with `diag(U)=0`.  Then `P=VV^T` is a
rank-`n/2` projection and every row `v_i` of `V` has squared norm `1/2`.
Gaussian linearization gives

```math
\mathbb E_x\exp\!\left({\beta\over2}x^{\mathsf T}Ux\right)
=e^{-\beta n/2}\,
 \mathbb E_g\prod_{i=1}^n
 \cosh\!\left(\sqrt{2\beta}\,v_i^{\mathsf T}g\right).       \tag{6.1}
```

Use `N(0,sI_(n/2))` as a trial law in the Gibbs variational formula on the
right.  Since every marginal `v_i^Tg` has variance `s/2`, this proves the
frame-independent lower bound

```math
{1\over n}\log\mathbb E_x
 \exp\!\left({\beta\over2}x^{\mathsf T}Ux\right)
\ge \mathcal L(\beta)
:=\sup_{s>0}\left\{-{\beta\over2}
 +\mathbb E\log\cosh(\sqrt{\beta s}\,G)
 -{s-1-\log s\over4}\right\},                             \tag{6.2}
```

where `G` is standard Gaussian.  Thus the isotropic trial is a valid joint
inequality, not a scalar-channel decomposition.  The same bound applies to
`-U`, and hence to normalized `cosh` pressure.  It nevertheless misses the
conference value by a fixed amount.  Expanding

```math
\mathbb E\log\cosh(\sqrt aG)
={a\over2}-{a^2\over4}+{a^3\over3}-{17a^4\over24}+O(a^5)
```

and solving the saddle equation gives `s=1+2 beta+2 beta^2+O(beta^4)` and

```math
\mathcal L(\beta)
={\beta^2\over4}-{5\beta^4\over24}+O(\beta^6),
\qquad
\psi(\beta)-\mathcal L(\beta)
={\beta^4\over12}+O(\beta^6).                              \tag{6.3}
```

Hence an isotropic Gaussian/Brascamp--Lieb trial cannot prove the polar
involution target.  A successful variational law must retain genuinely
frame-dependent information or use an additional dense-flat stability
argument.
