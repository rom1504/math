# High-temperature pressure is stable in normalized Frobenius distance

Status: **proved pressure-continuity lemma**.  In the strict operator-norm
high-temperature regime, normalized `cosh` pressure is Lipschitz in nuclear
norm and therefore continuous in normalized Frobenius norm.  This permits
distributed edits touching linearly many coordinates; it is strictly more
flexible than the now-falsified exceptional-coordinate puncturing proposal.

## 1. Statement

For a real symmetric `n` by `n` interaction matrix `X`, put

```math
Z_n(X)=2^{-n}\sum_{x\in\{-1,1\}^n}
 \exp\!\left({1\over2}x^{\mathsf T}Xx\right),
\qquad
\overline Z_n(X)={Z_n(X)+Z_n(-X)\over2}.
\tag{1.1}
```

Fix `kappa<1/2`.  There is a finite constant `K_kappa`, independent of
dimension, such that whenever

```math
\|X\|_{\rm op},\|Y\|_{\rm op}\le\kappa,
\tag{1.2}
```

one has

```math
\boxed{
\left|\log\overline Z_n(X)-\log\overline Z_n(Y)\right|
\le {K_\kappa\over2}\|X-Y\|_* .}
\tag{1.3}
```

Here `||.||_*` is nuclear norm.  Since
`||E||_*<=sqrt(n)||E||_F`, (1.3) implies

```math
{1\over n}\left|\log\overline Z_n(X_n)
 -\log\overline Z_n(Y_n)\right|\longrightarrow0
\tag{1.4}
```

for any two uniformly strict-high-temperature sequences satisfying

```math
{1\over\sqrt n}\|X_n-Y_n\|_{\rm F}\longrightarrow0.
\tag{1.5}
```

## 2. Dimension-free covariance input

Bauerschmidt and Bodineau,
[*A very simple proof of the LSI for high temperature spin systems*](https://arxiv.org/abs/1712.03676),
prove a dimension-free functional inequality for Ising interactions after a
diagonal shift makes the interaction positive semidefinite.  In the form used
explicitly in the proof of Corollary 2.10 of
Fan--Misiakiewicz--Wang--Wen,
[*Dynamical mean-field limit and replica-symmetric free energy for the
orthogonally-invariant SK model*](https://arxiv.org/abs/2607.10102), it gives:
for every `kappa<1/2`, every symmetric `J` with `||J||_op<=kappa`, and every
external field, the Ising law with interaction `J` satisfies a Poincare
inequality whose constant depends only on `kappa` (and on the field bound).

At zero field, global spin flip makes every one-spin mean zero.  Thus, if

```math
C_J=\mathbb E_J[xx^{\mathsf T}],
```

the Poincare inequality applied to `v^T x` gives the uniform covariance
bound

```math
v^{\mathsf T}C_Jv=\operatorname{Var}_J(v^{\mathsf T}x)
\le K_\kappa\|v\|_2^2,
\qquad\hbox{hence}\qquad
\|C_J\|_{\rm op}\le K_\kappa.                         \tag{2.1}
```

The strict margin in (1.2) is essential to keeping `K_kappa` independent of
`n`; no claim is made at the boundary `1/2`.

## 3. Interpolation proof

Let `E=X-Y` and `J_s=Y+sE`.  Convexity of operator norm gives
`||J_s||_op<=kappa` for `0<=s<=1`.  Differentiating the ordinary normalized
partition function yields

```math
{d\over ds}\log Z_n(J_s)
={1\over2}\mathbb E_{J_s}[x^{\mathsf T}Ex]
={1\over2}\operatorname{Tr}(E C_{J_s}).                 \tag{3.1}
```

Nuclear/operator duality and (2.1) give

```math
\left|{d\over ds}\log Z_n(J_s)\right|
\le {K_\kappa\over2}\|E\|_* .                          \tag{3.2}
```

Integration proves the same bound for `log Z_n(X)-log Z_n(Y)`.  Apply it
again to `-X,-Y`.  If two positive summands each change by at most a
multiplicative factor `exp(delta)`, their sum does as well; using (1.1)
therefore proves (1.3).  Finally, nuclear norm is at most `sqrt(n)` times
Frobenius norm, proving (1.4)--(1.5).

## 4. Research consequence

The correct high-temperature regularization metric is substantially weaker
than maximum-entry power control.  To transfer a known pressure limit from a
regular sequence `Y_n` to a possibly traffic-irregular sequence `X_n`, it is
enough to find `Y_n` with the same strict operator-norm margin and normalized
Frobenius distance `o(1)`.  The edit may touch a linear number of coordinates
and may create linearly many exceptional power entries.

This does not by itself construct the regular endpoint.  In particular,
suppose `beta>0` and

```math
{1\over\sqrt n}\|X_n^2-\beta^2I\|_{\rm F}\longrightarrow0.
\tag{4.1}
```

Functional calculus with

```math
\big||\lambda|-\beta\big|
\le {1\over\beta}|\lambda^2-\beta^2|
```

shows that the spectral polar matrix

```math
U_n=\beta\,\operatorname{sgn}(X_n),
\qquad U_n^2=\beta^2I,                                  \tag{4.2}
```

may be chosen so that `||X_n-U_n||_F/sqrt(n)->0`.  When both sequences stay
strictly below the `1/2` operator threshold, (1.4) proves that they have the
same limiting pressure.  Thus Frobenius-near conference structure reduces
the pressure problem to a symmetric involution.

That involution need not satisfy the entrywise delocalization hypotheses of
the available free-energy theorem; eigenvalues alone still do not determine
Boolean pressure.  The new positive target is therefore a **global
Frobenius regularization** to an endpoint with a verified pressure law, or a
direct pressure theorem for the polar involutions arising from flat hollow
sign matrices, rather than deletion of `o(n)` bad coordinates.
