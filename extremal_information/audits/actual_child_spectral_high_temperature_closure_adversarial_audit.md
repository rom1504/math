# Adversarial audit: spectral high-temperature closure

**Object audited:**
[`drafts/actual_child_spectral_high_temperature_closure.md`](../drafts/actual_child_spectral_high_temperature_closure.md)

**Verdict:** **PASS.**  The spectral good event, covariance-to-Frobenius
gradient estimate, convex extension, exceptional negative-moment estimate,
and localized-spike construction have the stated factors.  No theorem
modification is required.  The result remains deliberately confined to a
strict small-`beta` interval and does not establish the unproved spectral
premise for exact minimizers.

## 0. Typical-gradient abstraction and phase dichotomy

**Verdict on SH.0: PASS.**  In the actual log-partition setting the global
bound `|partial_(ij)L|<=t` makes the supporting-plane extension polynomially
bounded even if the displayed hypotheses do not separately bound `R_N`.
Thus the exponentially small exceptional set gives
`E|L-g|^2=o(1)`.  The assumptions `L>=0` and `E L<=hN` are exactly what is
needed for

```math
E[e^{-\lambda(L-EL)};K_N^c]
\le e^{-(c-\lambda h)N};
```

the strict restriction `lambda<c/h` is therefore correct.  The good-set MGF
gives `log E exp[-lambda(L-EL)]<=C_0lambda^2R_N^2+o(N)`, and division by
`lambda` proves (SH.0b).

The response conversion has the right scale:

```math
R_N^2={\beta^2\over N}
 \sup_{B\in K_N}\|E_{\nu_B}[\tau XY^T]\|_F^2.
```

Hence a gain at least `eta N` forces (SH.0g) on each carrier separately.
The quantifier “every carrier” is valid provided the same fixed `c,h` and
`lambda<c/h` are used, as stated.  Conversely response norm `o(N)` gives
`R_N=o(sqrt(N))` and excludes linear gain.  Finally, (SH.0e) follows by
expanding the Frobenius square with two independent replicas; there is no
missing absolute value or orientation factor because the resulting signed
expectation equals a squared norm and is automatically nonnegative.

## 1. Spectral event and gradient bound

For the block matrix in (SH.2),

```math
\|S_\epsilon(B)\|_{op}
\le \max(\|A\|_{op},\|D\|_{op})+\|B\|_{op}.
```

The off-diagonal block operator has norm exactly `||B||op`, so (SH.3)
combined with the standard rectangular Rademacher tail gives (SH.7)--(SH.8)
with an exponential `e^{-c_sN}` exceptional probability.

For either Gibbs sign, the high-temperature covariance theorem gives a full
spin covariance `C` with `||C||op<=K_kappa`.  Its rectangular cross block
therefore obeys

```math
\|C_{12}\|_F\le\sqrt{\min(m,n)}\,\|C_{12}\|_{op}
\le K_\kappa\sqrt{\min(m,n)}.
```

Differentiating the symmetric block interaction contributes exactly `t`,
not `t/2` or `2t`, because the two symmetric off-diagonal entries cancel the
quadratic-form factor `1/2`.  Taking the convex mixture of the two Gibbs signs
preserves the bound, proving (SH.9).

## 2. Convex extension and two-sided concentration

`L_epsilon` is a log-sum-exp of affine bridge functionals and hence convex.
The set `K` is convex because it is an operator-norm sublevel set of an
affine map.  The supremum of tangent affine minorants based at points of `K`
therefore agrees with `L_epsilon` on `K` and has Euclidean Lipschitz constant
at most `beta K_kappa`.  Standard convex-Lipschitz concentration on a
Rademacher product cube is two-sided about a median (and hence yields a
centered subgaussian MGF after changing universal constants), so (SH.10) is
valid for positive and negative real `z`.

Both `L_epsilon` and this extension have polynomial range on the finite
cube.  Multiplying the squared range by `e^{-c_sN}` proves
`E|L_epsilon-g|^2=o(1)`, which is enough to transfer the variance bound and,
by ANOVA orthogonality, the cross-row bound.

## 3. Exceptional negative moment

The annealed child competitor gives

```math
F_m(t)+F_n(t)+mn\log\cosh t
\le {N(N-1)\over2}\log\cosh t
\le {\beta^2(N-1)\over4}.
```

The harmless extra `log 2` covers the symmetrized orientation split, so
(SH.11) is valid.  Since every symmetrized normalized partition function is
at least one, `L_epsilon>=0`.  Consequently the bad set contributes at most
the quantity in (SH.12).  Choosing
`lambda beta^2/4<c_s/2` makes it exponentially negligible.  On the good set,
the negative MGF of `g-Eg` is bounded by (SH.10), and the `o(1)` difference
of means is harmless.  Jensen supplies the lower bound one.  Thus (SH.5)
and the nonnegative split (SH.6) follow.

## 4. Localized spike

Editing a `k`-vertex principal block changes at most `k(k-1)/2` signs; the
`2t_n` one-edge oscillation therefore gives `t_nk(k-1)`, matching (SH.18).
For `1/2<alpha<3/4` this is `o(n)`.  The edited principal compression is
`J_k-I_k`, of spectral norm `k-1`; compression cannot increase operator
norm, hence the full matrix has norm at least `k-1`.  SH.2 is therefore a
valid scalable ceiling for arguments based only on `o(n)` pressure
near-optimality.  It does not claim, and does not imply, that an exact
minimizer itself contains such a spike.
