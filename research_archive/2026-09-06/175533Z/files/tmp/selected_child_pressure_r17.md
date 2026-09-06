# Wave 17 Route 2: selected-child pressure and the exact `A_9` wall

## Status

- **Verified:** the zero-temperature one-deletion excess identity, the
  finite-temperature cavity factorization, both centered drift identities,
  the conditional constant-shortfall comparison theorem, and the
  fixed-orientation Finner inequality below.
- **Verified finite obstruction:** for the exact order-nine minimizer, every
  one-vertex cavity reward is strictly less than `3` at every inverse
  temperature `beta>0`.  This simultaneously defeats the simple
  `kappa>=q_r-q_{r-1}` sufficient condition and the convergence-relevant
  root-scale-centered zero-error criterion.
- **Numerical:** the cavity identity and selected normalized pressure were
  evaluated on all named `A_2,A_4,A_5,A_6,A_7,A_8,A_9` examples and every
  principal restriction of `A_9`.  A bounded additive shortfall is compatible
  with all of these data, but no asymptotic theorem is proved.
- **Open target:** prove a uniform constant-shortfall cavity selection lemma
  on restrictions of global minimizers.  This would imply convergence.

## 1. Exact zero-temperature deletion algebra

Let `B` have order `r`, let `C_i=B[-i]`, and put

```math
d_i=Q(B)-Q(C_i),
\qquad
\varepsilon(B)=Q(B)-q_r.
```

Then, exactly,

```math
\boxed{
\varepsilon(C_i)
=\varepsilon(B)+(q_r-q_{r-1})-d_i.
}
\tag{SC1}
```

Consequently

```math
\boxed{
e(C_i)-e(B)
=\frac{q_r-q_{r-1}-d_i}{(r-1)^{3/2}}
+\varepsilon(B)
\left[(r-1)^{-3/2}-r^{-3/2}\right].
}
\tag{SC2}
```

Thus selection by terminal excess is exactly selection by a large decrement;
there is no hidden second scalar.

There is a useful degeneracy refinement.  Let `g(C_i)` be the number of
oriented projective absolute grounds of `C_i`, and let `N_B(t)` count oriented
projective states of parent deficit at most `t`.  If `b_i` is the deleted row,
then a child ground `(sigma,y)` has two parent extensions with deficits

```math
d_i\mp2\sigma b_i^{\mathsf T}y.
```

Parent optimality gives `2|b_i^T y|<=d_i`; the extension map is injective.
Hence

```math
\boxed{N_B(2d_i)\ge2g(C_i).}
\tag{SC3}
```

In particular, a flat step `d_i=0` forces `g(B)>=2g(C_i)`.  For the certified
`A_9`, `g(A_9)=25`, all nine children have `d_i=0`, and their ground counts
are `2` or `4`.  For the child obtained by deleting zero-based vertex seven,
`Q=24=q_8+4`, `g=4`, and every order-seven child has
`Q=22=q_7+4`, recovering the flat `A_9 -> A_8 -> A_7` wall.  The normalized
terminal excess rises on the second step even though the raw excess stays
four.

## 2. Exact finite-temperature cavity identity

Let `Omega_r` be the `2^r` oriented projective states `(sigma,x)`, and set

```math
e_B(\sigma,x)=\sigma x^{\mathsf T}Bx,
\qquad
\mathcal Z_\beta(B)
=2^{-r}\sum_{\omega\in\Omega_r}e^{\beta e_B(\omega)},
\qquad
P_\beta(B)=\beta^{-1}\log\mathcal Z_\beta(B).
```

The elementary pressure bounds are

```math
Q(B)-\frac{r\log2}{\beta}
\le P_\beta(B)\le Q(B).
\tag{SC4}
```

For a fixed deletion `i`, let `mu_{beta,C_i}` be the Gibbs law proportional
to `exp(beta e_{C_i}(sigma,y))`, and define its cavity reward

```math
\kappa_{\beta,i}(B)
=\frac1\beta\log
\mathbb E_{\mu_{\beta,C_i}}
\cosh\!\left(2\beta b_i^{\mathsf T}y\right).
```

Every child state has two parent extensions.  Summing those extensions gives
the exact factorization

```math
\boxed{
\mathcal Z_\beta(B)
=\mathcal Z_\beta(C_i)
\mathbb E_{\mu_{\beta,C_i}}
\cosh\!\left(2\beta b_i^{\mathsf T}y\right),
\qquad
P_\beta(B)-P_\beta(C_i)=\kappa_{\beta,i}(B)\ge0.
}
\tag{SC5}
```

Moreover `kappa_{beta,i}->d_i` as `beta->infinity`.  Thus (SC5) is exactly a
soft deletion decrement, not an average statement.

Centering by the finite minimum gives

```math
\boxed{
[P_\beta(C_i)-q_{r-1}]-[P_\beta(B)-q_r]
=(q_r-q_{r-1})-\kappa_{\beta,i}(B).
}
\tag{SC6}
```

If

```math
D_\beta(B)=\sum_{\omega\in\Omega_r}
e^{-\beta[Q(B)-e_B(\omega)]},
```

then the normalized centered pressure is

```math
\boxed{
\frac{P_\beta(B)-q_r}{r^{3/2}}
=e(B)+\frac{\log D_\beta(B)}{\beta r^{3/2}}
-\frac{\log2}{\beta\sqrt r}.
}
\tag{SC7}
```

This is a precise strengthening of (10.531): terminal excess plus a positive
near-ground entropy correction, minus the vanishing order potential
`V_beta(r)=log(2)/(beta sqrt(r))`.  A deterministic child makes (SC7)
nonincreasing exactly when

```math
\kappa_{\beta,i}(B)
\ge q_r-q_{r-1}
+[P_\beta(B)-q_r]
\left[1-\left(1-\frac1r\right)^{3/2}\right].
\tag{SC8}
```

The simpler condition `kappa>=q_r-q_{r-1}` is sufficient along a path on
which the normalized centered pressure starts nonpositive and never rises.
At a deterministic landing order `N`, (SC7) would then give
`e(B_N)<=log(2)/(beta sqrt(N))`.  This still controls only terminal excess,
so by itself it does not imply convergence.

There is no nontrivial parent-average identity hidden here:

```math
\frac1r\sum_i\kappa_{\beta,i}(B)
=P_\beta(B)-\frac1r\sum_iP_\beta(C_i).
```

A deterministic selected child uses `max_i kappa_i`; a randomized child uses
the corresponding average.  Neither has a known universal lower bound at the
required scale.

## 3. Exact `A_9` obstruction at every temperature

For `A_9`, all children have `Q=24`.  Put `s=e^{-beta}` and use the
unnormalized deficit polynomial

```math
Z_B(s)=\sum_{\omega}s^{Q(B)-e_B(\omega)}.
```

Exact enumeration gives

```math
Z_{A_9}(s)
=25+60s^8+111s^{16}+120s^{24}
+111s^{32}+60s^{40}+25s^{48}.
\tag{SC9}
```

There are only three child-polynomial types.  For each of the nine children,
the exact checker converts

```math
2Z_{C_i}(s)-s^3Z_{A_9}(s)
```

to its degree-51 Bernstein basis and finds every coefficient nonnegative,
with the zeroth coefficient strictly positive.  Therefore it is strictly
positive for `0<=s<1`.  Since the normalized state counts are `512` and
`256`, this proves

```math
\boxed{
\kappa_{\beta,i}(A_9)<3
\quad\text{for every }i\text{ and every }\beta>0.
}
\tag{SC10}
```

This is stronger than a numerical temperature scan.  Because
`q_9-q_8=4`, (SC10) falsifies the error-free sufficient criterion after
(SC8) for every fixed `beta`.  It also falsifies the more useful root-scale
zero-error criterion below: for

```math
\alpha=\frac{q_9}{9^{3/2}}=\frac89,
\qquad
\alpha[9^{3/2}-8^{3/2}]
=24-\frac{128}{9}\sqrt2>3.
```

Thus finite temperature does not dissolve the flat wall; it only turns it
into a bounded cavity shortfall.

## 4. The corrected root-scale pressure target

Center instead at the root scale `alpha`:

```math
H_{\alpha,\beta}(B)
=P_\beta(B)-\alpha r^{3/2}.
```

Equation (SC5) gives the exact drift

```math
\boxed{
H_{\alpha,\beta}(C_i)-H_{\alpha,\beta}(B)
=\alpha[r^{3/2}-(r-1)^{3/2}]
-\kappa_{\beta,i}(B).
}
\tag{SC11}
```

This centering does address the desired scalar comparison.  Let `A` be an
order-`n` global minimizer, set `alpha=q_n/n^{3/2}`, and suppose one can
select a deletion path down to order `m` such that, at every order `r`,

```math
\boxed{
\kappa_{\beta,i}(B_r)
\ge\alpha[r^{3/2}-(r-1)^{3/2}]-K
}
\tag{SC12}
```

for fixed `beta>0` and a universal constant `K`.  Telescoping (SC11), using
`P_beta(A)<=Q(A)=alpha n^{3/2}`, and then (SC4), proves

```math
\boxed{
\frac{q_m}{m^{3/2}}
\le\frac{q_n}{n^{3/2}}
+\frac{K(n-m)}{m^{3/2}}
+\frac{\log2}{\beta\sqrt m}.
}
\tag{SC13}
```

For `m>=rho n`, the error is `O_{rho,K,beta}(m^{-1/2})`.  Such comparisons
are tail-summable over geometric order windows and would prove convergence
via (10.529).  Thus (SC12), unlike terminal excess alone, is a genuine
sufficient lemma.

The `A_9` certificate proves that `K=0` is false at every temperature and in
fact forces `K>24-(128/9)sqrt(2)-3`.  It does not obstruct a universal
constant `K`.  The finite census is consistent with bounded `K`, but that is
only numerical evidence.

A fixed positive `beta` is the clean useful scale: the entropy error in
(SC13) is already `O(m^{-1/2})`, so no zero-temperature limit is needed.
Taking `beta=b/sqrt(r)` with fixed `b` leaves a nonvanishing entropy error;
taking `b->infinity` changes the temperature between parent and child and
introduces an uncontrolled pressure drift unless beta is frozen over each
comparison window.  In addition, the zero-temperature limit returns the
flat decrement wall.  Fixed beta avoids both issues.

## 5. Finner/Shearer check

For a fixed orientation `sigma`, let

```math
Z_B^\sigma(\beta)
=\mathbb E_xe^{\beta\sigma x^{\mathsf T}Bx}.
```

Since

```math
\sum_i x_{-i}^{\mathsf T}C_ix_{-i}
=(r-2)x^{\mathsf T}Bx,
```

Finner's inequality on the complements of single coordinates gives

```math
\boxed{
Z_B^\sigma(\beta)
\le
\prod_{i=1}^r
Z_{C_i}^\sigma\!\left(\beta\frac{r-1}{r-2}\right)^{1/(r-1)}.
}
\tag{SC14}
```

This forces at least one child pressure to be large; selected descent needs
one child pressure to be small, equivalently a large cavity reward.  It is
therefore the wrong direction.  The orientation variable cannot simply be
included in the same fractional cover because it belongs to all `r` factors,
whose weights sum to `r/(r-1)>1`.  SC14 supplies no two-sided selected-child
bound.

## 6. Checkers

- `tmp/selected_child_r17.py`: exact restriction norms, excesses, ground
  counts, and the named flat chain.
- `tmp/check_selected_child_soft_r17.py`: all-energy cavity-factorization
  audit on the named minimizers and every restriction of `A_9`; floating
  point is used only for log-sum-exp values.
- `tmp/prove_a9_soft_wall_r17.py`: exact integer deficit polynomials and
  rational Bernstein coefficients proving (SC10).

The sharp ledger-ready conclusion is: **the exact soft cavity identity is
verified; the simple `q`-gap condition and the convergence-relevant
root-scale zero-error condition are both falsified by `A_9` at every fixed
temperature; the constant-shortfall root-scale criterion (SC12) is a new
open target which, if proved, would yield the required adaptive comparison
tail and convergence.**
