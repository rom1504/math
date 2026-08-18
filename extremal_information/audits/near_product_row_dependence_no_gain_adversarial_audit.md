# Adversarial audit: near-product cross-row dependence no-gain theorem

**Frozen source:**
`extremal_information/drafts/near_product_row_dependence_no_gain.md`

**SHA-256:**
`ac9475dee0090a8f98bcdb9bedae7c2731c5e67f75f5619294bd6eec65c170a8`

**Verdict:** **PASS.**  The row-level transport constant, total-correlation
identity, nuclear-rank estimate, supporting-line direction and gradient
normalization, exceptional-event split, and linear-information conclusion
are correct.  In particular, no mismatch or pressure-comparison cost is
needed on the irregular event.  No repair is required.

This audit treats the companion projected sharp-edge statement NP.7 as a
proved input, as requested.

## 1. Total correlation and the exact row-transport constant

Because every row marginal of `q_r` is `mu_r`, direct expansion on the
finite alphabet gives

```math
D(q_r\|\mu_r^{\otimes r})
=-H(q_r)-E_{q_r}\log\prod_i\mu_r(C_i)
=\sum_iH((C_r)_{i,*})-H(C_r).
```

Thus NP.4 is an identity in nats, not an inequality or an independence
assumption.

For a function with row-Hamming Lipschitz constant one, changing the `i`th
row changes its value by at most one.  The bounded-difference exponential
lemma therefore gives exactly

```math
\log E_{\mu^{\otimes r}}e^{s(\phi-E\phi)}
\le {s^2\over8}\sum_{i=1}^r1^2={s^2r\over8}.
```

Entropy duality yields

```math
E_q\phi-E_{\mu^{\otimes r}}\phi
\le {\tau_r\over s}+{sr\over8}.
```

The optimizer is `s=sqrt(8tau_r/r)` (with the zero-entropy case obtained by
continuity), and the optimized value is

```math
\sqrt{r\tau_r/2}.
```

Finite-space Kantorovich--Rubinstein duality for the row-Hamming cost turns
this into the `W_1` distance and supplies an optimal coupling with the exact
two stated marginals.  Applying the estimate to `-phi` is harmless, though
the one-sided dual formula already suffices.  NP.11 has no missing factor of
two.

## 2. Nuclear cost of whole-row mismatches

If `K` whole rows differ, `D=C-B` has rank at most `K` and at most `rK`
nonzero entries, each of magnitude two.  Hence, pointwise,

```math
\|D\|_F^2\le4rK,
\qquad
\|D\|_*\le\sqrt{\operatorname{rank}D}\,\|D\|_F
\le2K\sqrt r.
```

Taking expectations and substituting NP.11 gives

```math
E\|C-B\|_*
\le2\sqrt r\sqrt{r\tau_r/2}
=\sqrt2\,r\sqrt{\tau_r}.
```

Thus `tau_r=o(r)` gives `o(r^(3/2))` expected nuclear cost exactly as
claimed.  Dependence among the paired row mismatches is irrelevant.

## 3. Exceptional covariance projection cost

The estimate NP.28 can be reconstructed directly.  Since conditioning a
uniform Rademacher row on an event of probability at least `p_0` gives the
quadratic-form domination

```math
\Sigma_r\preceq p_0^{-1}I,
```

one has

```math
E\|B_rP_r\|_F^2
=r\operatorname{tr}(P_r\Sigma_r)
\le {rk_r\over p_0}.
```

Because `rank(BP)<=k_r`,

```math
E\|B_rP_r\|_*
\le\sqrt{k_r}\,(E\|B_rP_r\|_F^2)^{1/2}
\le k_r\sqrt{r/p_0}=O_{p_0}(r).
```

The nuclear-triangle split NP.27 is therefore valid and introduces only a
normalized `O(r^(-1/2))` pressure penalty.

## 4. Supporting line: direction and gradient normalization

On

```math
G_r=\{\|B_rV_r\|_{op}\le(2+\delta)\sqrt r\},
```

the full scaled interaction with bridge `B_rV_r` has norm at most

```math
{\beta\over\sqrt2}
 (\sqrt{1-1/r}+2+\delta)<\kappa<1/2
```

for all large `r`.  Let

```math
D=C_r-B_rV_r,
\qquad
g(s)=f_r(B_rV_r+sD).
```

The pressure is a log-sum-exp of affine functions of the bridge, so `g` is
globally convex.  At the regular base only, the audited covariance bound
gives

```math
|g'(0)|
\le {K_\kappa\over2}
 \left\|{\beta\over\sqrt{2r}}
 \begin{pmatrix}0&D\\D^T&0\end{pmatrix}\right\|_*
={K_\kappa\beta\over\sqrt{2r}}\|D\|_*.
```

The equality uses the two copies of every singular value in the symmetric
dilation.  Convexity then gives, in the correct direction,

```math
f_r(C_r)=g(1)\ge g(0)+g'(0)
\ge f_r(B_rV_r)-{K_\kappa\beta\over\sqrt{2r}}\|D\|_*.
```

Neither the dependent endpoint nor any intermediate point needs an
operator bound.  NP.26 has the correct sign, scale, and factor.

## 5. The irregular event requires no mismatch payment

This is the most delicate logical point.  Put

```math
Z_C=(h_\beta-f_r(C_r)/r)_+,
\qquad
Z_B=(h_\beta-f_r(B_rV_r)/r)_+,
```

and let

```math
L={K_\kappa\beta\over\sqrt2\,r^{3/2}}
  \|C_r-B_rV_r\|_*.
```

On `G_r`, NP.26 and `(a+b)_+<=a_++b` for `b>=0` give

```math
Z_C\le Z_B+L.
```

On `G_r^c`, pressure nonnegativity gives simply

```math
0\le Z_C\le h_\beta.
```

Therefore the required expectation inequality is

```math
E Z_C
\le E Z_B+E L+h_\beta P(G_r^c).
\tag{EA.1}
```

The nuclear penalty in (EA.1) may safely be bounded *without* its indicator;
there is no need to compare `C_r` and `B_rV_r` or pay their mismatch on
`G_r^c`.  The first term tends to zero by NP.19 and the last by NP.7.
Substituting NP.14 and NP.28 in the middle term gives exactly

```math
K_\kappa\beta\sqrt{\tau_r/r}
+O_{p_0,\beta}(r^{-1/2}).
```

This proves NP.22 with a deterministic `epsilon_r->0` depending only on the
fixed fibre sequence, `beta`, and `p_0`, not on the joint law `q_r`.

## 6. Quantitative information inference

If the left side of NP.22 is at least a fixed `gamma>0`, then for all large
orders on that subsequence

```math
K_\kappa\beta\sqrt{\tau_r/r}\ge\gamma-\epsilon_r,
```

and squaring gives

```math
\tau_r\ge
\left({\gamma-o(1)\over K_\kappa\beta}\right)^2r.
```

Thus the `Omega(r)` conclusion has the advertised dependence and no hidden
loss.  If instead `E f_r(C_r)/r<=h_beta-gamma`, the pointwise inequality

```math
(h_\beta-z)_+\ge h_\beta-z
```

puts the shortfall expectation above `gamma`, so the same inference applies.

The two stress-test entropy calculations are also consistent: conditioning
on one balanced parity costs `log 2`, while the mixture of two disjoint
product halves has density `2^(r-1)` on either monochromatic half and hence
total correlation `(r-1)log 2`.

## 7. Precise repair

None required.
