# Sublinear cross-row information cannot create a favorable conference-pressure phase

**Status.** Task-local theorem draft.  The new ingredient is a transport
argument on the alphabet of whole rows.  It extends the constant-density
row-product theorem to arbitrary dependence between rows, provided the total
correlation is `o(r)`.  The proof is conditional only on the regular-product
bulk statement isolated below; the companion sharp-edge theorem supplies
that statement for every constant-density centrally symmetric fibre.

## 1. Setup

Let `U_r` be uniform on `{+-1}^r`, let `E_r=-E_r`, and suppose

```math
p_r=U_r(E_r)\ge p_0>0.
\tag{NP.1}
```

Put `mu_r=U_r(.|E_r)` and

```math
\Sigma_r=\mathbb E_{R\sim\mu_r}RR^T.
\tag{NP.2}
```

Let `P_r` be the spectral projection of `Sigma_r` onto eigenvalues outside
`[1-r^{-1/4},1+r^{-1/4}]`, and write `V_r=I-P_r`.  The constant-density
Fourier argument gives

```math
k_r:=\operatorname{rank}P_r=O_{p_0}(\sqrt r).
\tag{NP.3}
```

Let `q_r` be an arbitrary law on exact-sign `r by r` bridges `C_r` such
that every row marginal is `mu_r`.  Its total correlation, in nats, is

```math
\tau_r
=D(q_r\|\mu_r^{\otimes r})
=\sum_{i=1}^r H((C_r)_{i,*})-H(C_r).
\tag{NP.4}
```

For a conference signing and fixed orientation, use the pressure

```math
f_r(C)=\log\left[2^{-2r}\sum_{x,y}
\cosh\left\{{\beta\over\sqrt{2r}}
\big(H_A(x)+\epsilon H_A(y)+x^TCy\big)\right\}\right]
\tag{NP.5}
```

and put

```math
h_\beta=2\psi(\beta/\sqrt2)+{\beta^2\over4}.
\tag{NP.6}
```

The only regularity input needed in the new argument is the following
product-row statement.  If `B_r` has independent rows with law `mu_r`, then

```math
\|B_rV_r\|_{op}\le(2+o_{\Pr}(1))\sqrt r.
\tag{NP.7}
```

The companion sharp-edge theorem proves (NP.7).  Keeping it visible as a
hypothesis separates the new cross-row argument from that imported
random-matrix step.

## 2. Whole-row transport

Equip the row alphabet `E_r` with the discrete metric

```math
\rho(a,b)=\mathbf 1_{a\ne b}
\tag{NP.8}
```

and `E_r^r` with its product Hamming metric

```math
d_{\rm row}(C,B)=\#\{i:C_{i,*}\ne B_{i,*}\}.
\tag{NP.9}
```

### Lemma NP.1 (total correlation transports only few whole rows)

There is a coupling `(C_r,B_r)` with

```math
C_r\sim q_r,
\qquad
B_r\sim\mu_r^{\otimes r},
\tag{NP.10}
```

such that, with `K_r=d_row(C_r,B_r)`,

```math
\boxed{
\mathbb E K_r\le\sqrt{{r\tau_r\over2}}.}
\tag{NP.11}
```

**Proof.**  A one-Lipschitz function on `(E_r^r,d_row)` has oscillation at
most one when one row symbol is changed.  Under the product measure
`mu_r^{otimes r}`, bounded differences gives

```math
\log\mathbb E\exp\{s(\phi-\mathbb E\phi)\}
\le {s^2r\over8}.
\tag{NP.12}
```

Entropy duality followed by optimization in `s` gives

```math
\mathbb E_{q_r}\phi-\mathbb E_{\mu_r^{\otimes r}}\phi
\le\sqrt{{r\tau_r\over2}}.
\tag{NP.13}
```

Applying this also to `-phi` and using finite-space
Kantorovich--Rubinstein duality proves the same bound for the row-Hamming
`W_1` distance.  An optimal coupling exists on the finite alphabet and
gives (NP.10)--(NP.11).  `square`

The coupling has the exact requested marginals.  In particular, `B_r` is
an exact independent product sample, not merely a matrix whose individual
rows have the right laws.  The pairs `(C_{i,*},B_{i,*})` need not be
independent across `i`; no such property is used.

### Lemma NP.2 (the mismatch is low-nuclear-cost)

Under the coupling in Lemma NP.1,

```math
\boxed{
\mathbb E\|C_r-B_r\|_*
\le\sqrt2\,r\sqrt{\tau_r}.}
\tag{NP.14}
```

**Proof.**  Pointwise, `C_r-B_r` has at most `K_r` nonzero rows, so its rank
is at most `K_r`.  Exact signs give

```math
\|C_r-B_r\|_F^2\le4rK_r.
\tag{NP.15}
```

Therefore

```math
\|C_r-B_r\|_*
\le\sqrt{K_r}\|C_r-B_r\|_F
\le2K_r\sqrt r.
\tag{NP.16}
```

Take expectations and use (NP.11).  `square`

In particular,

```math
\tau_r=o(r)
\quad\Longrightarrow\quad
\mathbb E\|C_r-B_r\|_*=o(r^{3/2}).
\tag{NP.17}
```

This is stronger for the present purpose than a Frobenius coupling: the
dependent endpoint can be operator-irregular, because it will be restored
by convexity rather than placed inside the high-temperature ball.

## 3. The regular product base has the iid pressure floor

### Lemma NP.3 (projected product pressure)

Assume (NP.7).  For every fixed

```math
0<\beta<{\sqrt2\over6},
\tag{NP.18}
```

one has

```math
\mathbb E\left[
\left(h_\beta-{f_r(B_rV_r)\over r}\right)_+
\right]\longrightarrow0.
\tag{NP.19}
```

**Proof.**  Couple each row of `B_r` to a uniform Rademacher row `W_r` by
the constant-density Hamming transport.  The companion reduction proves

```math
\mathbb E\|(B_r-W_r)V_r\|_F=O_{p_0}(r^{3/4})=o(r),
\tag{NP.20}
```

and

```math
\mathbb E\|W_rP_r\|_*=O_{p_0}(r)=o(r^{3/2}).
\tag{NP.21}
```

By (NP.7) and the standard iid Bernoulli edge, `B_rV_r`, `W_rV_r`, and
`W_r` lie with probability `1-o(1)` in a common strict
high-temperature operator ball.  The audited high-temperature stability
bound compares `B_rV_r` to `W_rV_r` at cost `o(r)` by (NP.20), and compares
`W_rV_r` to `W_r` at cost `o(r)` by (NP.21).  The uniform iid conference
pressure theorem supplies the floor `h_beta r+o(r)` for `W_r`.

On the exceptional operator event the positive normalized shortfall is at
most `h_beta`, since every cosh pressure is nonnegative.  This proves
(NP.19).  `square`

## 4. Near-product row dependence cannot lower pressure

### Theorem NP.4 (quantitative total-correlation barrier)

Assume (NP.1), (NP.7), and (NP.18).  There is a deterministic sequence
`epsilon_r -> 0`, depending on the fibre sequence and on `beta,p_0`, such
that every joint exact-sign law `q_r` with row marginals `mu_r` satisfies

```math
\boxed{
\mathbb E_{q_r}\left[
\left(h_\beta-{f_r(C_r)\over r}\right)_+
\right]
\le
\epsilon_r+K_\kappa\beta
\sqrt{{\tau_r\over r}}.}
\tag{NP.22}
```

Here `kappa<1/2` is any fixed common-ball parameter chosen strictly above
`beta(3+delta)/sqrt(2)` for some fixed

```math
0<\delta<{1\over\sqrt2\beta}-3,
\tag{NP.23}
```

and `K_kappa` is the constant in the audited high-temperature covariance
bound.

Consequently,

```math
\boxed{
\tau_r=o(r)
\quad\Longrightarrow\quad
\mathbb E_{q_r}\left[
\left(h_\beta-{f_r(C_r)\over r}\right)_+
\right]\longrightarrow0.}
\tag{NP.24}
```

**Proof.**  Use Lemma NP.1 to couple `C_r` to the independent product
sample `B_r`.  Let `G_r` be the event on which

```math
\|B_rV_r\|_{op}\le(2+\delta)\sqrt r.
\tag{NP.25}
```

By (NP.7), `Pr(G_r^c)=o(1)`.  On `G_r`, the audited high-temperature
covariance bound is available at the regular base `B_rV_r`.  Convexity of
pressure along the affine line from that base to `C_r` gives the global
supporting-line inequality

```math
f_r(C_r)
\ge f_r(B_rV_r)
-{K_\kappa\beta\over\sqrt{2r}}
\|C_r-B_rV_r\|_*.
\tag{NP.26}
```

No regularity of `C_r`, or of any intermediate point, is needed.  By the
nuclear triangle inequality,

```math
\|C_r-B_rV_r\|_*
\le\|C_r-B_r\|_*+\|B_rP_r\|_*.
\tag{NP.27}
```

The same conditioned-row calculation as in the companion reduction gives

```math
\mathbb E\|B_rP_r\|_*
\le k_r\sqrt{r/p_0}=O_{p_0}(r).
\tag{NP.28}
```

It follows from (NP.14) that, after division by `r`, the expected pressure
penalty in (NP.26) is at most

```math
{K_\kappa\beta\over\sqrt2\,r^{3/2}}
\left(\sqrt2\,r\sqrt{\tau_r}+O_{p_0}(r)\right)
=K_\kappa\beta\sqrt{{\tau_r\over r}}+O_{p_0,\beta}(r^{-1/2}).
\tag{NP.29}
```

On `G_r^c`, the positive normalized shortfall is at most `h_beta`.  Combine
(NP.19), (NP.26)--(NP.29), and `Pr(G_r^c)=o(1)` to obtain (NP.22).  The
specialization (NP.24) is immediate.  `square`

The proof uses the regularity of a product sample only as a base-point
property.  It never claims that the dependent bridge is close in operator
norm or itself lies in a high-temperature ball.

## 5. Linear-information consequence

The quantitative form makes the information conclusion explicit.  If for
some `gamma>0` and all large `r` in a subsequence,

```math
\mathbb E_{q_r}\left[
\left(h_\beta-{f_r(C_r)\over r}\right)_+
\right]\ge\gamma,
\tag{NP.30}
```

then (NP.22) implies along that subsequence

```math
\boxed{
\tau_r\ge
\left({\gamma-o(1)\over K_\kappa\beta}\right)^2 r
=\Omega_{\beta,\gamma}(r).}
\tag{NP.31}
```

In particular, lowering expected pressure by `gamma r` also forces this
conclusion, because

```math
(h_\beta-z)_+\ge h_\beta-z.
\tag{NP.32}
```

Thus **linear total correlation between rows is necessary** for a positive
linear favorable pressure phase when the common row fibre has constant
density.  This matches the natural entropy scale: a latent global bit that
selects between two distinct product row laws can already cost `Theta(r)`
relative entropy, whereas identifying whole rows or forcing repeated rows
typically costs much more.

## 6. Stress tests and limitations

1. **Exact marginals are preserved.**  Kantorovich--Rubinstein gives a joint
   coupling with first marginal exactly `q_r` and second marginal exactly
   `mu_r^{otimes r}`.  It does not replace `q_r` by a rowwise-independent
   approximation.  Each row of either endpoint therefore has exactly the
   declared fibre law.

2. **High-order row events are allowed.**  The argument treats each row as
   one alphabet symbol.  Full-parity fibres, block-majority-parity fibres,
   and arbitrary higher-order Fourier constraints are covered once the
   product regular-bulk theorem (NP.7) is available.

3. **A genuinely high-order cross-row constraint does not falsify the
   theorem.**  Let `z:E_r -> {+-1}` be balanced under `mu_r`, and condition
   `mu_r^{otimes r}` on

   ```math
   \prod_{i=1}^r z(C_{i,*})=1.
   ```

   Every individual row marginal remains exactly `mu_r`, while the rows
   obey an order-`r` parity constraint.  The total correlation is only
   `log 2=o(r)`, so NP.24 applies.  This test confirms that the argument is
   not silently using bounded-order row dependence.

4. **The `o(r)` threshold is sharp for this proof architecture.**  Split
   `E_r` into two equal halves with conditional laws `nu_+` and `nu_-`, and
   take

   ```math
   q_r={1\over2}\nu_+^{\otimes r}
       +{1\over2}\nu_-^{\otimes r}.
   ```

   Each row marginal is `mu_r`, but the common latent half is retained by
   all rows and

   ```math
   D(q_r\|\mu_r^{\otimes r})=(r-1)\log2.
   ```

   Thus a single reusable global latent bit already sits exactly at the
   `Theta(r)` boundary.  At this scale (NP.22) permits a constant normalized
   pressure change; the theorem neither constructs nor rules out a
   favorable linearly correlated law.

5. **Sublinear total correlation is not entrywise proximity.**  The optimal
   coupling may replace all `r` signs in each of `o(r)` rows.  What matters
   is that this difference has `o(r)` row rank and hence
   `o(r^(3/2))` nuclear cost.

6. **Product regularity remains an essential input.**  Without (NP.7), the
   support bound cannot be invoked at `B_rV_r`.  The companion
   constant-density sharp-edge theorem is what makes the present conclusion
   unconditional for the declared fibre class.  If that imported theorem
   were unavailable, NP.4 would remain a conditional extension theorem.

7. **The temperature endpoint is excluded.**  The proof requires a strict
   common high-temperature ball, so it covers fixed
   `beta<sqrt(2)/6`, not equality or larger inverse temperature.

8. **The result is one-sided.**  It excludes a lower pressure phase.  It
   does not assert convergence of the full pressure rate and does not
   exclude a dependent law from raising pressure.

9. **Vanishing row density is not covered.**  If `p_r` tends to zero, the
   Fourier rank, transport, and nuclear estimates change scale and the
   conclusion need not follow from the present proof.

## 7. Frontier movement

The product restriction is not the true boundary of the no-gain theorem.
After the regular product base is established, arbitrary cross-row
dependence of **sublinear total correlation** is harmless: it can be coupled
to that base by changing `o(r)` whole rows, and the changes are restored at
subcritical nuclear cost.

The remaining information-theoretic frontier is therefore at speed `r`.
A favorable constant-density bridge law must carry `Omega(r)` total
correlation between its rows.  The next discriminating question is whether
some linearly correlated law can genuinely exploit this allowance, or
whether a stronger response-specific transport inequality raises the
necessary information scale beyond linear.
