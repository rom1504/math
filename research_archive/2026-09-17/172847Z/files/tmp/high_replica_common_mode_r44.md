# Wave 44 memo: high replicas cannot bypass the inverse restriction gap

Status: the stitching identity and pressure upper bound below are
**Verified**. The finite data are exact integer-deficit enumerations followed
by floating exponentiation. No tracked file was edited.

## 1. Exact common-center partition

Fix

```math
d_S(z)=Q(A[S])-|z_S^{\mathsf T}A[S]z_S|,\qquad
w_S(y)=e^{-\theta d_S(y)},\qquad
q(z)=\mathbb E_{S\sim U_m}w_S(z_S).
```

Use full-cube spins here. All quantities are invariant under global negation,
so full-cube and projective expectations agree. Define

```math
\zeta_S=2^{-m}\sum_yw_S(y),\qquad
\mu_S(y)=\frac{w_S(y)}{2^m\zeta_S}.
\tag{R44.1}
```

For a fixed tuple `boldsymbol S=(S_1,...,S_r)`, let
`U=bigcup_a S_a` and independently sample `Y_a~mu_(S_a)`. Let `Cons` mean
that the local words agree exactly on every overlap. On `Cons`, stitch them
on `U` and extend by uniform spins outside `U`, obtaining `Z`. For the old
hard row law `nu_2=U_z(.|mathcal C_2)`, direct summation gives

```math
\boxed{
\mathbb E_{z\sim\nu_2}\prod_{a=1}^r w_{S_a}(z_{S_a})
=\frac{2^{rm-|U|}\prod_a\zeta_{S_a}}{U_z(\mathcal C_2)}
\Pr\{\mathrm{Cons},\ R_2(Z)\le2n(n-1)\}.
}
\tag{R44.2}
```

Consequently

```math
\boxed{
\mathbb E_{z\sim\nu_2}q(z)^r
=\mathbb E_{S_1,\ldots,S_r}
\left[
\frac{2^{rm-|\cup_aS_a|}\prod_a\zeta_{S_a}}{U_z(\mathcal C_2)}
\Pr\{\mathrm{Cons},\ R_2(Z)\le2n(n-1)\}
\right].
}
\tag{R44.3}
```

This exposes a genuine soft-label common-mode and its hard row correction.
It is an exact normalization, not by itself a lower bound.

## 2. Centered signed energy is subproject at a low-row center

Put

```math
E(z)=z^{\mathsf T}Az,\qquad
C_S(z)=z_S^{\mathsf T}A[S]z_S,\qquad
X_S(z)=C_S(z)-p_2E(z).
```

Integrating the slice tail (10.1124) gives, uniformly for
`|lambda|<=c_0/||A||_op`,

```math
\log\mathbb E_{U_m}e^{\lambda X_S(z)}
\le C\left[
\log n+|\lambda|\frac{q_n}{n}
+\lambda^2\{R_2(z)+n^2\}
\right].
\tag{R44.4}
```

The `log n` is the harmless slice-conditioning cost. At

```math
\theta=\frac{bL_0}{H}=\Theta(n^{-3/4-c}),\qquad
L_0=n^{3/4-c},\qquad0<c<1/4,
```

one has `theta||A||_op=O(n^(-c))`. Hence

```math
\sup_{z\in\mathcal C_2,\ \sigma=\pm1}
\log\mathbb E e^{\theta\sigma X_S(z)}
=O(\log n+n^{1/2-2c})=o(L_0).
\tag{R44.5}
```

For the weaker project cap `R_2(z)=O(n^(9/4-c))`, the same bound is

```math
O(\log n+n^{3/4-3c})=o(L_0).
\tag{R44.6}
```

Thus (10.1152) uses the old hard class, but the obstruction below also holds
for every cut admitted by the larger project cap in (10.795).

## 3. Strict pressure contains a stronger restriction edge

Every principal signing satisfies `Q(A[S])>=q_m`, and
`e^(theta|x|)<=e^(theta x)+e^(-theta x)`. Therefore (R44.5)--(R44.6) give,
uniformly for every eligible center,

```math
\begin{aligned}
q(z)
&=\mathbb E_Se^{-\theta Q(A[S])+\theta|C_S(z)|}\\
&\le e^{-\theta q_m}\sum_{\sigma=\pm1}
e^{\theta\sigma p_2E(z)}\mathbb E_Se^{\theta\sigma X_S(z)}\\
&\le2\exp\{-\theta q_m+\theta p_2q_n+o(L_0)\}.
\end{aligned}
\tag{R44.7}
```

Recall

```math
H=(p^{3/2}-p_2)q_n+t,\qquad
G_{n,m}=q_m-p^{3/2}q_n,\qquad
K_z^{\rm abs}(\theta)=\theta H+\log q(z).
```

Substitution proves the central **Verified inverse theorem**:

```math
\boxed{
K_z^{\rm abs}(\theta)\le\theta\{t-G_{n,m}\}+o(L_0).
}
\tag{R44.8}
```

It holds uniformly on either row class. Consequently, for every `r>=1`,
including `r=Theta(n/L_0)`, the high-replica pressure obeys

```math
\boxed{
\mathcal P_r(b)\le\theta\{t-G_{n,m}\}+o(L_0).
}
\tag{R44.9}
```

The last step merely transfers the new pointwise bound through a power mean;
it is not the substantive max approximation already recorded in (10.1152).

If `mathcal P_r(b)>=aL_0` for fixed `a>0`, then
`theta=bL_0/H` and (R44.9) already force

```math
\boxed{
G_{n,m}\le t-\frac abH+o(H).
}
\tag{R44.10}
```

This is a constant-leading negative restriction gap, much stronger than the
power-saving upper edge sought from (10.795). In the near-equality regime
`G_(n,m)=t+o(H)`, every eligible center has only
`K_z^abs=o(L_0)`, so no fixed positive margin is possible.

Therefore high-replica strict positive pressure cannot independently create
the desired restriction edge: its premise already contains the stronger
conclusion (R44.10). This retires (10.1152) as a proof mechanism. It does
**not** falsify the hard exceptional-center tail (10.967) or the bare tail
(10.795), because fixed positive pressure was only a sufficient soft-to-hard
condition, never a necessary one.

## 4. Optional principal-pressure diagnosis

Let

```math
Z_Q(\theta)=\mathbb E_Se^{-\theta Q(A[S])},\qquad
\Pi_Q(\theta)=\log\mathbb E_S
e^{-\theta(Q(A[S])-\mathbb E Q(A[S]))}.
```

Hölder with exponents `1+epsilon` and `(1+epsilon)/epsilon`, taking
`epsilon=Cn^(-c)`, combines (R44.4) with `Q>=0` to give

```math
K_z^{\rm abs}(\theta)
\le\theta\{t-\overline G_{n,m}\}+\Pi_Q(\theta)+o(L_0),
\tag{R44.11}
```

where
`overline G_(n,m)=E_SQ(A[S])-p^(3/2)q_n` is the circular mean excess in
(10.970). The old-C2 Hölder error is
`O(n^(3/4-2c)+n^(1/2-c)+log n)` and remains `o(L_0)` at the project cap.
Since `Pi_Q<=theta(E Q-q_m)`, this recovers (R44.8).

This identifies the only possible scalar residual as principal-norm lower
tail. The self-bound estimate (10.1135) controls its mgf only at
`O(n^(1-2c))`, larger than `L_0`; (R44.11) is an inverse diagnosis, not a
retry of ordinary moments.

## 5. Finite independent checks

`high_replica_common_mode_r44_check.py` checks (R44.2) on a nontrivial
three-selector tuple, checks `E_(nu_2)q^3` by both sides of the tuple
expansion, and verifies the exact finite Hölder inequality separately for
both signs and every center.

At `t=0` and finite normalization `theta H=1`:

| signing, `m` | old `C_2` projective centers | `max R_2` | `max K_z` | `G_(n,m)` | `Pi_Q` |
|:---|---:|---:|---:|---:|---:|
| `A_6,5` | `32/32` | 30 | 0.817719 | 0.392742 | 0 |
| `A_8,6` | `128/128` | 72 | 0.203325 | 1.581048 | 0.125249 |
| `A_9,7` | `256/256` | 128 | 0.041747 | 4.537547 | 0.296423 |

These small orders are not in the asymptotic regime where the error in
(R44.8) is negligible; their positive finite pressures do not contradict the
theorem.

For the hard event `d_S(z)<=H`, exact `A_8,m=6` has 22 empty selector-pair
intersections out of 378 and 1700 empty triple intersections out of 3276.
Exact `A_9,m=7` has 362/630 empty pairs and 6070/7140 empty triples. Yet both
have positive finite best-center pressure. This falsifies a uniform hard
pairwise-Helly or all-list compatibility premise, not the soft partition
(R44.3), an averaged overlap theorem, or an asymptotic hard tail.

## 6. Frontier recommendation

Retire high-replica strict pressure (10.1152). Preserve the hard
exceptional-center tail (10.967) and bare tail (10.795). Any exceptional-
center continuation must work directly at the hard boundary or use a soft
functional without a fixed positive pressure margin; replicas cannot evade
the inverse restriction gap.
