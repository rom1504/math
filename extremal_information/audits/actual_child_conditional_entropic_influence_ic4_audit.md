# Audit of IC.4: conditional-entropic influence

**Object.**  Theorem IC.4 in
[`../drafts/actual_child_interaction_curvature_dichotomy.md`](../drafts/actual_child_interaction_curvature_dichotomy.md).

**Verdict: PASS.**  Every sign and normalization is correct.  The
`s`-dependent conditional log-partition is a legitimate comparison because
the product-space entropy inequality is applied at each fixed `s`; it is
never differentiated.  In fact there is an exact nonnegative gap identity
that proves `E_s<=A_s` directly.

## 1. Conditional normalization and KL direction

Fix `s>0`, a row `i`, and `B_-i`.  Put

```math
M_i(B_{-i})=E_{r_i}e^{-s h(R_i,B_{-i})}
           =e^{-s h_{i,s}(B_{-i})}.                 \tag{IC4A.1}
```

Since `q_s` has density `e^{-sh}/E_r e^{-sh}` relative to the product law
`r`, conditioning gives exactly

```math
{dq_s^{i\mid-i}\over dr_i}
={e^{-sh}\over M_i}
=e^{-s(h-h_{i,s})}.                                 \tag{IC4A.2}
```

Thus the sign in IC.27 is negative.  If `a=s(h-h_(i,s))`, then

```math
E_{q_s^{i\mid-i}}e^a=1,
```

and therefore

```math
\begin{aligned}
E_{q_s^{i\mid-i}}\tau(a)
 &=1-E_{q_s^{i\mid-i}}a-1\\
 &=-sE_{q_s^{i\mid-i}}(h-h_{i,s})\\
 &=D(q_s^{i\mid-i}\Vert r_i).                     \tag{IC4A.3}
\end{aligned}
```

This verifies both the positive sign inside `tau` and the forward KL
direction in IC.28.

## 2. Why the `s`-dependent comparison is legitimate

For a positive random variable `Y` and any positive comparison constant
`c`, the conditional entropy variational inequality is

```math
\operatorname{Ent}(Y)
\le E\{Y\log(Y/c)-Y+c\}.                            \tag{IC4A.4}
```

Its equality choice is `c=EY`.  Conditional on `B_-i`, take

```math
Y=e^{-sh},
\qquad c=M_i=e^{-s h_{i,s}}.
```

The comparison is independent of row `i`, as required, and it attains
equality in (IC4A.4).  Entropy tensorization on the product law `r` hence
gives

```math
\operatorname{Ent}_r(e^{-sh})
\le E_r e^{-sh}\sum_i
 \tau\{s(h-h_{i,s})\}.                              \tag{IC4A.5}
```

Although `h_(i,s)` depends on `s`, (IC4A.5) is an algebraic inequality at
each fixed `s`.  The subsequent proof differentiates only `K(s)`, not the
comparison.  No term `partial_s h_(i,s)` is missing.

## 3. Normalization and integration constant

Let `M(s)=E_r e^{-sh}` and retain the centered cumulant

```math
K(s)=\log M(s)+sE_rh.
```

Then

```math
{\operatorname{Ent}_r(e^{-sh})\over M(s)}
=-sE_{q_s}h-\log M(s)
=sK'(s)-K(s).                                       \tag{IC4A.6}
```

Dividing (IC4A.5) by `M(s)`, conditioning on `B_-i`, and using (IC4A.3)
gives

```math
sK'(s)-K(s)
\le\sum_iE_{q_s}D(q_s^{i\mid-i}\Vert r_i)
=s^2\mathcal E_s.                                   \tag{IC4A.7}
```

Consequently

```math
\left({K(s)\over s}\right)'\le\mathcal E_s.
```

Since `K(0)=K'(0)=0`, `K(s)/s ->0`.  Integrating from zero to `lambda`
therefore yields exactly

```math
\mathcal J=K(\lambda)
\le\lambda\int_0^\lambda\mathcal E_s\,ds.          \tag{IC4A.8}
```

The prefactor `lambda` is correct.

## 4. The zero-tilt limit

On the finite cube, conditional exponential tilting gives uniformly for
fixed `N`

```math
D(q_s^{i\mid-i}\Vert r_i)
={s^2\over2}\operatorname{Var}_{r_i}
 (h\mid B_{-i})+O(s^3).                             \tag{IC4A.9}
```

Also the outer marginal of `q_s` tends to `r_-i`.  Hence

```math
\lim_{s\downarrow0}\mathcal E_s
={1\over2}\sum_iE_{r_{-i}}
 \operatorname{Var}_{r_i}(h\mid B_{-i}),            \tag{IC4A.10}
```

which is IC.22.  No linear term survives.

## 5. Exact comparison with `A_s`

The claimed inequality can be strengthened.  For any row-independent
comparison `c_i(B_-i)`, direct use of (IC4A.2) gives

```math
E_{q_s^{i\mid-i}}\tau\{s(h-c_i)\}
=D(q_s^{i\mid-i}\Vert r_i)
 +\tau\{s(h_{i,s}-c_i)\}.                           \tag{IC4A.11}
```

Taking `c_i=bar h_i`, averaging over the `q_s` marginal of `B_-i`, summing,
and dividing by `s^2` proves the exact identity

```math
\boxed{
\mathcal A_s-\mathcal E_s
={1\over s^2}\sum_iE_{(q_s)_{-i}}
 \tau\{s(h_{i,s}-\bar h_i)\}\ge0.}                \tag{IC4A.12}
```

Thus `E_s<=A_s` is correct, and `h_(i,s)` is genuinely the optimal
conditional comparison for this entropy bound.

## 6. Converse scale and scope

If `J>=eta N`, (IC4A.8) implies

```math
\sup_{0\le s\le\lambda}\mathcal E_s
\ge {\eta N\over\lambda^2}.                        \tag{IC4A.13}
```

By continuity, a witnessing parameter may be chosen in `(0,lambda)` (the
integral cannot be large solely at the endpoint `s=0`).  For that parameter,

```math
\sum_iE_{q_s}D(q_s^{i\mid-i}\Vert r_i)
=s^2\mathcal E_s
\ge {\eta s^2N\over\lambda^2}.                     \tag{IC4A.14}
```

This verifies IC.24--IC.25.  It does not give an unscaled `Omega(N)` raw
conditional KL unless the witness `s` is bounded away from zero; the source
states this limitation correctly.

Finally, IC.4 remains a sufficient canonical-certificate criterion.  It
does not reverse the inequality, identify the best row-product shadow, or
prove the influence small for optimized children.  Its valid advance is to
replace an arbitrary row comparison by the exact conditional entropic one
and thereby isolate the weakest tensorization witness available along this
specific escort path.
