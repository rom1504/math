# When can coarse row features certify extensive actual-child dependence?

Status: **rigorous follow-on theorem note**.  This note asks the scalable
question forced by the finite certificates in
[`actual_child_row_product_global_certificate.md`](actual_child_row_product_global_certificate.md):
when can a deterministic row-feature image prove
`I_lambda^leftarrow=Omega(N)` rather than merely a positive finite gap?

The answer has an exact information ceiling and two sufficient amplification
mechanisms.  It also identifies the additional theorem that an actual-child
argument would have to prove.

## 1. Coarse reverse-product information

Let `q` be a law on `product_(i=1)^m X_i`, let
`phi_i:X_i -> Y_i` be rowwise maps, and write

```math
 Q=(\phi_1,\ldots,\phi_m)_\#q,
 \qquad
 J(Q)=\inf_{P=\otimes_iP_i}D(P\Vert Q).                \tag{CA.1}
```

The data-processing theorem GC.2 says

```math
 \mathcal I^{\leftarrow}(q)
 :=\inf_{p=\otimes_i p_i}D(p\Vert q)\ge J(Q).          \tag{CA.2}
```

Thus `J(Q)` is precisely the amount of directed dependence exposed by the
declared row features.

## 2. An exact feature-budget ceiling

**Theorem CA.1 (coarse-atom ceiling).**  For every `Q`,

```math
 \boxed{
 0\le J(Q)
 \le-\log\max_yQ(y)
 \le\log|\operatorname{range}(\phi_1,\ldots,\phi_m)|
 \le\sum_{i=1}^m\log|Y_i|.}                           \tag{CA.3}
```

In particular, if the total coarse feature budget satisfies

```math
 \sum_i\log|Y_i|=o(N),                                \tag{CA.4}
```

then **no** GC.2 certificate using those features can prove
`I^leftarrow(q)=Omega(N)`.  If at most `r` rows have nontrivial features and
`|Y_i|<=K`, then `J(Q)<=r log K`.

*Proof.*  A point mass `delta_y` is a product law.  Choosing an atom of
maximum `Q` mass gives

```math
 J(Q)\le D(\delta_y\Vert Q)=-\log\max_yQ(y).
```

The largest atom has mass at least the reciprocal of the number of possible
coarse atoms.  The remaining bounds are immediate. `square`

This is the relevant distinction between a fixed **total** feature list and
a fixed number of features **per row**.  A fixed total list is incapable of
an extensive certificate.  One bit on each of `Theta(N)` rows has
`exp(Theta(N))` coarse atoms and can, in principle, certify a linear gap.

## 3. Exact and approximate block amplification

Partition the row indices into disjoint blocks `B in mathcal B`.  Let `Q_B`
be laws on the corresponding coarse coordinates and put
`Q_0=tensor_(B in mathcal B)Q_B`.

**Theorem CA.2 (one-sided block amplification).**

1. If `Q=Q_0`, then

   ```math
   \boxed{J(Q)=\sum_{B\in\mathcal B}J(Q_B).}           \tag{CA.5}
   ```

2. More generally, if

   ```math
   Q(y)\le e^{\varepsilon_N}Q_0(y)
   \quad\hbox{for every }y,                            \tag{CA.6}
   ```

   then

   ```math
   \boxed{J(Q)\ge\sum_{B\in\mathcal B}J(Q_B)
                    -\varepsilon_N.}                  \tag{CA.7}
   ```

Consequently, bounded-size blocks with
`sum_B J(Q_B)>=cN` and `epsilon_N=o(N)` expose
`I^leftarrow(q)>=cN-o(N)`.

*Proof.*  A row-product `P` factors across the blocks, and

```math
 D(P\Vert Q_0)=\sum_BD(P_B\Vert Q_B).
```

Minimization separates, proving (CA.5).  Under (CA.6),

```math
 D(P\Vert Q)
 =D(P\Vert Q_0)-\mathbb E_P\log{Q\over Q_0}
 \ge D(P\Vert Q_0)-\varepsilon_N.
```

Minimize over row products. `square`

A convenient sufficient form of (CA.6) is a coarse-potential decomposition

```math
 \log Q(y)=\sum_B\log Q_B(y_B)+R_N(y)-c_N
```

with `osc(R_N)=o(N)`: normalization then gives (CA.6) with an `o(N)`
exponent.  This is a genuine bounded-holonomy condition.  It does not require
the complete original bridge pressure table once the row maps and the
bounded-block potentials are known.

## 4. Contraction amplification without block independence

There is a second mechanism.  On the coarse alphabet let `f=-log Q`, let
`C` be its row-rectangle matrix, and consider

```math
 D(P\Vert Q)=\mathbb E_Pf-\sum_iH(P_i).                \tag{CA.8}
```

**Theorem CA.3 (contractive exposed-dependence certificate).**  Suppose a
product law `R=tensor_iR_i` satisfies the coordinate mean-field equations for
(CA.8).  If

```math
 \lambda_{\max}(C)<4,                                 \tag{CA.9}
```

then

```math
 \boxed{J(Q)=D(R\Vert Q).}                            \tag{CA.10}
```

Hence a checkable lower bound `D(R||Q)>=cN` proves
`I^leftarrow(q)>=cN`.

*Proof.*  This is GC.1 with `lambda=1`; the additive constants converting
`-H(P_i)` to `D(P_i||U_i)` do not change minimizers. `square`

CA.3 is useful only when `Q` has a compressed potential presentation from
which the fixed point, rectangle matrix, and cross-entropy can be computed.
Writing down an arbitrary table with `product_i|Y_i|` entries would merely
replace one exponential landscape by another.

## 5. The ceiling is sharp at one bit per row

Let `m` be even and, on each adjacent pair, put

```math
 Q_\rho(a,b)={1+\rho ab\over4},
 \qquad a,b\in\{-1,1\},
 \qquad 0<|\rho|<\tanh1.                              \tag{CA.11}
```

Set `Q=tensor_(j=1)^(m/2)Q_rho`.  The uniform product is the unique product
minimizer for each pair by CA.3: its sole rectangle oscillation is
`4|atanh rho|<4`.  Therefore

```math
 \boxed{
 J(Q)=-{m\over4}\log(1-\rho^2)=\Theta(m).}            \tag{CA.12}
```

Thus one binary feature on every row is enough for a linear certificate.
The obstruction in CA.1 is the **total** feature budget, not bounded alphabet
per row.

## 6. Implication for the actual optimized-child program

The finite `N=8,9` Walsh certificates prove that actual escorts expose
nonzero dependence in very small row alphabets.  They do not show
amplification.  The weakest concrete scalable replacement is now:

> **Actual-child coarse amplification lemma.**  Find fixed-alphabet row maps
> for comparable splits and either
>
> 1. bounded row blocks whose local reverse-product gaps sum to `cN`, while
>    the full coarse law obeys the one-sided domination (CA.6) with
>    `epsilon_N=o(N)`; or
> 2. a compressed coarse potential with a product fixed point `R_N`,
>    `lambda_max(C_N)<4-o(1)`, and `D(R_N||Q_N)>=cN`.

Either statement, together with GC.2, proves
`I_lambda^leftarrow>=cN-o(N)` for the actual optimized children.  Both retain
only constant alphabet per row plus a subextensive compatibility remainder.
They are therefore strictly smaller than the complete `2^(mn)` bridge
response table.

CA.1 also prevents a false next step: a bounded number of global parity or
overlap observables can never establish the required extensive resource,
regardless of how well chosen they are.  Any successful coarse certificate
must expose `Omega(N)` total feature bits (or an equivalent growing alphabet)
and prove that their local dependence survives composition rather than being
cancelled by an extensive compatibility term.
