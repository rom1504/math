# Recurrence scope of the actual-child negative-overlap bound

Status: **rigorous implication and direction audit**.  This note propagates
the negative-window overlap estimate through every existing variational and
basin inequality.  It identifies a sign obstruction that is easy to miss:
small raw overlap controls the distance from the full escort to a row-product
shadow, but it does not show that either law reaches the composition target.
Conversely, a positive overlap floor gives no lower bound on product gain.

The conclusion is exact.  A power-saving overlap estimate would provide a
Hammersley-summable *conditional* error once target reach is independently
known.  It cannot itself create the missing almost-subadditive recurrence.

## 1. Three distinct gaps

Fix a comparable split `m+n=N`, raw temperature `t=beta/sqrt(N)`, actual
contracted-temperature minimizing children, and one orientation.  Let

```math
L(B)=\log\overline Z_N(A,\epsilon D,B;t),
\qquad
V_\lambda=-{1\over\lambda}\log E_Ue^{-\lambda L},       \tag{RC.1}
```

and let `V_lambda^row` be the same Gibbs variational problem restricted to
independent complete bridge rows.  For a declared composition target `T_N`,
define

```math
\begin{aligned}
G_N&=E_UL-V_\lambda^{\rm row},
   &&\text{fair-to-product gain},\\
I_N&=\lambda(V_\lambda^{\rm row}-V_\lambda),
   &&\text{reverse product dependence},\\
\Delta_N&=V_\lambda^{\rm row}-T_N,
   &&\text{product target excess}.
\end{aligned}                                             \tag{RC.2}
```

All three are different quantities.  The exact identities

```math
\boxed{
\Delta_N=(V_\lambda-T_N)+{I_N\over\lambda}
        =(E_UL-T_N)-G_N}                                  \tag{RC.3}
```

show both possible interfaces.  Controlling `G_N` says only that the
product target excess is close to the fair target excess.  Controlling
`I_N` transfers a separately known target estimate from the full escort to
the product shadow.

## 2. What an overlap power saving really gives

Theorem 37.50 proves for the raw actual-child path

```math
\begin{aligned}
0\le G_N
 &\le C_{\rm LS}\lambda t^2mn\widehat\rho_N^-(\lambda),\\
0\le I_N
 &\le C_{\rm LS}\lambda^2t^2mn
                         \widehat\rho_N^-(\lambda).       \tag{RC.4}
\end{aligned}
```

On comparable splits, `t^2mn=Theta_(beta)(N)`.  Hence

```math
\widehat\rho_N^-(\lambda)=O(N^{-\alpha})
\quad\Longrightarrow\quad
G_N=O(N^{1-\alpha}),\qquad I_N=O(N^{1-\alpha}).             \tag{RC.5}
```

If one **additionally** knows target reach at rate

```math
V_\lambda\le T_N+E_N,                                  \tag{RC.6}
```

then (RC.3)--(RC.5) give the exact product-target implication

```math
\boxed{
\Delta_N\le E_N+C_{\rm LS}\lambda\beta^2{mn\over N}
                 \widehat\rho_N^-(\lambda).}             \tag{RC.7}
```

Thus `E_N=O(N^(1-gamma))` and overlap exponent `alpha>0`
give a basin/product defect

```math
O\bigl(N^{1-\min\{\alpha,\gamma\}}\bigr).                \tag{RC.8}
```

Theorem 37.19 then supplies an `exp[-O(N)]` family of bridges at that
pressure accuracy.  A balanced-tree merge has summable error density for
every positive exponent.  More generally, define the dyadic-shell envelope

```math
\bar\epsilon(r)=\sup_{r\le M\le2r}\epsilon_M.
```

A rate `N epsilon_N` is usable when

```math
\lim_{k\to\infty}\sum_{j\ge0}\bar\epsilon(2^jk)=0.     \tag{RC.9}
```

If `epsilon_N` is nonincreasing, the shell envelope can equivalently be
replaced by its left endpoint.  Without such regularity, values only on one
exact dyadic ray do not control nearby orders in an all-order balanced tree.

This is the exact conditional recurrence consequence of an overlap bound.

## 3. Why this is not an independent recurrence

Condition (RC.6) already implies

```math
\boxed{
P_N(\beta)\le\min_B L(B)\le V_\lambda
             \le T_N+E_N.}                              \tag{RC.10}
```

The first inequality is parent minimization and the second is the elementary
soft-min bound.  Therefore target reach at a Hammersley-summable rate already
is the desired scalar pressure recurrence; overlap adds basin abundance or
row-product structure, not a sharper scalar recurrence.

Without (RC.6), neither line of (RC.4) upper-bounds `V_lambda-T_N` or
`Delta_N`: (RC.3) leaves the fair excess `E_UL-T_N` and the full target
excess completely uncontrolled.  Consequently raw overlap decay alone
cannot move Level 5 to Level 6.

The converse direction is also unavailable.  Equations (RC.4) are upper
bounds, so a lower bound on `widehat rho_N^-` gives no lower bound on `G_N`
or `I_N`.  The fixed-projective rank-one family of Theorem 37.53 has

```math
\widehat\rho_N^-(\lambda)\longrightarrow c_{\beta,\lambda}>0,
\qquad I_N=0,                                         \tag{RC.11}
```

because every negative tilt is exactly row product.  Hence the actual-law
positive floor in Theorems 37.52 and 37.56 is a genuine falsifier of the
proposed decay theorem, but not a favorable-dependence or recurrence
certificate.

## 4. Exact campaign verdict

The proposed SML

```math
\widehat\rho_N^-(\lambda)=O(N^{-\alpha})
```

would have the correct numerical scale *after* an independent target-reach
theorem, but its inequality direction cannot establish target reach.  It is
also false above the strong-channel threshold on every actual child law.
The only permissible continuation is therefore directional: subtract the
row-explainable support response and control a renormalized product error,
or build a rare-event/renormalization theorem which directly bounds the
target excess in the correct direction.
