# Internal-block replacement in a global minimizer

Scratch memo; no tracked-file claim is made here.

Write

```math
q_m=\min\{Q(B):B\text{ is an order-}m\text{ signing}\},
\qquad q_0=q_1=0.
```

## 1. One-block replacement

Let `A` be an order-`n` minimizer, so `Q(A)=q_n`, and split

```math
A=\begin{pmatrix}B&C\\ C^{\mathsf T}&D\end{pmatrix},
\qquad |B|=s.
```

Choose any order-`s` minimizer `B_*`.  For oriented states put

```math
\delta_{B_*}(\sigma,x)=q_s-\sigma x^{\mathsf T}B_*x,
\qquad
\delta_D(\sigma,y)=Q(D)-\sigma y^{\mathsf T}Dy.
```

Both deficits are nonnegative.  The exact coupled response is

```math
\Xi_{B_*\to D}(C)
=\max_{\sigma,x,y}
\left[
2\sigma x^{\mathsf T}Cy
-\delta_{B_*}(\sigma,x)
-\delta_D(\sigma,y)
\right].
```

For the replacement signing
`A_*=[[B_*,C],[C^T,D]]` one has the identity

```math
Q(A_*)=q_s+Q(D)+\Xi_{B_*\to D}(C).
```

Since `A` is globally minimizing, `Q(A_*)>=q_n`; hence

```math
\boxed{
\Xi_{B_*\to D}(C)\ge q_n-q_s-Q(D).
}
```

This holds for every choice, switching, and global negation of the
lower-order minimizer.  It is the strongest direct statement of the
replacement comparison: all terms discarded from it are nonnegative
endpoint deficits.

Define the ordinary all-successor layer capacity

```math
\mathcal L_{B\to D}
=\max_{\sigma,y}
\left[
2\|Cy\|_1-\delta_D(\sigma,y)
\right]_+.
```

Maximizing the cross term over `x` and dropping `delta_{B_*}` gives
`Xi<=mathcal L`, and therefore

```math
\boxed{
\mathcal L_{B\to D}
\ge [q_n-q_s-Q(D)]_+.
}
```

If `d_B=q_n-Q(D)` is the principal-block decrement, this is simply

```math
\mathcal L_{B\to D}\ge[d_B-q_s]_+.
```

Let

```math
e_B=Q(B)-q_s,
\qquad
\alpha=Q(B)+Q(D)-Q(A).
```

Thus `alpha` is the absolute-orientation additivity defect of this split.
The same inequality is

```math
\boxed{
\mathcal L_{B\to D}\ge[e_B-\alpha]_+.
}
```

At an endpoint split, separate `P/N` superadditivity gives

```math
\alpha
\le\frac{I(B)+I(D)-I(A)}2.
```

Consequently the endpoint-only weakening is

```math
\mathcal L_{B\to D}
\ge
\left[
e_B-\frac{I(B)+I(D)-I(A)}2
\right]_+.
```

At an active orientation-compatible zero cut (`alpha=0`), the whole
internal optimality excess `Q(B)-q_s` must therefore reappear in a genuine
successor layer.  This is non-vacuous and is not the cross-only comparison
of Section 10.55.1.

There is also an orientation-refined necessary inequality.  If
`(P_*,N_*)` are the two oriented extrema of `B_*` and

```math
\mathcal L_+^D=\max_y[2\|Cy\|_1-(P(D)-y^TDy)]_+,
\quad
\mathcal L_-^D=\max_y[2\|Cy\|_1-(N(D)+y^TDy)]_+,
```

then

```math
q_n\le
\max\{P_*+P(D)+\mathcal L_+^D,
       N_*+N(D)+\mathcal L_-^D\}.
```

The same inequality with `P_*` and `N_*` interchanged follows by replacing
with `-B_*`.

## 2. Replacing both children

Replace both `B` and `D` by arbitrary minimizers `B_*` and `D_*`.  Define

```math
\Xi_{B_*,D_*}(C)
=\max_{\sigma,x,y}
\left[
2\sigma x^TCy
-\delta_{B_*}(\sigma,x)
-\delta_{D_*}(\sigma,y)
\right].
```

The same exact identity and global comparison give

```math
\boxed{
\Xi_{B_*,D_*}(C)
\ge q_n-q_s-q_{n-s}
=e_B+e_D-\alpha.
}
```

This holds for every pair of lower-order minimizers.  In particular,

```math
2\|C\|_{\infty\to1}\ge q_n-q_s-q_{n-s}.
```

At an active cut `alpha=0`, the cross block must recover the sum of both
internal optimality excesses while also paying the two minimizer endpoint
deficits.  This is the sharp coupled version of the one-block statement.

## 3. Simultaneous laminar replacements

Let

```math
U_t=H_t\sqcup U_{t+1},
\qquad
A_t=A[U_t]
=\begin{pmatrix}D_t&C_t\\ C_t^T&A_{t+1}\end{pmatrix}
```

be a nested tower, with `A_0=A`.  Put `h_t=|H_t|`, `m_L=|U_L|`, and

```math
d_t=Q(A_t)-Q(A_{t+1}).
```

Choose arbitrary order-`h_t` minimizers `G_t` and an order-`m_L`
minimizer `G_L`.  Define the hybrid replacement tower recursively by

```math
\widetilde A_L=G_L,
\qquad
\widetilde A_t
=\begin{pmatrix}G_t&C_t\\ C_t^T&\widetilde A_{t+1}\end{pmatrix},
\qquad
\widetilde Q_t=Q(\widetilde A_t).
```

At layer `t` define

```math
\widetilde\Xi_t
=\max_{\sigma,x,y}
\left[
2\sigma x^TC_ty
-\bigl(q_{h_t}-\sigma x^TG_tx\bigr)
-\bigl(\widetilde Q_{t+1}
       -\sigma y^T\widetilde A_{t+1}y\bigr)
\right]
```

and

```math
\widetilde{\mathcal L}_t
=\max_{\sigma,y}
\left[
2\|C_ty\|_1
-\bigl(\widetilde Q_{t+1}
       -\sigma y^T\widetilde A_{t+1}y\bigr)
\right]_+.
```

Exactly,

```math
\widetilde\Xi_t
=\widetilde Q_t-q_{h_t}-\widetilde Q_{t+1},
\qquad
\widetilde\Xi_t\le\widetilde{\mathcal L}_t.
```

The first identity telescopes.  Since `widetilde A_0` is another order-`n`
signing, global minimality gives `widetilde Q_0>=q_n`, and hence

```math
\boxed{
\sum_{t<L}\widetilde{\mathcal L}_t
\ge q_n-q_{m_L}-\sum_{t<L}q_{h_t}.
}
```

Equivalently, because the original decrements telescope,

```math
q_n-q_{m_L}-\sum_tq_{h_t}
=\sum_t(d_t-q_{h_t})+[Q(A_L)-q_{m_L}].
```

If the terminal block is left unchanged, the terminal excess is omitted:

```math
\sum_t\widetilde{\mathcal L}_t
\ge q_n-Q(A_L)-\sum_tq_{h_t}
=\sum_t(d_t-q_{h_t}).
```

There is a single-witness version.  The off-diagonal edges of the block
partition occur exactly once in the triangular mosaic

```math
2\sum_{t<L}\sigma x_t^TC_ty_{t+1}.
```

Expanding `Q(widetilde A_0)` gives

```math
\max_{\sigma,(x_t),y_L}
\left[
2\sum_t\sigma x_t^TC_ty_{t+1}
-\sum_t\delta_{G_t}(\sigma,x_t)
-\delta_{G_L}(\sigma,y_L)
\right]
\ge q_n-q_{m_L}-\sum_tq_{h_t}.
```

Thus one common nested successor state charges the full scalar partition
defect.  Dropping the nonnegative deficits yields the weaker raw exposure
consequence

```math
2\sum_t\|C_ty_{t+1}\|_1
\ge q_n-q_{m_L}-\sum_tq_{h_t}
```

for some common state.

This common witness does **not** come from choosing the maximizer in each
`widetilde Xi_t` independently.  It follows separately by expanding the
single global variational formula for `Q(widetilde A_0)`:

```math
Q(\widetilde A_0)-\sum_{t<L}q_{h_t}-q_{m_L}
=\max_{\sigma,(x_t),y_L}
\left[
2\sum_t\sigma x_t^TC_ty_{t+1}
-\sum_t\delta_{G_t}(\sigma,x_t)
-\delta_{G_L}(\sigma,y_L)
\right].
```

Thus the orientation and all nested tail spins really are common.

## 4. Exact remaining obstruction

The telescope is for the hybrid successors `widetilde A_{t+1}`.  It does
not imply

```math
\sum_t\mathcal L_t(A_{t+1})
\ge q_n-q_{m_L}-\sum_tq_{h_t}
```

for the original successor deficits.  Replacing the later diagonal blocks
can completely change the near-ground layers and endpoint tree of the tail.
Likewise, at a proper descendant `U_t`, global minimality compares a
replacement in `H_t` against its full complement in `V`, not merely against
the sibling `U_{t+1}`; cross edges to earlier blocks cannot be discarded.

Therefore the precise next lemma needed to combine replacement with the
verified `K<=4` path cover is a hybrid-to-original layer-stability estimate,
or a simultaneous choice of the lower-order minimizers for which the hybrid
layer rewards are dominated by the original endpoint-chain layers plus a
summable error.  Global minimality alone gives the exact hybrid telescope
above and no direct comparison of the two deficit landscapes.

A precise sufficient condition is

```math
\sup_{\sigma,y}
\left[
\delta_{A_{t+1}}(\sigma,y)
-\delta_{\widetilde A_{t+1}}(\sigma,y)
\right]_+
\le\varepsilon_t.
```

It implies immediately

```math
\widetilde{\mathcal L}_t
\le\mathcal L_t+\varepsilon_t.
```

Hence a choice of replacement minimizers with
`sum_t epsilon_t=o(n^(3/2))` would transfer the exact hybrid telescope to
the original layer capacities.  A stronger but simpler sufficient condition
is uniform control of the changes in both `Q(A_{t+1})` and every oriented
Boolean energy of the successor.

## 5. Exact failure of monotone hybrid-to-original stability

Even inside a global minimizer and at an endpoint split, replacing a child
by another minimizer can increase the layer capacity sharply.  Consider the
order-six minimizer

```math
A=\begin{pmatrix}
0&1&1&1&1&1\\
1&0&-1&-1&1&1\\
1&-1&0&1&-1&1\\
1&-1&1&0&1&-1\\
1&1&-1&1&0&-1\\
1&1&1&-1&-1&0
\end{pmatrix},
\qquad Q(A)=q_6=10.
```

One endpoint pair has the split `T={0,1,2,3,4}`, `S={5}`.  Its cross row is

```math
C=(1,1,1,-1,-1),
```

and the original five-block is

```math
D=\begin{pmatrix}
0&1&1&1&1\\
1&0&-1&-1&1\\
1&-1&0&1&-1\\
1&-1&1&0&1\\
1&1&-1&1&0
\end{pmatrix},
\qquad Q(D)=q_5=8.
```

Direct enumeration gives

```math
\mathcal L(C,D)=2.
```

Replace `D` by the equally optimal order-five signing

```math
D'=\begin{pmatrix}
0&-1&-1&1&1\\
-1&0&1&1&-1\\
-1&1&0&-1&1\\
1&1&-1&0&-1\\
1&-1&1&-1&0
\end{pmatrix},
\qquad Q(D')=q_5=8.
```

Then

```math
\boxed{\mathcal L(C,D')=10>2=\mathcal L(C,D).}
```

For `D'`, orientation `sigma=-1` and
`y=(-1,-1,-1,1,1)` are an exact negative endpoint and give cross exposure
ten with zero deficit.  For `D`, exhaustive evaluation of its sixteen
projective states in both orientations gives maximum payoff two.  The full
replacement signing has `Q=18`, compared with `Q(A)=10`.

Thus neither monotonicity of layer capacities nor a harmless bound on
`Q(widetilde A_0)-q_n` follows merely because every replacement block is
itself minimizing.  Any useful stability theorem must choose the minimizers
compatibly with the fixed cross mosaic; arbitrary lower-order minimizers are
already ruled out at order six.
