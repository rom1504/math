# Joint global hybrid selection: exact variational form, scale wall, and an order-nine obstruction

Scratch memo for Wave 13, routes 6/7.  The exact finite claims are checked by
`tmp/global_hybrid_selection_r13.py`.  No tracked file is edited.

## 1. Exact pure selection problem

Let an order-`n` global minimizer be partitioned into vertex blocks
`V_1,...,V_k`, of orders `m_i`.  Write

```math
A=C+\bigoplus_iD_i,
```

where `C` has zero diagonal blocks and retains every original cross block.
Thus, for `z=(z_i)`,

```math
c(z):=z^{\mathsf T}Cz
=2\sum_{i<j}z_i^{\mathsf T}C_{ij}z_j.
```

Put

```math
\mathcal M_m=\{G:\ G\text{ is an order-}m\text{ signing and }Q(G)=q_m\}.
```

For `\mathbf G=(G_i)\in\prod_i\mathcal M_{m_i}`, let
`A_{\mathbf G}=C+\bigoplus_iG_i`.  The desired global hybrid excess is the
finite pure minimax value

```math
\boxed{
E(A,\mathcal P)
=\min_{\mathbf G\in\prod_i\mathcal M_{m_i}}
  \max_{\sigma\in\{\pm1\},\ z_i\in\{\pm1\}^{m_i}}
  \left\{
  \sigma c(z)+\sum_i\sigma z_i^{\mathsf T}G_iz_i-q_n
  \right\}.
}
\tag{R13.1}
```

Since `A` is globally minimizing, `E(A,\mathcal P)\ge0`.  This is a
`min max`, not a `max min`: the same tuple of block minimizers must satisfy
all global states.

Equivalently, with

```math
\delta_G(\sigma,z)=q_m-\sigma z^{\mathsf T}Gz,
\qquad
\Phi_{\mathbf G}
=\max_{\sigma,z}
\left\{\sigma c(z)-\sum_i\delta_{G_i}(\sigma,z_i)\right\},
```

one has

```math
\boxed{
E(A,\mathcal P)
=\min_{\mathbf G}\Phi_{\mathbf G}+\sum_iq_{m_i}-q_n.
}
\tag{R13.2}
```

For the original blocks, defining their deficits with `Q(D_i)` instead of
`q_{m_i}` gives

```math
\Phi_{\mathbf D}=q_n-\sum_iQ(D_i).
```

Consequently the exact common-mosaic identity, now minimized over all joint
choices, is

```math
\boxed{
\min_{\mathbf G}(\Phi_{\mathbf G}-\Phi_{\mathbf D})
=\underbrace{\sum_i(Q(D_i)-q_{m_i})}_{E_{\rm int}}
 +\underbrace{E(A,\mathcal P)}_{E_{\rm hyb}}.
}
\tag{R13.3}
```

This separates the compulsory internal excess from the genuinely new global
hybrid excess.

There is an exact mixed LP relaxation.  Let `\Omega` be the finite set of
oriented projective global states and set

```math
f(\mathbf G,\omega)
=\sigma c(z)+\sum_i\sigma z_i^{\mathsf T}G_iz_i-q_n.
```

Finite minimax gives

```math
\boxed{
\begin{aligned}
E_{\rm mix}
&:=\min_{\mu\in\Delta(\prod_i\mathcal M_{m_i})}
       \max_{\omega}\mathbb E_\mu f(\mathbf G,\omega)\\
&=\max_{p\in\Delta(\Omega)}
\left[
\mathbb E_p(\sigma c(z)-q_n)
+\sum_i\min_{G_i\in\mathcal M_{m_i}}
          \mathbb E_p\,\sigma z_i^{\mathsf T}G_iz_i
\right]
\le E(A,\mathcal P).
\end{aligned}
}
\tag{R13.4}
```

The dual distribution identifies relevant shared states, but it does not by
itself round to one deterministic tuple.  Signed-permutation orbit
restrictions are obtained by replacing each `\mathcal M_{m_i}` by the desired
orbit; distinct switching/permutation classes must not silently be omitted.

## 2. A deterministic theorem and the first nontrivial scale

Let `s=\max_i m_i` and `H=A_{\mathbf G}-A=\bigoplus_i(G_i-D_i)`.  Block
subadditivity is immediate but worth recording with the correct factor:

```math
Q(H)
=\max_z\left|\sum_i z_i^{\mathsf T}(G_i-D_i)z_i\right|
\le\sum_iQ(G_i-D_i).
```

The triangle inequality and the entrywise bound for an order-`m_i` signing
give

```math
Q(G_i-D_i)
\le q_{m_i}+Q(D_i)
\le2m_i(m_i-1).
```

As `Q(A)=q_n` and `Q(A+H)\le Q(A)+Q(H)`, **every** joint choice of
minimizers satisfies

```math
\boxed{
0\le Q(A_{\mathbf G})-q_n
\le\sum_i\{q_{m_i}+Q(D_i)\}
\le2\sum_i m_i(m_i-1)
\le2n(s-1).
}
\tag{R13.5}
```

There is a sharper instance-dependent Hamming version.  If `d_i` is the
number of undirected edges on which `D_i` and `G_i` differ, then one changed
edge contributes `\pm4` to `z^{\mathsf T}(G_i-D_i)z`, in the present doubled
normalization.  Hence

```math
\boxed{
Q(A_{\mathbf G})-q_n
\le4\sum_i d_i,
\qquad
E(A,\mathcal P)
\le4\sum_i\operatorname{dist}_{\rm Ham}
       (D_i,\mathcal M_{m_i}).
}
\tag{R13.5a}
```

Thus a total orbit/minimizer Hamming distance `o(n^{3/2})` is another exact
sufficient condition.  No theorem currently relates the scalar excess
`Q(D_i)-q_{m_i}` to this much stronger distance, so (R13.5a) is a target,
not a closure of the critical-scale problem.

Thus `E(A,\mathcal P)=o(n^{3/2})` is automatic whenever
`\sum_i m_i^2=o(n^{3/2})`, in particular for `s=o(\sqrt n)`.  No
randomization is needed in that range.

This does not yet harvest leading-order information.  Indeed

```math
0\le E_{\rm int}
\le\sum_i m_i(m_i-1)
\le n(s-1).
\tag{R13.6}
```

So below `\sqrt n`, both the error and the entire internal excess that the
replacement identity can expose are already `o(n^{3/2})`.  For roughly
equal blocks of order `n^\alpha`, the deterministic error is
`O(n^{1+\alpha})`; the first nontrivial scale is
`\alpha=1/2`.  A useful critical-scale theorem must improve the `O(ns)`
bound, not merely change its constant.

## 3. What signed-switching randomization can rigorously say

Fix one minimizer `H_i\in\mathcal M_{m_i}` in each block and independently
apply a uniform vertex switching `S_i`; permutations may also be included but
are unnecessary for the following estimate.  Put

```math
G_i=S_iH_iS_i,
\qquad
B=\sum_iq_{m_i},
\qquad
V=\sum_iq_{m_i}^2.
```

For every fixed oriented global state `\omega=(\sigma,z)`, the variables

```math
X_{i,\omega}=\sigma z_i^{\mathsf T}G_iz_i
```

are independent, centered, and lie in `[-q_{m_i},q_{m_i}]`.  Centering is
exact because uniform vertex switching makes `S_iz_i` a uniform Rademacher
spin and every diagonal is zero.  Hoeffding therefore gives

```math
\Pr\left\{\sum_iX_{i,\omega}>r\right\}
\le\exp\left(-\frac{r^2}{2V}\right)
\qquad(r>0).
\tag{R13.7}
```

For a proposed excess `t`, define the cross margin

```math
r_\omega=q_n+t-\sigma c(z)
```

and the only states that are not deterministically harmless,

```math
\Omega_t=\{\omega:0<r_\omega<B\}.
```

If every margin is positive and

```math
\boxed{
\sum_{\omega\in\Omega_t}
\exp\left(-\frac{r_\omega^2}{2V}\right)<1,
}
\tag{R13.8}
```

then some deterministic tuple of signed-switching orbit representatives has
`Q(A_{\mathbf G})\le q_n+t`.  This is the honest low-slack-state union bound:
states with `r_\omega\ge B` cannot violate the target, and no enumeration of
the other states is needed.

Let `\Delta_C=[Q(C)-q_n]_+`.  Taking `t=\Delta_C+u` makes every margin at
least `u`.  If `N_t=|\Omega_t|`, the cruder sufficient condition is

```math
u>\sqrt{2V\log N_t}.
\tag{R13.9}
```

For equal blocks of order `s`, the standard random-signing upper bound
`q_s=O(s^{3/2})` gives `V=O(ns^2)`, hence fluctuation cost

```math
O\!\left(s\sqrt{n\log N_t}\right).
\tag{R13.10}
```

At the critical scale `s\asymp\sqrt n`, this is `o(n^{3/2})` precisely if
the relevant profile entropy satisfies `\log N_t=o(n)`; one also needs
`\Delta_C=o(n^{3/2})`.  More generally, at `s=n^\alpha` the condition is
`\log N_t=o(n^{2-2\alpha})`.  With the full `2^{\Theta(n)}` state set at
the critical scale, the bound is only `O(n^{3/2})`.  Thus two genuinely new
inputs are required: a small cross-only overshoot and a subcritical entropy
bound for the near-top cross profiles.  Uniform orbit averaging alone proves
neither.

## 4. Exact order-nine global selector wall

The order-nine minimizer in ledger (10.298), in the displayed vertex order,
is partitioned as

```math
I=\{1,2,3\},\qquad J=\{4,5,6,7,8,9\}.
```

Exact evaluation gives

```math
Q(A)=q_9=24,
\qquad
(Q(D_I),Q(D_J))=(6,14),
\qquad
(q_3,q_6)=(6,10),
\qquad
Q(C)=24.
\tag{R13.11}
```

Despite zero cross-only overshoot, joint optimal-block replacement cannot
preserve the parent optimum:

```math
\boxed{
\min_{G_3\in\mathcal M_3,\ G_6\in\mathcal M_6}
Q\!\left(C\oplus_{\rm diag}(G_3,G_6)\right)
=28.
}
\tag{R13.12}
```

Here is a compact six-state certificate.  Set

```math
x_a=(1,-1,1),\quad x_b=(1,1,1),
```

```math
y_a=(1,1,1,1,1,-1),\quad
y_b=(1,1,1,1,-1,-1).
```

The cross quadratic has values `\pm24` on `(x_a,\pm y_a)`, `\pm20` on
`(x_a,\pm y_b)`, and `\pm20` on `(x_b,\pm y_b)`.  Put

```math
a=x_a^{\mathsf T}G_3x_a,
\quad b=x_b^{\mathsf T}G_3x_b,
\quad u=y_a^{\mathsf T}G_6y_a,
\quad v=y_b^{\mathsf T}G_6y_b.
```

If the full hybrid norm were at most 24, the two signs of those three
projective pairs would force simultaneously

```math
a+u=0,
\qquad |a+v|\le4,
\qquad |b+v|\le4.
\tag{R13.13}
```

All eight labelled order-three minimizers have the following two-state
profiles (multiplicities in the second column):

| `(a,b)` | count |
|---|---:|
| `(-6,2)` | 1 |
| `(-2,-2)` | 2 |
| `(-2,6)` | 1 |
| `(2,-6)` | 1 |
| `(2,2)` | 2 |
| `(6,-2)` | 1 |

All 384 labelled order-six minimizers have:

| `(u,v)` | count |
|---|---:|
| `(-10,-6)` | 60 |
| `(-10,10)` | 12 |
| `(-6,-10)` | 60 |
| `(-6,6)` | 60 |
| `(6,-6)` | 60 |
| `(6,10)` | 60 |
| `(10,-10)` | 12 |
| `(10,6)` | 60 |

The tables immediately contradict (R13.13).  When `a=\pm2`, the equality
`a+u=0` is impossible.  For `(a,b)=(-6,2)`, equality leaves
`(u,v)=(6,-6)` or `(6,10)`; respectively the second or third inequality
fails by giving 12.  The `a=6` case is symmetric.  Thus every one of the
`8\cdot384=3072` joint choices has a violation among these six global
states.  Since every order-nine signing energy is a multiple of four, the
norm is at least 28.

Equality is attained, for example, by

```math
G_3=
\begin{pmatrix}
0&-1&-1\\
-1&0&-1\\
-1&-1&0
\end{pmatrix},
\qquad
G_6=
\begin{pmatrix}
0&-1&-1&1&-1&1\\
-1&0&-1&1&1&-1\\
-1&-1&0&-1&-1&-1\\
1&1&-1&0&-1&-1\\
-1&1&-1&-1&0&1\\
1&-1&-1&-1&1&0
\end{pmatrix}.
\tag{R13.14}
```

Direct enumeration gives hybrid energy histogram

```text
-28:5, -20:20, -12:47, -4:49, 4:63, 12:50, 20:21, 28:1.
```

Exactly 672 of the 3,072 pairs attain norm 28.

The checker obtains every labelled minimizer, not merely one orbit sample:
it switches an arbitrary signing to first-row-positive gauge, exhausts the
`2^{\binom{m-1}{2}}` gauge matrices, retains every gauge minimizer, and then
expands all `2^{m-1}` switchings.  Hence the profile tables cover all of
`\mathcal M_3` and `\mathcal M_6`; permutations are already included among
the labelled matrices.

For this example the exact decomposition (R13.3) is

```math
E_{\rm int}=(6-6)+(14-10)=4,
\qquad E_{\rm hyb}=28-24=4,
```

```math
\Phi_{\mathbf D}=24-6-14=4,
\qquad
\min\Phi_{\mathbf G}=28-6-10=12,
```

so the common-mosaic response rises by eight, exactly `4+4`.

This is not a one-state or local-layer obstruction.  If the minimizing
blocks are allowed to depend on the state, exact enumeration gives

```math
\max_z\min_{G_3,G_6}
\left|c(z)+z_I^{\mathsf T}G_3z_I+z_J^{\mathsf T}G_6z_J\right|=8,
```

far below 24; simultaneous satisfaction of six shared global states is the
issue.  Also, for the tail-local functional of (10.520),
`J_C(D_6)=30` and `\min_{G_6}J_C(G_6)=30` (the corresponding layers are
16 and 20), while the reverse head-local minimum is 26 rather than the
original 30.  These separately optimized local numbers neither add to nor
determine the global pure minimax value.

## 5. Small-order audit and exact scope

The same checker found zero best hybrid excess in:

- every bipartition of every minimizing gauge representative through order
  six (case counts by order 2--6: `1,6,42,180,372`);
- all 3,240 order-seven minimizing gauges for one canonical split of each
  size `1+6`, `2+5`, `3+4` (switching/permutation invariance makes this cover
  every labelled bipartition);
- the displayed order-eight minimizer for splits `2+6`, `3+5`, `4+4`;
- the displayed order-nine minimizer for its `4+5` split.

The `3+6` order-nine wall is therefore the first obstruction found in this
audit.  It proves that joint selection need not preserve exact global
optimality, even with all minimizing classes available and `Q(C)=q_n`.

It does **not** falsify the desired asymptotic statement
`E(A,\mathcal P)=o(n^{3/2})`: its excess is the fixed value four, and no
globally minimizing blow-up of the witness has been proved.  The strongest
rigorous asymptotic conclusion here remains (R13.5), plus the conditional
entropy criterion (R13.8).  At the critical block scale, the next concrete
target is to prove both `Q(C)-q_n=o(n^{3/2})` and subexponential entropy of
the cross profiles within the needed margin, or to construct a genuine
critical-scale family violating one of those requirements.
