# Wave 48A: every row-good fixed-density complement column has no-saving mass

## Result

The entire fixed-density row-good complement-column implementation proposed
in `STEERING.md` is incompatible with the already verified fixed-slice upper
tail (10.1022), once its variables are matched.  Complement activity alone,
with no positive slack assumption and no scalar optimality, is a
`Theta(q_n)` deviation of the centered principal-selector payoff.  Hence any
complement column with the project row bound has at most
`exp{-Omega(n^(3/4))}` mass.

A secondary, weaker consequence combines scalar optimality and arithmetic
high-slack retention to show that column mass times conditional high-slack
mass pays the same unsaved exponent.  The direct active-event obstruction is
stronger and should be the strategic conclusion.

This is a rigorous obstruction to the structured complement implementation,
not to the bare arbitrary-cut target.  It uses a common fixed cut throughout
and therefore says nothing about cuts outside the complement-incidence
columns.

## Setup

Let `A` be an exact order-`n` minimizer and `q=Q(A)=q_n`.  Let

```math
p=m/n\in[p_0,p_1]\subset(1/\sqrt2,1),\qquad
p_2=(m)_2/(n)_2,qquad a_{n,m}=2p_2-1.
```

For an oriented cut `d=(sigma,x)`, write

```math
E=\langle A,d\rangle,qquad
c_A(S,d)=\langle A[S],d[S]\rangle,qquad
X_d(S)=c_A(S,d)-p_2E.
```

Retain

```math
B_S=-A+2P_SAP_S,\qquad
L_S(d)=\langle B_S,d\rangle=2c_A(S,d)-E,
```

and the active column and its slack

```math
\mathcal I_d=\{S:L_S(d)\ge q\},\qquad
\alpha_d=U_m(\mathcal I_d),\qquad
\Gamma_S=L_S(d)-q\quad(S\in\mathcal I_d).
```

Choose the fixed noncentral switching layer from (10.1185)--(10.1190), with
`k/n -> kappa` and

```math
\theta=\frac{(n-2k)^2-n}{n(n-1)}>1/3
```

for all sufficiently large `n`.  Let `d` minimize

```math
F_\lambda(d)=-\log\alpha_d+\lambda R_2(d),\qquad
\lambda=b_\kappa n^{-3/2},
```

over nonempty columns, with the fixed sufficiently small `b_kappa` used in
(10.1190).  Define the conditional high-slack fraction

```math
\beta_d=\Pr_{S\mid\mathcal I_d}
\{\theta\Gamma_S>(1-\theta)q\}.
```

The conclusion below is uniform when `p` stays in the displayed compact
window and `k/n` stays in a compact subwindow of the noncentral range used by
(10.1190).

## Exact omitted-block identity

Let `T=S^c`, and let `d^T=(sigma,D_Tx)` be the physical-spin switch on `T`.
The factorization (10.1206), with the ordered matrix pairing, gives

```math
L_S(d)=\langle A,d^T\rangle-2c_A(T,d).
```

Indeed, `A^{E(T)}` differs from `A` by `-2A[T]`, while flipping every spin
in `T` leaves the internal cut matrix `d[T]` unchanged.  Thus, with

```math
\Delta_T=q-\langle A,d^T\rangle\ge0,
```

one has the exact identity

```math
\boxed{\Gamma_S=-\Delta_T-2c_A(T,d).}
```

Consequently

```math
-c_A(T,d)=\frac{\Gamma_S+\Delta_T}{2}\ge\frac{\Gamma_S}{2},
\qquad
Q(A[T])\ge\frac{\Gamma_S}{2}.
```

In particular, every arithmetic high-slack incidence has

```math
Q(A[T])>\frac{1-\theta}{2\theta}q.
```

This is the exact promised mapping to omitted-block principal norm.  The
important additional fact is that the same incidence is also a fixed-cut
large deviation on the retained block, to which (10.1022) applies.

## High slack is a constant-scale principal-selector deviation

The identities above give

```math
X_d(S)=\frac{L_S(d)-a_{n,m}E}{2}.
```

In fact, no high slack is needed.  On the full active event
`S in I_d`, one has `L_S(d)>=q`.  Since `p_0>1/sqrt(2)`, one has
`a_(n,m)>=0` for all sufficiently large `n`; also `E<=q` and
`p_2<=p^2<=p_1^2`.  It follows pointwise that

```math
\boxed{
S\in\mathcal I_d
\quad\Longrightarrow\quad
X_d(S)
\ge\frac{(1-a_{n,m})q}{2}
=(1-p_2)q
\ge(1-p_1^2)q.}
```

This mapping uses only fixed-density complement activity.  In particular,
it applies to every oriented cut, whether or not it is a scalar optimizer.

Apply the verified fixed-slice tail (10.1022) with
`u=(1-p_1^2)q`.  Its small centering correction is

```math
b_0=(p^2-p_2)q=\frac{p(1-p)}{n-1}q=o(q).
```

The exact-minimizer bounds `q>=c_0 n^(3/2)` and
`||A||_op<=sqrt(2q)` therefore give, uniformly in the compact density
window,

```math
\boxed{
\alpha_d
\le 3(n+1)\exp\left\{-c\min\left(
\frac{q^2}{R_2(d)},\ n,\ n^{3/4}
\right)\right\}.}
```

Consequently, for every fixed `c_1>0`, the project row estimate

```math
R_2(d)=O(n^{9/4-c_1})
```

forces

```math
\boxed{-\log\alpha_d=\Omega(n^{3/4}).}
```

Indeed `q^2/R_2(d)=Omega(n^(3/4+c_1))`, while the other two exponents are
`Omega(n)` and `Omega(n^(3/4))`.  This contradicts the project-saving mass
requirement

```math
-\log\alpha_d=O(n^{3/4-c_1}).
```

Thus **no fixed-density complement-incidence cut can satisfy the mass/row
pair needed for (10.795)**.  This conclusion requires neither positive
slack, conditional retention, scalar optimality, nor the arithmetic switching
theorem.

For comparison, on the smaller high-slack event one gets the strict
threshold below.

On the high-slack event, `L_S(d)>q/theta`.  Since `p_0>1/sqrt(2)`, one has
`a_(n,m)>=0` for all sufficiently large `n`; also `E<=q` and
`p_2<=p^2<=p_1^2`.  Hence

```math
X_d(S)>
\frac{(1/\theta-a_{n,m})q}{2}
\ge (1-p_1^2)q=:\delta_0q.
```

The final constant is positive because `theta<=1` and `p_1<1`.  Therefore
the unconditional high-slack set `J_d` obeys

```math
U_m(J_d)=\alpha_d\beta_d,
\qquad
J_d\subseteq\{S:X_d(S)>\delta_0q\}.
```

The same fixed-slice tail then yields constants depending only on the fixed
windows such that

```math
\boxed{
\alpha_d\beta_d
\le Cn\exp\left\{-c\min\left(
\frac{q^2}{R_2(d)},\ n,\ n^{3/4}
\right)\right\}.}
```

Here `q^2/n^2>=c n`, and
`q/||A||_op>=sqrt(q/2)>=c n^(3/4)`.

## The arithmetic theorem bounds the same cut's row

For completeness, the strict inequality in the high-slack definition does
not hide a vanishing non-polynomial Cantelli cost.  Put

```math
\delta_S(k)=\theta(q+\Gamma_S)-q.
```

All energies are integers and `theta` has denominator `n(n-1)`.  Thus
`delta_S(k)>0` implies

```math
\delta_S(k)\ge 1/[n(n-1)].
```

Moreover every switched energy has absolute value at most
`Q(B_S)<=q+2Q(A[T])<=3q<=3n(n-1)`, so its variance is at most `9n^4`.
The Cantelli factor in (10.1186) is consequently at least `c n^(-8)` on
every high-slack incidence.  Averaging over the active column gives

```math
\overline p_k(d)\ge c n^{-8}\beta_d.
```

Writing `g_d=-log beta_d` (and `g_d=+infinity` if `beta_d=0`), the verified
scalar-optimality row theorem (10.1190) therefore gives

```math
\boxed{
R_2(d)\le C_\kappa
\{n^2+n^{3/2}(g_d+O(\log n))\}.}
```

## No-saving product theorem

Let `h_d=-log alpha_d`.  If `g_d>=n^(3/4)`, the desired conclusion is
immediate.  Otherwise the preceding row bound gives

```math
R_2(d)=O(n^{9/4}).
```

Before splitting into cases, the two boxed estimates give the more
quantitative tradeoff

```math
h_d+g_d+O(\log n)
\ge c\min\left\{
\frac{n^3}{n^2+n^{3/2}(g_d+O(\log n))},\ n,\ n^{3/4}
\right\}.
```

This form makes explicit that improving the conditional retention exponent
is exactly what forces the same scalar optimizer into a row regime where the
fixed-slice deviation has the unsaved exponent.

Insert this in the fixed-slice upper tail.  All three exponents there are
`Omega(n^(3/4))`, so

```math
e^{-(h_d+g_d)}=\alpha_d\beta_d
\le Cn e^{-c n^{3/4}}.
```

After absorbing `log n`, both cases prove the **Verified consequence of
(10.1022) and (10.1190)**

```math
\boxed{
-\log(\alpha_d\beta_d)=h_d+g_d\ge c_* n^{3/4}
}
```

for all sufficiently large `n`.  The constant `c_*>0` may depend on the
compact density windows and on `kappa`, but not on `n`, `A`, or the scalar
optimizer.

Equivalently, there is no `c>0` for which the same scalar optimizer can obey

```math
-\log\alpha_d=O(n^{3/4-c}),\qquad
\beta_d\ge\exp\{-O(n^{3/4-c})\}.
```

Indeed those two statements give
`h_d+g_d=O(n^(3/4-c))=o(n^(3/4))`, contradicting the boxed product wall.

## Scope and research judgment

1. The direct active-event theorem falsifies the entire fixed-density
   row-good complement-column implementation in the current `STEERING.md`,
   not merely its arithmetic mass/high-slack package.  Moving away from the
   near diagonal does not avoid the no-saving exponent.
2. The omitted-block norm implication is exact, but it does not create a
   useful dichotomy: the subset with large omitted norm is precisely trapped
   inside a fixed-cut selector large deviation.  Adaptive principal grounds
   of `A[T]` cannot change that fixed-cut fact.
3. This does **not** falsify the bare arbitrary-cut condition (10.795), direct
   profile-conditioned counting, or non-complement favorable events whose
   threshold is only the permitted `O(n^(3/2-c))` loss rather than complement
   activity at `L_S>=q`.
4. A different row-regularization mechanism cannot repair a fixed-density
   complement column: the mass obstruction applies to every cut once its row
   is at project scale.  A complement successor would have to weaken the
   activity threshold itself or leave compact fixed densities, and would then
   need a new proof that the resulting selectors remain favorable.
