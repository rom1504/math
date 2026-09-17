# Wave 33: soft center kernels and the exact common-core capacity

Status: derivations below are proved algebraically; the finite-dimensional
primal/dual identity and the hard/soft sandwich are also checked in
`tmp/soft_center_r33_check.py`.

## 1. Setup

Use the notation of (10.925)--(10.926).  Thus `Omega` is the finite selector
set, `nu=nu_2` is uniform on the low-row centers `C_2`, and

```math
a_z(S)\in\mathbb Z_{\ge0}
```

is the exact minimum projective distance from `z` to a favorable full
completion for `S`.  Fix an integer batch size `s>=2`, a selector law `w`,
and iid `S_1,...,S_s~w`.  Write

```math
X_z=\sum_{j=1}^s a_z(S_j),\qquad
K_\lambda(z,S)=e^{-\lambda a_z(S)},\qquad
q_{\lambda,w}(z)=\sum_Sw_SK_\lambda(z,S).
```

Independence gives the exact identity

```math
L_\lambda(w):=\mathbb E_{z\sim\nu}\mathbb E[e^{-\lambda X_z}\mid z]
=\mathbb E_\nu q_{\lambda,w}(z)^s.                 \tag{S1}
```

## 2. The lower-tail direction and the sharp subtraction

Let `D=floor(C_D k_0)`.  Since `X_z` is integer, pointwise

```math
\mathbf 1_{\{X_z\le D\}}
\ge
\frac{e^{-\lambda X_z}-e^{-\lambda(D+1)}}
     {1-e^{-\lambda(D+1)}}.                         \tag{S2}
```

For `X_z<=D`, the right side is at most one; for `X_z>=D+1`, it
is nonpositive.  Consequently

```math
\boxed{
\Pr_{z,S_1,\ldots,S_s}\{X_z\le D\}
\ge
\max\left\{0,
\frac{L_\lambda(w)-b_\lambda}{1-b_\lambda}\right\},\qquad
b_\lambda=e^{-\lambda(D+1)}.}
                                                               \tag{S3}
```

The truncation at zero and the subtraction are essential.  The usual
Laplace/Markov implication has the *opposite* direction:

```math
\Pr\{X_z\le D\}\le e^{\lambda D}L_\lambda(w).       \tag{S4}
```

Moreover the positive branch of (S3) is sharp given only `L_lambda`: for
`L_lambda in [b_lambda,1]`, a law supported at `0` and `D+1` attains
equality.  Below `b_lambda`, zero is the best universal lower bound.  Thus no
argument using only this one Laplace moment can remove the tail term.

Taking the infimum over selector laws in (S3) is legitimate because its
right side is increasing in `L_lambda`:

```math
\inf_w\Pr\{X_z\le D\}
\ge
\max\left\{0,
\frac{V_{s,\lambda}^{s}-b_\lambda}{1-b_\lambda}\right\},
\qquad
V_{s,\lambda}:=\inf_{w\in\Delta(\Omega)}
\|K_\lambda w\|_{L^s(\nu)}.                         \tag{S5}
```

## 3. Exact convex/minimax formulation

Put `q=s/(s-1)`.  Positivity lets one restrict the dual of the `L^s` norm
to nonnegative functions.  Finite-dimensional minimax therefore gives

```math
\boxed{
V_{s,\lambda}
=\max_{\substack{h\ge0\\\|h\|_{L^q(\nu)}\le1}}
  \min_{S\in\Omega}
  \mathbb E_{z\sim\nu}\!\left[h(z)e^{-\lambda a_z(S)}\right].}
                                                               \tag{S6}
```

Indeed,

```math
\begin{aligned}
\inf_{w\in\Delta(\Omega)}\|K_\lambda w\|_s
&=\inf_w\sup_{h\ge0,\ \|h\|_q\le1}
  \mathbb E_\nu[h(K_\lambda w)]\\
&=\sup_h\inf_w\sum_Sw_S\mathbb E_\nu[hK_{\lambda,S}]\\
&=\sup_h\min_S\mathbb E_\nu[hK_{\lambda,S}].
\end{aligned}
```

This is an exact fractional common-center capacity.  At a primal minimizer
`w_*`, put `u=K_lambda w_*` and `V=||u||_s`.  When `V>0`, an optimal dual
witness is

```math
h_*(z)=u(z)^{s-1}/V^{s-1}.                            \tag{S7}
```

The KKT conditions say

```math
\mathbb E_\nu[h_*K_{\lambda,S}]\ge V
```

for every selector, with equality for every `S` in the support of `w_*`.

There is also an exact Renyi decomposition.  Let
`Z_w=E_nu q_{lambda,w}` and, when `Z_w>0`, let `P_w` have density
`q_{lambda,w}/Z_w` relative to `nu`.  Then

```math
L_\lambda(w)
=Z_w^s\exp\{(s-1)D_s(P_w\|\nu)\}.                    \tag{S8}
```

Thus (S5)--(S6), rather than a first moment alone, is the relevant
order-`s` soft-cover functional.

## 4. A strictly weaker signing-specific certificate

Equations (S5)--(S6) isolate the following sufficient statement about the
*actual favorable fibers of the signing*.  For fixed constants `A<Lambda`,
take

```math
\lambda=\Lambda\frac{rL_0}{D+1}.
```

It suffices to construct one nonnegative center weight `h` (depending on the
target signing and target pair, but not on `S` or `w`) such that

```math
\|h\|_{s/(s-1)}\le1,
\qquad
\min_S\mathbb E_{\nu_2}
   [h(z)e^{-\lambda a_z(S)}]\ge e^{-A TL_0}.          \tag{S9}
```

Since `s=ceil(r/T)`, exactly

```math
r\le sT<r+T.
```

Thus (S9) gives `V^s>=e^{-A sT L_0}`, whereas the subtraction is
`b_lambda=e^{-Lambda rL_0}`.  The exact tail-margin condition is

```math
\Lambda r>A sT.                                      \tag{S9a}
```

For fixed `Lambda>A`, it holds eventually because `T/r=o(1)`.  In fact the
ratio of the subtraction to the main term then tends to zero exponentially.
Equation (S5) gives an `e^{-O(rL_0)}` hard-event lower bound, and (10.926)
proves the partial collision.

This common-core certificate is genuinely weaker than uniform normalized
degrees.  As a hard-kernel toy (or the `lambda->infinity` limit when the
cost is zero exactly on `C`), suppose all columns are the indicator of one
shared center set `C` of `nu`-mass `delta`.  Then

```math
\min_S\mathbb E K_S=\delta,
\qquad
V_{s,\lambda}=\delta^{1/s},
\qquad
h=\delta^{-(s-1)/s}\mathbf1_C.                       \tag{S10}
```

For `delta=e^{-A rL_0}`, the normalized-degree/Jensen route sees only
`e^{-A rL_0}` per selector and loses a factor `s` in the exponent, while
(S6) sees

```math
\delta^{1/s}=e^{-A(r/s)L_0}\ge e^{-A TL_0}
```

and pays for the shared core only once.  Since `T/r=o(1)`, this is
`e^{-(A+o(1))TL_0}`.  This is precisely the cross-selector correlation
missing from (10.928).  The exact soft kernel has small nonzero tails rather
than indicator columns; the toy demonstrates strict separation of the
high-power capacity from normalized degrees, not an additional claim about
the geometry of `a_z(S)`.

However, (S9) is not currently a theorem about minimizers.  It is a clean
new target whose proof must use global signing minimality to force a common
weighted set of low-row centers.  The trace bound `|C_2|>=2^{n-1}` alone
does not imply it.

## 5. The soft scale is nearly hard

Tail subtraction forces

```math
\lambda(D+1)>A sT L_0,
\qquad r\le sT<r+T,
\qquad\text{hence}\qquad
\lambda=\Omega(rL_0/D)=\Omega(r\log n).               \tag{S11}
```

At the Wave 32 scales this is

```math
\lambda=\Omega(n^{1/4+c_0}(\log n)^2).
```

For the equal-share radius `d=D/s`, one has

```math
\lambda d=\Theta(TL_0).                               \tag{S12}
```

Thus the kernel has to resolve essentially the same tiny favorable
neighborhood as the hard moment.  A bound on `L_lambda` without a strict
constant advantage over `e^{-lambda(D+1)}` proves nothing.  Softening is
therefore a useful convex dualization, not by itself a geometric smoothing
mechanism.

The precise one-kernel falsifier is

```math
V_{s,\lambda}^{s}\le e^{-\lambda(D+1)},             \tag{S13}
```

which makes the sharp lower bound (S5) equal to zero.  At the project level,
falsifying an adaptively chosen kernel requires a sequence of target pairs
such that, on each target pair, **every** admissible `lambda` capable of a
tail margin fails uniformly, for example, with `A_n` denoting the admissible
set,

```math
\sup_{\lambda\in\mathcal A_n}V_{s,\lambda}
=e^{-\omega(TL_0)}.                                  \tag{S14}
```

Such a sequence prevents an adaptive choice of `lambda` from reaching the
required exponent.  The weaker quantifier "for every `lambda` there exists
a bad target pair" would not suffice, because the kernel may be selected
after seeing the target pair.

## 6. Independent-label wall survives the exact dual

Use the deterministic independent-label system underlying (10.930), with
uniform selector law `w_U`.  Fix a minimum admissible kernel scale

```math
\lambda_{\min}=\Lambda_{\min}\frac{rL_0}{D+1},
\qquad
d_{\max}=\left\lfloor\frac{\alpha n}{2\lambda_{\min}}\right\rfloor.
```

Here `lambda_min->infinity`, so `d_max=o(n)`.  Apply the Chernoff construction
once at this largest radius.  It supplies one deterministic assignment
(depending on `lambda_min`, but not on the later adaptive choice of
`lambda>=lambda_min`) such that, for a fixed `alpha>0`,

```math
\max_z w_U\{S:a_z(S)\le d_{\max}\}\le e^{-\alpha n}. \tag{S15}
```

For integer distances,

```math
e^{-\lambda a_z(S)}
\le\mathbf1_{\{a_z(S)\le d\}}+e^{-\lambda(d+1)}.
```

For every `lambda>=lambda_min`, use `d=d_max`.  Then
`lambda(d_max+1)>=lambda_min(d_max+1)>alpha n/2`.
Equations (S15) and the last display give, simultaneously for every such
`lambda`,

```math
q_{\lambda,w_U}(z)\le e^{-\alpha n}+e^{-\alpha n/2}
\le2e^{-\alpha n/2},
```

uniformly in `z`, and hence

```math
\boxed{V_{s,\lambda}\le2e^{-\alpha n/2}
=e^{-\Omega(n)}=e^{-\omega(TL_0)}.}                  \tag{S16}
```

Here `TL_0=o(n)` because `T<=n^eta` and `eta<c_0`.  Thus one abstract
independent-label assignment defeats every adaptive kernel choice
`lambda>=lambda_min`; even the exact weighted dual certificate fails.  No
choice of a concentrated dual weight `h` rescues it.  As before, those labels
need not arise as favorable fibers of one globally minimizing signing, so
(S16) identifies the indispensable signing-specific input rather than
falsifying (S9).

## 7. Surviving target

The weakest exact condition visible to a single soft kernel is the capacity
margin

```math
\boxed{
V_{s,\lambda}^{s}\ge e^{-O(rL_0)}
\quad\text{and}\quad
e^{-\lambda(D+1)}=o(V_{s,\lambda}^{s}).}            \tag{S17}
```

By (S5) it is sufficient, and the subtraction in (S3) is information-
theoretically sharp given only that transform.  Formula (S6) converts (S17)
exactly into a common weighted-center theorem.  This formulation is much
weaker than uniform normalized degrees, but the independent-label wall and
the large required `lambda` show that only a new global-minimality mechanism
can establish it.
