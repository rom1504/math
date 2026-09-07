# Wave 46B: normalized common-core escape and the exact coarea target

The statements explicitly marked **Verified** below were derived from the
definitions and checked independently in
`tmp/common_core_excess_r46_check.py`.  The minimizer-specific boundary
estimate is **Open**.  The scalable cylinder is an abstract slice family, not
an asserted exact-minimizer incidence family.

## 1. Incidence-size-biased escape is an exact boundary ratio

Let `Omega=binom([n],m)` have its uniform law, let `nu` be the fixed
selector-independent project-row center law, and put

```math
f_z(S)=\mathbf 1\{D_z(S)\le H\},\qquad
D_z(S)=Q(A[S])-|z_S^{\mathsf T}A[S]z_S|,
\qquad a_z=\langle f_z,1\rangle.
```

Write `J_1=E_z a_z>0`.  For the down-up kernel `K_ell`, set

```math
B_\ell(z)=\langle f_z,(I-K_\ell)f_z\rangle,
\qquad B_\ell=\mathbb E_z B_\ell(z),
\qquad \epsilon_\ell=B_\ell/J_1.
```

Self-adjointness and stochasticity give the **Verified exact identities**

```math
\begin{aligned}
B_\ell(z)
 &=\mathbb E_{S,T}f_z(S)(1-f_z(T))\\
 &=\frac12\mathbb E_{S,T}(f_z(S)-f_z(T))^2,\\
\rho_\ell
 &=\frac{\mathbb E_z\langle f_z,K_\ell f_z\rangle}{J_1}
 =1-\epsilon_\ell.
\end{aligned}
\tag{R46.1}
```

Here `S` is uniform and `T` is a `K_ell` child.  Thus `epsilon_ell` is
literally the probability of escape after first size-biasing `(z,S)` by the
hard incidence.  If

```math
b_z(R)=\Pr\{f_z(S)=1\mid R\subset S\},\qquad |R|=\ell,
```

then the common-core factorization also gives

```math
\boxed{
B_\ell=\mathbb E_{z,R}b_z(R)(1-b_z(R)),\qquad
\epsilon_\ell=
\frac{\mathbb E_{z,R}b_z(R)(1-b_z(R))}{J_1}.}
\tag{R46.2}
```

This is the requested incidence-size-biased boundary identity.  It includes
empty selector pairs without any convention.

## 2. Full-spectrum extraction and the missing factor in Wave 45

Let `P_j` be the Johnson harmonic projections and
`W_j(z)=||P_jf_z||_2^2`.  The eigenvalues are

```math
\lambda_j(\ell)=
\frac{(\ell)_j(n-m)_j}{(m)_j(n-\ell)_j},
\qquad \lambda_0=1,
```

with zero eigenvalues after the numerator vanishes.  Put

```math
\lambda=\lambda_1(\ell),\qquad
\delta=1-\lambda,\qquad
\Xi_\ell=
\mathbb E_z\sum_{j\ge2}(\lambda-\lambda_j)W_j(z)\ge0,
\qquad J_2=\mathbb E_z a_z^2.
```

Since `sum_(j>=1)W_j=a_z-a_z^2`, the complete spectrum yields the
**Verified exact accounting**

```math
\boxed{
B_\ell=\delta(J_1-J_2)+\Xi_\ell,
\qquad
\frac{J_2}{J_1}
=\frac{\rho_\ell-\lambda+\Xi_\ell/J_1}{1-\lambda}.}
\tag{R46.3}
```

Consequently

```math
\boxed{
\max_z a_z\ge\frac{J_2}{J_1}
\ge\frac{\rho_\ell-\lambda}{1-\lambda}
=1-\frac{\epsilon_\ell}{1-\lambda}.}
\tag{R46.4}
```

The bound `(rho-lambda)/(1-lambda)` is the correct one-step extraction; it
strictly sharpens `rho-lambda` in (10.1176).  It also explains exactly what
the higher Johnson levels can add: a certified lower bound on `Xi_ell`
improves (R46.4) by `Xi_ell/[J_1(1-lambda)]`.  Without such information the
basic bound is sharp.  For example, a one-coordinate star has only constant
and level-one components and attains equality.  The `A_6` hard fibers also
have `Xi_ell=0` at every audited scale.

For a desired hard degree `eta`, the clean sufficient criterion is therefore

```math
\boxed{
B_\ell\le(1-\lambda)(1-\eta)J_1,
\quad\text{equivalently}\quad
\rho_\ell-\lambda\ge(1-\lambda)\eta.}
\tag{R46.5}
```

Taking `eta=exp{-C L_0}` proves the required
`max_z a_z>=exp{-C L_0}`.  The universal spectral inequality points in the
opposite direction,

```math
B_\ell\ge(1-\lambda)(J_1-J_2),
```

so (R46.5) asks the favorable fibers to have almost the least boundary
compatible with their second moment; it is not a generic expansion estimate.

## 3. The near-core scale needs no powers

For `ell=m-s`, exactly

```math
\delta_s=1-\lambda_1(m-s)=\frac{sn}{m(n-m+s)}.
\tag{R46.6}
```

If `m/n -> p in (0,1)` and `s=(alpha+o(1))n/L_0`, then

```math
\delta_s=
\left(\frac{\alpha}{p(1-p)}+o(1)\right)L_0^{-1}.
\tag{R46.7}
```

Thus if the incidence-biased escape has the expansion

```math
\epsilon_{m-s}=\frac{\beta+o(1)}{L_0},
\qquad
\alpha_0=\frac{\alpha}{p(1-p)},
```

then (R46.4), at `t=1`, already gives

```math
\max_z a_z\ge1-\frac{\beta}{\alpha_0}+o(1).
\tag{R46.8}
```

In particular a strict rate gap `beta<alpha_0` gives a constant degree, not
merely a premise that powers can amplify.  For the actual project target it
is enough to prove the much finer one-sided estimate

```math
\epsilon_{m-s}
\le\delta_s\{1-\exp(-C L_0)\}.
\tag{R46.9}
```

**Falsified as a mechanism:** powers of `K_ell` are unnecessary once there is
a rate gap and cannot create one when there is not.  The substantive issue is
the single-step normalized boundary (R46.9).

## 4. An exact clipped-deficit coarea sufficient inequality

The boundary at one threshold can be replaced by a precise layer-cake target.
For `0<=h<=H`, define

```math
f_{z,h}(S)=\mathbf1\{D_z(S)\le h\},\quad
J(h)=\mathbb E_z\mathbb E_S f_{z,h}(S),\quad
B_\ell(h)=\mathbb E_z\langle f_{z,h},(I-K_\ell)f_{z,h}\rangle,
```

and let the clipped favorable slack be

```math
u_z(S)=(H-D_z(S))_+.
```

Because `D_z(S)>=0`, Fubini and reversibility give the **Verified exact
coarea identities**

```math
\boxed{
\int_0^H J(h)\,dh=\mathbb E_{z,S}u_z(S),
\qquad
\int_0^H B_\ell(h)\,dh
=\frac12\mathbb E_{z,S,T}|u_z(S)-u_z(T)|.}
\tag{R46.10}
```

Therefore the concrete inequality

```math
\boxed{
\frac12\mathbb E_{z,S,T}|u_z(S)-u_z(T)|
\le\delta_s(1-\eta)\mathbb E_{z,S}u_z(S)}
\tag{R46.11}
```

would imply that some `h<=H` has
`B_(m-s)(h)<=delta_s(1-eta)J(h)`.  Applying (R46.4) to that level gives one
project-row center with

```math
U_m\{S:D_z(S)\le h\}\ge\eta,
```

and hence the same lower bound at the original cap `H`.  With
`eta=exp{-C L_0}`, (R46.11) proves the hard exceptional-center target.
This is an exact sufficient theorem, not a consequence currently known from
minimality.  It also makes the difficulty visible: it asks the clipped
deficit to be nearly `L^1`-invariant at the Johnson spectral-gap scale.

## 5. What exact minimality and crude replacement smoothness really give

Exact minimality supplies only scoped inputs here.  Principal monotonicity
and the definition of `q_m` give

```math
q_m\le Q(A[S])\le q_n,\qquad D_z(S)\ge0.
\tag{R46.12}
```

If `|S|=|T|=m`, `r=|S\setminus T|=|T\setminus S|`, and `C=S cap T`, set

```math
L_r=2r(m-r)+r(r-1).
```

Extending a ground of `A[C]` randomly to either branch proves
`Q(A[C])<=Q(A[S]),Q(A[T])`, while the at most
`r(m-r)+binom(r,2)` new edges prove

```math
|Q(A[S])-Q(A[T])|\le L_r.
```

For the fixed-center energy, each of the two branches contributes at most
`L_r`, and therefore

```math
\boxed{
\big||z_S^TA[S]z_S|-|z_T^TA[T]z_T|\big|\le2L_r,
\qquad |D_z(S)-D_z(T)|\le3L_r\le6rn.}
\tag{R46.13}
```

Thus a selector favorable at threshold `H_0` keeps all selectors within
Johnson radius `s` favorable only after thickening the threshold by at most
`3L_s<=6sn`.  With `s=n/L_0`, this is `O(n^2/L_0)`.  It is `o(H)` for the
macroscopic cap `H=Theta(n^(3/2))` whenever `c<1/4`, but it is `o(t)` for the
actual saving tolerance `t=Theta(n^(3/2-c))` only when `c<1/8`.

This does not produce the required degree.  One Johnson ball of radius
`r=o(n)` has

```math
\frac{\sum_{j\le r}\binom mj\binom{n-m}j}{\binom nm}
=\exp\{-\Theta(n)+O(r\log(n/r))\}.
\tag{R46.14}
```

Even spending the full tolerance allows total radius only
`O(t/n)=O(n^(1/2-c))=o(n)`, so the mass remains `exp{-Theta(n)}`, much less
than `exp{-O(L_0)}`.  Iterating for the `Theta(L_0)` spectral mixing time
would instead accumulate `Theta(n^2)` crude slack.  **Falsified:** a proof
based only on (R46.13) and Johnson-ball growth.  This is the selector-space
version of the existing Lipschitz-ball wall, not a new composition
mechanism.

There is a sharper generic `L^2` comparison, but it still does not prove
(R46.11).  The Johnson eigenvalues obey the elementary factorial comparison

```math
1-\lambda_j(m-s)\le s\{1-\lambda_j(m-1)\}.
```

Combining it with the verified adjacent principal-cap estimate (10.1135)
gives

```math
\mathbb E_{S,T}(Q(A[S])-Q(A[T]))^2
\le8s\,\mathbb E_SQ(A[S]).
\tag{R46.15}
```

The fixed-center energy is a degree-two slice polynomial, so

```math
\mathbb E_{S,T}(C_z(S)-C_z(T))^2
\le2\{1-\lambda_2(m-s)\}\operatorname{Var}_{U_m}(C_z).
\tag{R46.16}
```

At project row,
`Var(C_z)=O(n^(9/4-c)+n^2)` by (10.1123).  Hence, for
`s=Theta(n/L_0)`, the signed-energy contribution to the left of (R46.11)
is `O(n^(3/4))`, but the cap contribution available from (R46.15) is only

```math
O(\sqrt{s n^{3/2}})=O(n^{7/8+c/2}).
```

Even in the optimistic case `E u=Theta(H)`, the right side of (R46.11) is
only `O(delta_s H)=O(n^(3/4+c))`.  The cap estimate misses by
`n^(1/8-c/2)` for every `c<1/4`; a smaller actual slack mean worsens it.
This is a limitation of the known bounds, not a counterexample to
(R46.11).  Minimality would have to align principal-cap motion with the
fixed-center field or directly suppress clipped boundary.

## 6. Nonnegative mixtures do not improve extraction

Let `K_w=sum_ell w_ell K_ell`, where `w_ell>=0` and `sum w_ell=1`.  All
Johnson kernels diagonalize in the same harmonic decomposition, so

```math
\lambda_w=\sum_\ell w_\ell\lambda_\ell,
\qquad \rho_w=\sum_\ell w_\ell\rho_\ell,
\qquad \delta_w=\sum_\ell w_\ell(1-\lambda_\ell).
```

Consequently the normalized excess is exactly

```math
\boxed{
\frac{\rho_w-\lambda_w}{1-\lambda_w}
=\sum_\ell
\frac{w_\ell(1-\lambda_\ell)}{\delta_w}
\frac{\rho_\ell-\lambda_\ell}{1-\lambda_\ell}.}
\tag{R46.17}
```

The same formula holds after adding the full-spectrum corrections
`Xi_ell/J_1`.  Thus a nonnegative mixture is a convex average of the
constituent normalized bounds and cannot beat the best single scale.  The
finite LP has a vertex optimizer.  A mixture may package estimates proved at
several scales, but it supplies no new mathematical extraction mechanism.
Adding any mass at the identity has zero effective weight in (R46.17).

## 7. Exact finite tests and a scalable obstruction to genericity

For the `A_6,A_8,A_9` hard fibers, the checker reconstructs every kernel,
the boundary/coarea identity (using cap four), the normalized extraction,
and a rational kernel mixture.  The sequences below are the basic normalized
lower bounds `(rho_ell-lambda_ell)/(1-lambda_ell)` for
`ell=0,...,m-1`:

| instance | `J_2/J_1` | normalized lower bounds |
|:--|:--|:--|
| `A_6,m=5` | `2/3` | `2/3,2/3,2/3,2/3,2/3` |
| `A_8,m=6` | `11/50` | `11/50,99/500,209/1250,121/1000,11/250,-11/100` |
| `A_9,m=7` | `157/864` | `157/864,5275/31752,989/6804,293/2520,659/9072,-1/13608,-55/378` |

The missing mass is exactly the nonnegative `Xi_ell` term in (R46.3).  The
nearest-core basic criterion is negative on `A_8,A_9`, while the true
`J_2/J_1` is positive.  This confirms that nearest-core escape without
high-level information can be strictly misleading.

There is also a rigorous scalable warning.  Fix `p in (0,1)`, let
`m=pn+o(n)`, let `L->infinity` with `L=o(n)`, choose a fixed set `W` of
size `r=(kappa+o(1))L`, and define the slice cylinder

```math
F_W=\{S:W\subset S\}.
```

Its density is

```math
a=\frac{(m)_r}{(n)_r}
=\exp\{\kappa L\log p+o(L)\}=\exp\{-\Theta(L)\}.
\tag{R46.18}
```

Take `s=(alpha+o(1))n/L`.  Conditional on `S in F_W`, if exactly `j`
marked vertices are omitted from the common core, all must be resampled.
Thus the exact retention is

```math
\rho=
\sum_j
\frac{\binom rj\binom{m-r}{s-j}}{\binom ms}
\frac{(s)_j}{(n-m+s)_j}.
\tag{R46.19}
```

The probability no marked vertex is omitted tends
`exp{-kappa alpha/p}`.  The expected number omitted stays bounded, while a
given omitted vertex is resampled with probability `O(1/L)`, so all terms
with `j>=1` contribute `o(1)`.  Therefore

```math
\rho\longrightarrow e^{-\kappa\alpha/p}<1,
\qquad
\lambda_1=1-\frac{\alpha+o(1)}{p(1-p)L}.
\tag{R46.20}
```

The spectral excess is eventually negative even though the family has
exactly the desired `exp{-Theta(L)}` mass.  The checker evaluates (R46.19)
with exact rationals along `n=L^4`, `p=1/2`, `r=L`, `s=n/L` and observes
retention tending `e^{-2}` while `L(1-lambda)->4`.

This is **not** a minimizer counterexample.  It rigorously falsifies the idea
that project-mass alone, generic slice concentration, or a low-complexity
description should force near-core excess.  A positive result must use
special minimizer geometry.

## 8. Frontier

- **Verified:** the exact escape/core-variance identities (R46.1)--(R46.2),
  normalized and full-spectrum extraction (R46.3)--(R46.5), near-core
  calibration, coarea reduction (R46.10)--(R46.11), mixture collapse, the
  scoped Lipschitz and `L^2` bounds, and all finite/scalable computations.
- **Falsified:** spectral powers as an amplification mechanism; nonnegative
  kernel mixtures as an improvement over the best scale; pure
  Lipschitz-ball growth; and any generic claim that an
  `exp{-Theta(L_0)}` slice family must have positive near-core excess.
- **Open:** prove for actual exact-minimizer hard fibers either the direct
  normalized boundary estimate (R46.9), the integrated clipped-deficit
  estimate (R46.11), or a lower bound on the high-level correction `Xi_ell`
  strong enough in (R46.3).  Exact minimality presently supplies none of
  these.

