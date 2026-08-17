# Absolute-overlap shells give low-cap physical response separation

Status: **PROVED AND INDEPENDENTLY AUDITED; PASS after repairs**.

This note isolates the projective datum missing from the signed first-moment
balance of Theorems 36.2 and 36.5.  It also records a simpler physical
compiler than full pinning: a repeated rank-one shore of width `Theta(sqrt n)`.
The conclusion is a contextual packing conditional on projective shell
packing.  No such packing is proved here for exact minimizers.

## 1. Notation

Let `E=binom(n,2)` and

```math
H_a(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad Q(a)=\max_x|H_a(x)|.
```

Write

```math
c(x)=(x_ix_j)_{i<j},
\qquad z^u=\sigma_uc(u),\quad \sigma_u\in\{\pm1\}.
```

Thus `Q(a)=max_z <a,z>`, where `z` ranges over the augmented cuts.  A family
`U` is `(d,gamma)`-projectively separated near the top if

```math
\langle a,z^u\rangle\ge Q(a)-d
\quad(u\in U),
\qquad
{ |\langle z^u,z^v\rangle|\over E}\le1-\gamma
\quad(u\ne v).                                      \tag{AO.1}
```

The absolute value in (AO.1) is essential because the physical parent cap
optimizes the global orientation.

## 2. A spherical two-peak inequality

### Lemma AO.1

Let `U,V,Y` be unit vectors, `c=|U dot V|`, and suppose
`c^2<=1-theta`, where `0<theta<=1`.  For all `alpha,lambda>0`,

```math
{\alpha\over2}(V\mathbin\cdot Y)^2
+\lambda|U\mathbin\cdot Y|
\le {\alpha\over2}+\lambda
-{\theta\over4}\min\{\alpha,\lambda\}.           \tag{AO.2}
```

#### Proof

Put `r=|V dot Y|` and `s=|U dot Y|`.  The response deficit is at least

```math
{1\over2}\min\{\alpha,\lambda\}
\{(1-r^2)+(1-s^2)\},
```

because `1-s >= (1-s^2)/2`.  The largest eigenvalue of the sum of the two
rank-one projectors `VV^T+UU^T` is `1+c`.  Hence the expression in braces
is at least `1-c`, and

```math
1-c={1-c^2\over1+c}\ge{\theta\over2}.
```

This proves (AO.2). `square`

For augmented cuts one has the exact identity

```math
\langle z^u,z^v\rangle
=\sigma_u\sigma_v{(u\mathbin\cdot v)^2-n\over2}.     \tag{AO.3}
```

Consequently (AO.1) implies, for `U=u/sqrt n,V=v/sqrt n`,

```math
1-|U\mathbin\cdot V|^2\ge\gamma(1-1/n).             \tag{AO.4}
```

## 3. Sparse-flip children

Fix `0<alpha<1` and put `p=alpha/sqrt n`.  For each `u in U`, let

```math
D_u=\{e:a_ez^u_e=-1\}.
```

Independently flip every edge in `D_u` with probability `p`, obtaining an
exact signing `b^u`.  For every augmented cut `z`,

```math
\mathbb E\langle b^u,z\rangle
=(1-p)\langle a,z\rangle+p\langle z^u,z\rangle.      \tag{AO.5}
```

Indeed `2 1_(D_u)z^u=z^u-a` coordinatewise.  Bernstein's inequality and a
union bound over at most `2^(n+1)` children and `2^n` augmented cuts give a
simultaneous realization satisfying

```math
\left|\langle b^u,z\rangle-
 \{(1-p)\langle a,z\rangle+p\langle z^u,z\rangle\}
\right|
\le \rho_n,
\qquad
\rho_n=C(\sqrt\alpha n^{5/4}+n),                       \tag{AO.6}
```

for one absolute `C`.  The same event can be intersected with

```math
|F_u|\le 2pE,
\qquad
Q(b^u)\le Q(a)+4pE\le Q(a)+2\alpha n^{3/2}.            \tag{AO.7}
```

The exponent in the Bernstein tail is `Omega(n)` after increasing `C`, so
the stated union bound is legitimate even for the complete augmented-cut
family.  If `alpha=alpha_n` varies, the same simultaneous statement for an
exponential family is asserted only when `alpha_n sqrt n->infinity`.
Formula (AO.7) is intentionally coarse.

## 4. The exact-sign free-shore compiler

For a query `u`, let `h=floor(lambda sqrt n)` and take the exact-sign cross
block

```math
B^u=u\mathbf1_h^T\in\{\pm1\}^{n\times h}.             \tag{AO.8}
```

Complete the shore by any common hollow signing `C_h`, and define

```math
P^{v|u}=
\begin{pmatrix}b^v&B^u\\(B^u)^T&C_h\end{pmatrix}.     \tag{AO.9}
```

Introduce the projective trust response

```math
T_b(u,h)=\max_{y,\tau}
 \{\langle b,\tau c(y)\rangle+h|u\mathbin\cdot y|\}. \tag{AO.10}
```

The shore spins are completely free.  Optimizing them when `C_h` is omitted
gives (AO.10) exactly, and therefore

```math
|Q(P^{v|u})-T_(b^v)(u,h)|\le Q(C_h)
\le {h(h-1)\over2}=O(n).                              \tag{AO.11}
```

Thus the shore has target-scale field strength but only subleading internal
cap.  Unlike a multi-port raw frame, every endpoint in (AO.8) lies on the
single scalar segment `[-hu,hu]`; no endpoint-balance theorem is needed.

## 5. Conditional physical packing theorem

### Theorem AO.2

Suppose `U` obeys (AO.1), let `b^u` be one simultaneous realization from
(AO.6)--(AO.7), and use the common query family (AO.8)--(AO.9).  Put

```math
\theta_n=\gamma(1-1/n),
\qquad
\Delta_n={\theta_n\over4}\min\{\alpha,h/\sqrt n\}. \tag{AO.12}
```

Then for every ordered pair `u ne v`,

```math
Q(P^{u|u})-Q(P^{v|u})
\ge
\left(\Delta_n-{\alpha\over n}\right)n^{3/2}
-(1-p)d-2\rho_n-2Q(C_h).                           \tag{AO.13}
```

For one fixed ordered comparison, the same bound requires the near-top
condition only for the target `u`; the decoy `v` may be an arbitrary
augmented cut satisfying the overlap condition.  The symmetric hypothesis
on all members of `U` is used to make every member targetable and hence to
obtain a packing.

In particular, if `alpha,lambda,gamma` are fixed, `d=o(n^(3/2))`, and
`|U|=K_n`, then

```math
Q(P^{u|u})-Q(P^{v|u})
\ge
\left[
 {\gamma\over4}\min\{\alpha,\lambda\}-o(1)
\right]n^{3/2}.                                     \tag{AO.14}
```

All matrices in (AO.9) are hollow exact signings of order
`n+O(sqrt n)` and have cap `O(n^(3/2))` whenever `Q(a)=O(n^(3/2))`.
Hence response accuracy below half the gap on this declared all-spins-free
query bank requires at least `log_2 K_n` bits.

A law version is immediate.  If a shell law has
`E|R(Z,Z')|<=1-gamma`, its support contains an ordered pair meeting (AO.1).
More quantitatively, if every projective overlap ball

```math
\{z':|R(z,z')|>1-\gamma\}
```

has law mass at most `beta`, greedy deletion yields a
`floor(1/beta)`-point `(d,gamma)` packing and hence that many physical
contextual states.  This ball-mass hypothesis, unlike signed barycentre
control, directly measures the missing projective geometry.

#### Proof

At the target state `(u,sigma_u)`, (AO.5)--(AO.6) give

```math
T_(b^u)(u,h)
\ge(1-p)(Q(a)-d)+pE+hn-\rho_n.                    \tag{AO.15}
```

For a cross child and `z=tau c(y)`, use `<a,z><=Q(a)` and

```math
|\langle z^v,z\rangle|
\le{(v\mathbin\cdot y)^2+n\over2}
```

to obtain

```math
T_(b^v)(u,h)
\le(1-p)Q(a)+\rho_n+{pn\over2}
+\max_y\left\{{p\over2}(v\mathbin\cdot y)^2
                +h|u\mathbin\cdot y|\right\}.       \tag{AO.16}
```

Relax `y/sqrt n` to the unit sphere.  Equations (AO.4) and (AO.2) bound the
last maximum by

```math
n^{3/2}\left\{{\alpha\over2}+{h\over\sqrt n}-\Delta_n\right\}.
```

Subtract (AO.16) from (AO.15).  Since

```math
pE-{pn\over2}-{\alpha\over2}n^{3/2}
=-\alpha n^{1/2},
```

this gives the trust-response version of (AO.13).  Equation (AO.11) costs
at most `2Q(C_h)` and proves the displayed bound.  Finally

```math
Q(P^{v|u})\le Q(b^v)+nh+Q(C_h)=O(n^{3/2}),
```

which proves the cap and order claims. `square`

The convenient balanced choice `lambda=alpha` gives leading gap
`alpha gamma/4` in (AO.14).

There is also a uniform vanishing-scale form.  Let `alpha=alpha_n->0`, take
`lambda_n=alpha_n`, and suppose

```math
\alpha_n\sqrt n\longrightarrow\infty,
\qquad d=o(\alpha_n n^{3/2}).                       \tag{AO.14a}
```

Then `h=floor(alpha_n sqrt n)` satisfies `h/sqrt n=(1-o(1))alpha_n`, while

```math
{\rho_n\over\alpha_n n^{3/2}}
=O\left({1\over\sqrt{\alpha_n}n^{1/4}}
        +{1\over\alpha_n\sqrt n}\right)=o(1).       \tag{AO.14b}
```

Consequently the right side of (AO.14) becomes

```math
\left({\gamma\over4}-o(1)\right)
\alpha_n n^{3/2}.                                    \tag{AO.14c}
```

The simultaneous union bound still covers an exponentially large `U`, since
its required deviation is exactly the `rho_n` in (AO.6).

The ordered comparison in the proof is strictly asymmetric: (AO.15) needs
only the **target** `u` to satisfy `<a,z^u> >= Q(a)-d`.  No energy assumption
on the decoy direction `v` is used in (AO.16).  Requiring every member of
`U` to be near-top is necessary only when every child must be the target of
its own query, as in the `K_n`-state packing.

### Corollary AO.3 (one unconditional physical bit in every bounded-cap halo)

Let `a` be any signing with `Q(a)=O(n^(3/2))`, and choose a ground augmented
cut `z^u`.  Choose a Boolean `v` with `|u dot v|<=1` and either orientation
for `z^v`.  Then

```math
{|\langle z^u,z^v\rangle|\over E}
\le {1\over n-1}.
```

For every fixed `0<alpha<1`, apply the construction with `lambda=alpha`.
There are two exact-sign children in the
`2 alpha n^(3/2)` edit-Lipschitz halo of `a` and an all-spins-free common
context of width `Theta(sqrt n)` such that

```math
Q(P^{u|u})-Q(P^{v|u})
\ge\left({\alpha\over4}-o(1)\right)n^{3/2}.        \tag{AO.21}
```

Thus a genuine exact minimizer always has one physical hidden response bit
at total parent cap `O(n^(3/2))`.  This does not give a growing packing:
to target `b^v` in a reverse or third query, `z^v` itself would need to be
near-top.  The absolute-overlap shell hypothesis is precisely what upgrades
this unconditional ordered bit to a reusable family.

The corollary also holds in a vanishing halo.  If
`alpha_n->0` and `alpha_n sqrt n->infinity`, both children of an exact base
lie in `N_n(2alpha_n)` and

```math
Q(P^{u|u})-Q(P^{v|u})
\ge\left({\alpha_n\over4}-o(\alpha_n)\right)n^{3/2}.
                                                               \tag{AO.21d}
```

Thus the physical bit survives down to every halo asymptotically wider than
the concentration threshold `n^(-1/2)`; its normalized gap vanishes with the
halo width, as it must in this construction.

## 6. What fractional balance does and does not imply

The signed first moment in Theorem 36.2 does **not** imply the absolute
overlap premise (AO.1).  Let `u=mathbf1`, let `v` differ from `u` on
`k=ceil(sqrt n)` coordinates, and put

```math
z=c(u),\qquad z'=-c(v),\qquad
\mu={1\over2}(\delta_z+\delta_(z')).              \tag{AO.17}
```

The two words agree exactly on the `k(n-k)` crossing edges.  Thus

```math
{1\over E}\sum_e|\mathbb E_\mu Z_e|
={k(n-k)\over E}=O(n^{-1/2}),
```

and for independent `Z,Z'~mu`,

```math
\mathbb E R(Z,Z')={k(n-k)\over E}=O(n^{-1/2}).       \tag{AO.18}
```

Nevertheless every pair in the support has absolute normalized overlap

```math
1\quad\hbox{or}\quad
\left|1-{2k(n-k)\over E}\right|=1-O(n^{-1/2}).      \tag{AO.19}
```

This is a countermodel to the **marginal implication only**.  It is not
asserted to be the positive thin shell of an exact minimizer.  Ruling out
precisely this nearly orientation-antipodal two-cap geometry is the new
near-minimizer-specific obligation.

Signed overlap plus positivity gives only the scale-sharp weak exclusion

```math
R(z,z')\ge-1+{2(Q(a)-d)\over E},                    \tag{AO.20}
```

because `<a,z>+<a,z'><=2 #\{e:z_e=z'_e\}=E(1+R)`.
For `Q(a)=Theta(n^(3/2))`, (AO.20) permits
`R=-1+Theta(n^(-1/2))`, exactly the regime in (AO.17)--(AO.19).

There is a useful exact dichotomy behind that countermodel.  Suppose a law
`mu` on the positive `d`-shell has

```math
E^{-1}\sum_e|\mathbb E_\mu Z_e|\le\delta              \tag{AO.21a}
```

and contains no pair with `|R(z,z')|<=1-gamma`.  Fix `z_0` in its support
and give every `z` the sign `s(z)` for which
`R(z,z_0)` has that sign.  Then

```math
|\mathbb E_\mu s(Z)|\le\delta+\gamma.                 \tag{AO.21b}
```

Indeed `|R(z,z_0)-s(z)|<gamma`, while the absolute value of the mean of
`R(Z,z_0)` is at most the left side of (AO.21a).  Thus, when
`delta+gamma<1`, both signs occur: the shell is forced into two approximately
balanced projective caps around `z_0` and `-z_0`.  Choosing one atom from
each cap gives `R(z_+,z_-)<-1+2gamma`.  Positivity of both shell atoms and
(AO.20) then force

```math
\boxed{\gamma>{Q(a)-d\over E}.}                          \tag{AO.21c}
```

Therefore fractional balance yields a precise alternative, but only at the
`Q/E=Theta(n^(-1/2))` scale: either there is projective separation, or there
is a balanced pair of nearly orientation-antipodal caps.  At that forced
scale Theorem AO.2 supplies only an `O(n)` response gap, not the required
`Theta(n^(3/2))` gap.  A genuine target-scale result must rule out the
two-cap alternative at a fixed projective scale (or find a different
compiler that exploits it).

## 7. Scope and assumption audit

1. **What is proved.**  Every bounded-cap signing has the one ordered
   physical collision in Corollary AO.3.  More generally, a constant
   projective packing in a vanishing-width
   positive shell compiles into the same-size physical contextual packing,
   with exact signs, all spins free, order `n+O(sqrt n)`, total cap
   `O(n^(3/2))`, and a fixed target-scale gap.
2. **What is not proved.**  Theorems 36.2 and 36.5 do not provide even three,
   or necessarily two, projectively separated shell atoms.  No growing
   physical packing for actual minimizers follows here.
3. **Information content.**  Given the witnesses, checking (AO.1) is
   polynomial-size and retains only their energies and pairwise projective
   overlaps, far less than the full Boolean landscape.  Proving that every
   exact minimizer supplies such a family may still be hard; this note does
   not assign it a positive assumption-distance score without an independent
   structural implication.
4. **No independent scalar channels.**  Each context uses one rank-one field
   segment and the comparison is made only after maximizing the complete
   child-plus-field response.  The spherical inequality preserves the
   tradeoff between the sparse-flip peak and the field peak before taking an
   absolute value.
5. **No convergence claim.**  Even an exponential conditional packing would
   be an incompressibility theorem for a declared query class, not a recurrence
   or an all-order transfer theorem for `M_n`.
