# Wave 53: exact scalar and truncated-energy certificates for the K-profile

## Status and research conclusion

- **Verified:** finite-dimensional infimal-convolution duality gives an exact
  clipped-coordinate formula for `K_{1,2}(beta,sqrt(H))`, with no Holmstedt
  constant.  It implies the sharp universal lower bound

  ```math
  K_{1,2}(\beta,\sqrt H)\ge
  \min\left\{\frac{\lVert\beta\rVert_2^2}
                    {\lVert\beta\rVert_\infty},
              \sqrt H\lVert\beta\rVert_2\right\}.
  ```

- **Verified:** in the signing variables, if `R=r+b_H`, `u=A[T,S]y`,
  `G=||u||_2^2`, and `J=||u||_infty`, then

  ```math
  G\ge R^2/H,\qquad G\ge JR
  ```

  is a constant-one sufficient coefficient certificate for (10.1276).  A
  head/truncated-energy version below is strictly more flexible when a few
  coordinates are exceptional.

- **Falsified as a pointwise implication:** high cross energy by itself, even
  when it passes the necessary energy threshold, does not force this
  K-profile on exact finite minimizers.  An exact stored `A8` state passes
  `G>=r^2/H` but has `K_H(u)<r` already with `b_H=0`.

- **Open:** no known exact-minimizer constraint forces a saved population of
  states satisfying the scalar or truncated certificate.  The result reduces
  (10.1277) to an explicit high-cross plus delocalization abundance theorem;
  it does not prove that theorem or convergence.

## 1. Exact dual and clipped-coordinate formula

For a finite real vector `x` and `t>0`, put

```math
K_{1,2}(x,t)=\inf_{x=a+b}
 \{\lVert a\rVert_1+t\lVert b\rVert_2\}.
```

The exact finite-dimensional dual is

```math
\boxed{
K_{1,2}(x,t)=\max\{\langle x,z\rangle:
 \lVert z\rVert_\infty\le1,\ \lVert z\rVert_2\le t\}.}
\tag{R53.1}
```

For the easy direction, every feasible `z` and every decomposition `x=a+b`
obey

```math
\langle x,z\rangle
\le\lVert a\rVert_1\lVert z\rVert_\infty
 +\lVert b\rVert_2\lVert z\rVert_2
\le\lVert a\rVert_1+t\lVert b\rVert_2.
```

The reverse direction is ordinary finite-dimensional convex duality: the
conjugates of `||.||_1` and `t||.||_2` are the indicators of, respectively,
the unit `l_infty` ball and the radius-`t` `l_2` ball; the conjugate of their
infimal convolution is the sum of those indicators.  There is no closure or
infinite-sequence issue here.

Take `t=sqrt(H)` and write `b_i=|x_i|`.  If the support has at most `H`
coordinates, (R53.1) gives `K_H(x)=||x||_1`.  Otherwise there is a unique
`lambda>0` such that

```math
N_\lambda(x):=\sum_i\min\{1,b_i^2/\lambda^2\}=H.
```

The optimizing dual vector is

```math
z_i=\operatorname{sgn}(x_i)\min\{1,b_i/\lambda\},
```

and hence the **exact clipped formula** is

```math
\boxed{
K_H(x)=J_\lambda(x):=
\sum_i\min\{b_i,b_i^2/\lambda\}.}
\tag{R53.2}
```

Thus Holmstedt is not needed for any constant in this attack.  For any
chosen cutoff with `N_lambda<=H`, the same clipped vector is feasible and
gives the rigorous lower certificate `K_H(x)>=J_lambda(x)`.

There is also a useful head/tail form.  For any coordinate set `I` of size
`s<H`, set

```math
E_{I^c}=\sum_{i\notin I}x_i^2,
\qquad M_{I^c}=\max_{i\notin I}|x_i|.
```

If

```math
M_{I^c}\sqrt{H-s}\le\sqrt{E_{I^c}},
```

then taking signs on `I` and the normalized `l_2` direction on its complement
in (R53.1) proves

```math
\boxed{
K_H(x)\ge \lVert x_I\rVert_1+
 \sqrt{(H-s)E_{I^c}}.}
\tag{R53.3}
```

At the active set from (R53.2), (R53.3) is an equality.  Unlike a global
`l_infty` estimate, it allows a small exceptional head.

## 2. The sharp bound from energy and maximum coordinate

Let

```math
E=\lVert x\rVert_2^2,\qquad M=\lVert x\rVert_\infty.
```

For nonzero `x`, choose

```math
z=\frac{x}{\max\{M,\sqrt{E/H}\}}.
```

It is feasible in (R53.1), so

```math
\boxed{
K_H(x)\ge
\frac{E}{\max\{M,\sqrt{E/H}\}}
=\min\left\{\frac EM,\sqrt{HE}\right\}.}
\tag{R53.4}
```

Both the constant and the two regimes are sharp given only `(E,M)`.  A flat
vector with `s` nonzero coordinates, all of magnitude `M`, has

```math
K_H(x)=
\begin{cases}
sM=E/M,&s\le H,\\
M\sqrt{Hs}=\sqrt{HE},&s\ge H.
\end{cases}
```

Thus no universally stronger function of only `E`, `M`, and `H` can replace
(R53.4).  Improvement requires truncated profile information such as
(R53.2)--(R53.3).

## 3. Exact map to the completion problem

Use precisely the notation of (10.1243)--(10.1244) and (10.1275)--(10.1277):

```math
u=A[T,S]y,\qquad \beta=2\sigma u,
\qquad G=\lVert u\rVert_2^2,
\qquad J=\lVert u\rVert_\infty.
```

The negative threshold is

```math
\boxed{
r=-\frac g{p_2}
=\frac{Q(A[S])-p^{3/2}q-t-(1-p_2)e}{p_2}>0,}
\tag{R53.5}
```

where `e=sigma y^T A[S]y`.  Use the actual Hanson--Wright cutoff from the
Wave 52 proof,

```math
u_H=\frac{H\log12+\log4}{c_{\rm HW}},\qquad
b_H=\lVert A[T]\rVert_F\sqrt{u_H}
    +\lVert A[T]\rVert_{\rm op}u_H,
\qquad R=r+b_H.
```

Homogeneity and (R53.4) give

```math
\frac12\mathcal K_H(\beta)=K_H(u)
\ge\min\{G/J,\sqrt{HG}\}.
\tag{R53.6}
```

Consequently the **constant-one scalar certificate** is

```math
\boxed{
G\ge\frac{R^2}{H},\qquad G\ge JR
\quad\Longrightarrow\quad
r+b_H\le\frac12\mathcal K_H(\beta).}
\tag{R53.7}
```

Equivalently, in the original beta variables `(E_beta,M_beta)`, the two
conditions are

```math
E_\beta\ge\frac{4R^2}{H},\qquad E_\beta\ge2M_\beta R.
```

The simpler pair `G>=R^2/H` and `J<=R/H` also suffices, but is unnecessarily
strong: the sharp consequence of (R53.4) permits `J<=G/R`.  The head/tail
certificate (R53.3) gives the stronger alternative

```math
\lVert u_I\rVert_1+\sqrt{(H-|I|)\lVert u_{I^c}\rVert_2^2}\ge R
```

whenever its displayed tail-box condition holds.

Combining (R53.7) with (10.1276) yields, without any independence assumption,

```math
\Pr_w\{L(w)+Q_T(w)\le-r\}\ge\tfrac12\,12^{-H}.
```

## 4. Target-scale interpretation and exact missing abundance lemma

Put

```math
H=\lceil n^{3/4-c}\rceil,\qquad T_n=n^{3/2-c},
\qquad 0<c<1/4.
```

On a far state `r>=omega_n T_n`, with `omega_n` tending to infinity,
`b_H=O(T_n)=o(r)`.  Hence (R53.7) asks, at its minimum-energy edge, for

```math
G\gtrsim\omega_n^2 n^{9/4-c},
\qquad J\lesssim R/H\asymp\omega_n n^{3/4}.
\tag{R53.8}
```

This identifies the role of delocalization exactly.  The energy scale is the
necessary scale (10.1278), strengthened by the diverging far-margin factor;
at that edge the cross energy must be spread over about `H` effective
coordinates rather than concentrated in a much smaller saturated head.

Using only the trivial `J<=m=O(n)` in (R53.7) instead demands

```math
G\gtrsim mR\asymp\omega_n n^{5/2-c}.
```

The exact-minimizer operator estimate gives only

```math
G\le\lVert A\rVert_{\rm op}^2\lVert y\rVert_2^2
\le2qm=O(n^{5/2}),
```

so this worst-coordinate implementation is feasible only for a sufficiently
slow `omega_n=o(n^c)` and operates near the global ceiling.  The truncated
certificate is therefore materially preferable.

Define the scalar-certified population

```math
\mathcal D_n(\omega)=\left\{(S,\sigma,y):
\begin{array}{l}
g<0,\quad r\ge\omega_nT_n,\\
G\ge(r+b_H)^2/H,\\
G\ge J(r+b_H)
\end{array}\right\}.
```

The exact remaining population theorem is

```math
\boxed{
\text{For some }\omega_n\uparrow\infty\text{ and }C_0<\infty,
\qquad \nu_m(\mathcal D_n(\omega))\ge e^{-C_0H}}
\tag{R53.9}
```

uniformly over the required exact minimizers, compact density window, and
tolerances.  Equation (R53.7) makes `D_n(omega)` a subset of the `G_n(omega)`
in (10.1277), so (R53.9) proves that open lemma and hence convergence.  One
may enlarge `D_n` by replacing its last two inequalities with (R53.3).

This quantifies the missing assertion: it is not another deterministic tail
bound, but an `exp{-O(H)}` lower bound on the *joint* population of far local
margins, high cross energy, and a nonsaturated cross-field profile.

## 5. Why known exact-minimizer constraints do not supply (R53.9)

The currently available pointwise constraints include

```math
J\le m,\qquad
2\lVert u\rVert_1=\lVert\beta\rVert_1\le q,
\qquad
e^2+4G+2k(k-1)\le q^2.
\tag{R53.10}
```

These are upper bounds.  They do not give `G>=R^2/H`, and even that energy
bound would not control how `G` is distributed among coordinates.  Also,
for fixed `S` and uniform local `y`,

```math
\mathbb E_yG=mk=\Theta(n^2),
```

which is below (R53.8).  A saved lower tail at the larger scale must use an
exact-minimizer-specific correlation with the far-margin event; an average
or a generic upper concentration estimate cannot prove it.

Port bounds can make `J=O(sqrt(n))` for certain maximal-selector child
grounds, but the local-state law in (10.1277) ranges over all `y`.  Applying
that child-ground estimate here would silently change the population and is
not valid.

There is an exact finite pointwise obstruction in the stored `A8` minimizer.
At

```math
n=8,\quad m=4,\quad H=3,\quad
S=\{0,1,2,5\},\quad \sigma=-1,
\quad y=(1,-1,-1,1),
```

one has

```math
q=20,\quad Q(A[S])=12,\quad e=4,\quad
u=(0,0,4,4),
```

and, at `t=0`,

```math
r=\frac{124-70\sqrt2}{3}=8.335016877\ldots.
```

Here `G=32`, `J=4`, and support `(u)=2<=H`, so

```math
\sqrt{HG}=\sqrt{96}>r,
\qquad K_H(u)=\lVert u\rVert_1=G/J=8<r.
\tag{R53.11}
```

Thus the necessary energy inequality passes but coordinate concentration
defeats the K-profile, even before adding the positive `b_H`.  Exhaustive
enumeration finds 750 negative-margin states in this `A8,m=4` instance; 252
pass the energy-only test and 12 of those fail the exact K-profile.  This is
an exact finite algebraic obstruction, not an asymptotic or saved-mass
counterexample to (10.1277).

For a scalable illustration, fix `0<eta<min(c,1/4)`, take even `m~pn`, and
put

```math
d=\lfloor n^{1/2-c+\eta}\rfloor,
\qquad u=(\underbrace{m,\ldots,m}_{d\text{ entries}},0,\ldots,0).
```

This cross field is rigorously realizable by a rectangular sign block
relative to the chosen `y`: use `d` rows aligned with `y`, and balanced rows
elsewhere.  Since `d=o(H)`,

```math
K_H(u)=dm,\qquad G=dm^2,\qquad J=m.
```

For the formal far margin `r=(1+delta)dm`, with fixed `delta>0`,

```math
\frac r{T_n}\asymp n^\eta\to\infty,
\qquad
\frac{GH}{r^2}\asymp\frac Hd\to\infty,
\qquad K_H(u)<r.
\tag{R53.12}
```

Moreover `2||u||_1=o(n^{3/2})` and `G=o(n^3)`, so this profile respects the
scales of the `l_1` cap and Parseval cap in (R53.10).  One can formally take
`q=Theta(n^{3/2})`, `e=0`, and `h=p_2r` to satisfy the scalar margin relation
while keeping `Q(A[S])=p^{3/2}q+t+h<q` for large `n`.

This is explicitly an **abstract cross-profile/scalar obstruction**, not a
claimed full exact-minimizer construction: the planted cross block and the
formal scalar assignments have not been embedded into an exact minimizing
sign matrix with the asserted local margin.  Its valid scope is to show that
high cross energy and all currently displayed scalar caps cannot by
themselves force the profile.  A theorem correlating exact minimizing
geometry with delocalization could still prove (R53.9).

## 6. Reproducibility and final judgment

`tmp/kprofile_scalar_r53_check.py` and its saved output
`tmp/kprofile_scalar_r53_check.out` verify deterministically:

- (R53.1)--(R53.4) on 5,100 seeded integer vectors and the sharp flat cases;
- the head/tail certificates on every feasible tested prefix;
- the exact stored `A8` witness (R53.11), including the radical value of `r`;
- exhaustive finite-state counts for the pointwise failure; and
- every exponent comparison in the abstract target-scale model.

**Final research judgment:** (R53.7) is a useful, exact, and checkable
high-cross plus delocalization reduction of (10.1277), while (R53.2)--(R53.3)
are the sharper implementation when exceptional coordinates matter.  The
deterministic coefficient problem is now closed at constant one.  The live
mathematical burden is the joint `exp{-O(H)}` exact-minimizer abundance
statement (R53.9); no known constraint proves it, and energy alone is
pointwise falsified as a substitute.
