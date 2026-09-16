# Wave 25 memo: priced exact-ground separation

## Status

The LP geometry, two-ground formula, exact-ground exchange identities,
normal-cone obstruction, falsification criterion, and finite audits below are
**proved**.  The desired active-density lower bound remains **Open**.  Exact
`A_6,A_8,A_9` calculations are checked by
`tmp/priced_ground_separation_r25_check.py`.  Nothing here proves
convergence.

## 1. Exact price geometry, including the negative bracket

Fix an exact order-`n` minimizer `A`, an active target order `m`, and a
tolerance `t`.  For every exact positive parent ground `g`, put

```math
u_g=U_m\{S:Q(A[S])-c_A(S,g)\le t\},
\qquad c_g=R_2(g).
```

For budget `C`, retain the dummy action `(u,c)=(0,0)`.  The verified LP and
dual are

```math
Z_*(C)=
\max_\nu\left\{\sum_g\nu_gu_g:
\sum_g\nu_gu_g(c_g-C)\le0\right\}
=\min_{\theta\ge0}\max_g
u_g[1+\theta(C-c_g)].
```

Because of the dummy, the last maximum is equivalently

```math
\boxed{
Z_*(C)=\min_{\theta\ge0}\max_{g\in\Gamma(A)}
\bigl[u_g(1+\theta(C-c_g))\bigr]_+.}
\tag{P.1}
```

The positive part is not an extra relaxation: maximizing over the dummy and
the real columns is exactly
`max(0,max_g u_g[1+theta(C-c_g)])`.  A real action with `u_g=0` is identical
to the dummy for both objective and constraint, regardless of its listed
cost, and cannot subsidize a positive-coverage action.

Thus for `theta>0` every ground with

```math
c_g\ge C+1/\theta
```

has a nonpositive bracket and contributes exactly zero.  It must not be used
as a positive learner.

There is a closed primal formula.  If `c_a<C<c_b`, define

```math
\boxed{
Z_{a,b}(C)=
\left[
\frac{c_b-C}{c_b-c_a}\frac1{u_a}
+\frac{C-c_a}{c_b-c_a}\frac1{u_b}
\right]^{-1},}
\tag{P.2}
```

with value zero if either coverage is zero.  Then

```math
\boxed{
Z_*(C)=\max\left\{
\max_{g:c_g\le C}u_g,
\max_{a,b:c_a<C<c_b}Z_{a,b}(C)
\right\}.}
\tag{P.3}
```

To prove this, first dispose of the case `Z=0`.  If `Z>0`, no optimizer uses
the dummy or a zero-coverage real action: deleting such mass and rescaling
the positive actions preserves the sign of the constraint and increases
the objective.  On the positive-coverage support, use the captured law
`mu_g=nu_g u_g/Z`.  Then

```math
\frac1Z=\sum_g\frac{\mu_g}{u_g},
\qquad \sum_g\mu_gc_g\le C.
```

Conversely, every such `mu` reconstructs
`Z=(sum_g mu_g/u_g)^{-1}` and `nu_g=Zmu_g/u_g`.  Hence one minimizes
`sum_g mu_g/u_g` under one cost moment.  An extreme point has at most two
entries.  Two entries both below budget cannot beat the better singleton,
two above budget are infeasible, and a binding crossing pair has exactly the
weights in (P.2).

This gives the sharp structural dichotomy:

- either one under-budget exact ground already has the required coverage;
- or one ground below `C` and one above `C` close by the harmonic interpolation
  (P.2).

Write `Gamma_+={g:u_g>0}`.  If `Gamma_+` is empty, then `Z_*(C)=0`.
Otherwise the exact lower endpoint behavior is

```math
\boxed{
\begin{aligned}
\min_{g\in\Gamma_+}c_g>C&\Longrightarrow Z_*(C)=0,\\
\min_{g\in\Gamma_+}c_g=C&\Longrightarrow
Z_*(C)=\max_{g:c_g=C}u_g.
\end{aligned}}
\tag{P.4}
```

The second line holds because no positive captured mass above `C` can be
balanced when no positive-coverage ground lies below `C`.  At a crossing
endpoint `C=c_a` or `C=c_b`, (P.2) reduces to the corresponding singleton,
which is why (P.3) uses strict crossing inequalities.

## 2. Exact separation/falsification criterion

For a proposed overlap `z>0`, strong duality gives the exact alternative

```math
\boxed{
Z_*(C)<z
\quad\Longleftrightarrow\quad
\exists\theta\ge0\ \ \forall g\in\Gamma(A):
\bigl[u_g(1+\theta(C-c_g))\bigr]_+<z.}
\tag{P.5}
```

For `theta>0`, this says explicitly:

```math
\begin{cases}
u_g<\dfrac{z}{1+\theta(C-c_g)},
&c_g<C+1/\theta,\\[2mm]
\text{no condition},&c_g\ge C+1/\theta.
\end{cases}
\tag{P.6}
```

Consequently the Wave 25 target is falsified for proposed scales if, along
infinitely many active-ratio pairs and for every eligible exact minimizer,
one can exhibit prices `theta_n` for which (P.5) holds with

```math
z_n=\exp\{-\omega(n^{3/4-c})\}.
```

Conversely, proving the target means defeating every such price.  Equations
(P.3) and (P.5) are equivalent primal and separating formulations; neither
silently assumes that both grounds are individually under budget.

## 3. What exact-ground exchange really gives

Use unordered edges `e={i,j}` throughout this section; the original
quadratic form counts both orientations, so an oriented positive ground has
unordered energy `M(A)=q_n/2`.  Gauge one ground `g=(sigma,x)` to the
all-one state and put

```math
w_{ij}=\sigma a_{ij}x_ix_j,
\qquad r_i(g)=\sum_{j\ne i}w_{ij}.
```

Let `h=(tau,y)` be another exact positive ground, set

```math
\varepsilon=\sigma\tau\in\{+1,-1\},
\qquad z_i=x_iy_i,
\qquad H=\{i:z_i=-1\},
```

and let

```math
b_g(H)=\sum_{i\in H,j\notin H}w_{ij}.
```

Comparing the two exact endpoint energies gives the extremal-shore identity

```math
\boxed{
b_g(H)=
\begin{cases}
0,&\varepsilon=+1,\\
q_n/2,&\varepsilon=-1.
\end{cases}}
\tag{P.7}
```

For the cross degree

```math
a_i(H)=\sum_{j:z_j\ne z_i}w_{ij},
```

the new row fields and row-square cost are exactly

```math
\boxed{
r_i(h)=\varepsilon[r_i(g)-2a_i(H)],
\qquad
R_2(h)-R_2(g)
=4\sum_i[a_i(H)^2-r_i(g)a_i(H)].}
\tag{P.8}
```

Both endpoint row vectors are nonnegative.  In particular, for equal
orientation `a_i<=r_i/2`, while for opposite orientation
`a_i>=r_i/2`; the total cross degrees are respectively zero and `q_n`.

The selector loss exchange is equally exact but has no useful sign.  For an
`m`-set `S`, let

```math
E_g(S)=\sum_{i<j\in S}w_{ij},
\qquad
b_g(H;S)=
\sum_{i\in H\cap S,j\in S\setminus H}w_{ij}.
```

Then

```math
\boxed{
\ell(S,h)-\ell(S,g)
=2(1-\varepsilon)E_g(S)+4\varepsilon b_g(H;S).}
\tag{P.9}
```

For equal orientation this is `4b_g(H;S)`.  Although the full cut in (P.7)
has weight zero, its restriction to `S` can have either sign.  Likewise the
quadratic expression in (P.8) can have either sign.  Therefore an exchange
which lowers `R_2` need not improve coverage, and an exchange which improves
coverage need not lower `R_2`.  Equations (P.7)--(P.9) are the complete
two-ground switching calculation; an additional minimizer-specific
correlation theorem is required to turn it into the priced lower bound.

There is a useful algebraic restatement of the cost:

```math
\boxed{R_2(g)=\|Ax\|_2^2=x^TA^2x,
\qquad |x^TAx|=q_n.}
\tag{P.10}
```

Thus row-square selection is optimization of the second quadratic form
`A^2` over Boolean maximizers of the first.  Cauchy--Schwarz gives only the
lower floor `R_2(g)>=q_n^2/n`; the known universal upper ceiling remains
`R_2(g)<=2(n-1)q_n`.  Neither reaches the required existence bound
`O(n^{9/4-c})`.

## 4. Why a first-order normal-cone proof cannot close

Let `M(A)=q_n/2` in unordered coordinates and write the oriented edge word
of `g=(sigma,x)` as `d^g_{ij}=sigma x_ix_j`.  The subdifferential of the
cut norm at `A` is the convex hull of these exact-ground edge words.  For
every probability law `nu` on those grounds and the corresponding
subgradient `G=E_nu d^g`,

```math
\boxed{
\sum_e a_eG_e
=\mathbb E_{g\sim\nu}\sum_ea_ed^g_e
=M(A)>0.}
\tag{P.11}
```

If the signing vertex `A` were certified as a minimizer of this convex norm
over the continuous edge box `[-1,1]^E` by that normal cone, first-order
optimality would require

```math
a_eG_e\le0\qquad(e\in E).
```

Indeed the box normal has `a_ev_e>=0`, while stationarity `G+v=0` requires
`a_eG_e<=0`.  This contradicts (P.11), whose factor is exactly `q_n/2`
because each unordered edge occurs once.  Hence no **first-order
continuous-box normal-cone certificate** supported on the exact-ground face
can certify signing minimality or force the desired priced column through
that relaxation.  A discrete edge-flip certificate must bring in changed
active maximizers on the edges where the old exact face has positive
alignment; those witnesses may be positive-deficit states at `A`.  This is
precisely the witness gap already visible in (10.733) and (10.768), now as a
formal normal-cone obstruction.

This is not a no-go theorem for discrete separation generally, nor for a
genuinely discrete second-order exchange theorem.  It rules out importing a
first-order continuous normal-cone certificate without the non-ground
witnesses which that certificate cannot contain.

## 5. Exact finite price audit

For one deletion and zero tolerance, the exact types are

| minimizer | `(u_g,c_g)` | multiplicity |
|:---:|:---:|---:|
| `A_6` | `(5/6,30)` | `12` |
| `A_8` | `(3/8,64)` | `8` |
| `A_9` | `(1/9,80)` | `2` |
| `A_9` | `(2/9,96)` | `12` |
| `A_9` | `(1/3,112)` | `4` |
| `A_9` | `(2/9,112)` | `3` |
| `A_9` | `(1/3,128)` | `4` |

For `A_6`, if `C<30` the price
`theta>=1/(30-C)` makes the only bracket nonpositive and `Z_*=0`; for
`C>=30`, `theta=0` is optimal and `Z_*=5/6`.  The analogous statement for
`A_8` has `(64,3/8)`.

For `A_9`, the lower captured hull has vertices

```math
(c,1/u)=(80,9),(96,9/2),(112,3).
```

At `C=88`, the optimal price is `theta=1/24`, where the `(80,1/9)` and
`(96,2/9)` columns both equal `4/27`.  At `C=104`, the optimal price is
`theta=1/40`, where the `(96,2/9)` and `(112,1/3)` columns both equal
`4/15`.  Every dominated or nonactive column lies below these maxima.  This
checks both the crossing-pair formula and the sign of the price bracket.

The checker additionally verifies (P.7)--(P.10) for every ordered pair of
exact grounds and every one-deletion selector in `A_6,A_8,A_9`.  These
finite minimizers obey the exact exchange identities, but they do not reveal
a monotone cost--coverage direction.  They are diagnostics, not asymptotic
evidence.

## 6. Surviving exact lemma and conclusion

The desired statement remains the following **Open target**: for some active
fixed-ratio window, target-specific exact minimizers, and fixed
`0<c<1/4`, every `theta>=0` admits an exact parent ground satisfying

```math
\boxed{
\bigl[u_g(1+\theta(C-R_2(g)))\bigr]_+
\ge\exp\{-O(n^{3/4-c})\},
\quad
C=O(n^{9/4-c}),
\quad t=O(n^{3/2-c}).}
\tag{P.12}
```

The strongest proved reduction is the exact dichotomy (P.3): it is enough to
find either one adequate under-budget ground or one crossing pair with the
harmonic coverage in (P.2).  Ground switching supplies the exchange formulas
(P.7)--(P.9), but not the correlation needed to choose that pair.  The exact
separating falsifier is (P.5)--(P.6).  First-order normal-cone methods cannot
bridge the gap by (P.11); a positive continuation must exploit genuinely
discrete exact-minimizer structure linking restricted zero/max shores to
row-square improvement.
