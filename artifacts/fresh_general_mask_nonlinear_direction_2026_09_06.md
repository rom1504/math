# Every high-value marked mask has a finite nonlinear edge direction

Date: 2026-09-06. This is a proved Gaussian response inequality. It
strengthens the finite-degree selection step without asserting the local
slack needed to improve every mask's actual Boolean certificate.

## 1. An inequality for arbitrary masks, not just central ones

Let U be the marked-tree Gaussian creation isometry, let `G0=U1`, and
let H be any jointly even measurable mask with `0<=H<=1`. Write

```
p=E H^2, mu=E H, W=U H,
F=sign(W)(1-H), J=E F W, a=E G0 F.
```

Assume `J>0`, so `p>0`. Isometry and boundedness give

```
Var(W)=p, Cov(G0,W)=mu,
p<=mu<=sqrt(p), E F^2=1-2mu+p<=1-p,
J^2<=p(1-p), 0<J<=1/2.
```

Project both G0 and F orthogonally off W in L2. Cauchy--Schwarz yields

```
a >= (mu/p)J
      - sqrt[(1-mu^2/p)(1-2mu+p-J^2/p)]
  >= J - sqrt[(1-p)(1-p-J^2/p)].                         (1)
```

All square roots have nonnegative arguments by the preceding covariance
inequalities. No Gaussian regression assumption about H itself occurs.

When `J>2/(3sqrt(3))`, rationalizing (1) gives

```
a >= [J^2-p(1-p)^2]
     / [p(J+sqrt((1-p)(1-p-J^2/p)))]
  >= (4/3)(J^2-4/27).                                  (2)
```

Indeed `p(1-p)^2<=4/27`, while the positive denominator is at most
`pJ+p(1-p)<=J+1/4<=3/4`. In particular,

```
J>=43/100  =>  a>=9923/202500>49/1000.                    (3)
```

This replaces the weaker central-mask-only estimate `a>=1/200` by a
uniform bound for every even mask attaining this certificate threshold.
The proof uses only the creation isometry and the two-dimensional L2
Gram matrix. It does not use scalar fixed-point coordinates, a particular
anchor family, or an assumption that H is a function of W.

## 2. A degree-417 witness requires no regularity of the mask

For any standard Gaussian G and any bounded random variable F with
`|F|<=1`, put `a=E GF>=a0>0`. Let D be an odd integer at least three,
let `c_D=E G^(D+1)=D!!`, and set

```
p_D(G)=G^D-c_D G,
||p_D||_2^2=(2D-1)!!-(D!!)^2.
```

This polynomial has only Hermite degrees 3,5,...,D. Moreover,

```
|E F p_D(G)| >= a0 c_D-E|G|^D,
E|G|^D <= sqrt[E G^(D-1) E G^(D+1)] = c_D/sqrt(D).
```

Consequently, with normalized probabilists' Hermite polynomials h_r,

```
sum_(3<=r<=D, r odd) |E F h_r(G)|^2
 >= [(a0-D^(-1/2))_+ D!!]^2
    / [(2D-1)!!-(D!!)^2].                               (4)
```

For the masks in (3), take `G=G0`, `a0=49/1000`, and `D=417`.
The right side is strictly positive because `417*49^2=1001217>10^6`.
Thus all such masks have a quantitatively nonzero nonlinear edge
projection in the same finite set of odd degrees through 417.

One deliberately loose numerical form is also certified:

```
sqrt(sum_(3<=r<=417, r odd) |E F h_r(G0)|^2) > 7*10^(-68),
max_(3<=r<=417, r odd) |E F h_r(G0)| > 4*10^(-69).
```

Indeed `1/sqrt(417)<48971/10^6`, leaving margin `29/10^6` below
`49/1000`. The elementary recurrence
`(2D-1)!!/(D!!)^2<=2^(D-1)` bounds the residual denominator in (4),
so the projection norm exceeds `29/(10^6*2^208)`. There are 208 selected
degrees and `sqrt(208)<15`. All rational comparisons are replayed by
`computations/fresh_general_mask_degree_certificate.py`. These are not
claimed sharp constants or an additional evaluated original lower bound.

This is much stronger than a finite-degree argument requiring a uniform
regularity or threshold-noise bound. It follows from a separating odd
monomial and Cauchy--Schwarz, not from a compactness assumption about
the unknown mask. The independent monomial audit also gives general
degree and mass bounds in
`fresh_finite_hermite_monomial_audit_2026_09_06.md`.

The separating polynomial p_D is not the perturbation itself. For the
proved joint weighted transport identity, form the coefficient-aligned
response `h=sum_r f_r h_r/s`, with `f_r=E F h_r(G0)` and
`s^2=sum_r f_r^2`. This preserves the positive Schur-mixture weights
`f_r^2/s^2`; replacing it by an arbitrary polynomial would not.

## 3. What remains missing beyond the scalar family

For each fixed finite approximation, the all-odd weighted identity and
Schur normalization turn (4) into a nonzero jointly evaluated response
direction. They do not automatically give a feasible improvement of
the original mean vector: the latter requires inexpensive local slack
near `W=0`, where the old odd response can be softened at quadratic cost.

The scalar central-mask family has that slack uniformly by its proved
covariance nondegeneracy and Gaussian density bounds. Therefore (3)--(4)
can replace its compact finite-degree selection without altering any
other step of the already proved uniform escape theorem.

For arbitrary even masks, neither positive total slack nor (3) locates
that slack near `W=0`. Uniform local slack at the threshold `J>=.43`
has not been proved. Nor has this note shown a strict improvement over
the supremum of all arbitrary-mask certificates. No additional numerical
lower bound, matching upper construction, or convergence theorem follows
from this note alone.

Independent verification: algebra checked (1)--(3), all normalization and
rationalization constants; literature independently derived (4) for general
bounded F, without any measurability or parity shortcut. The director
reconstructed both arguments.

## 4. A genuine creation-space slack counterexample at lower values

The high-value qualification cannot simply be removed. Define the smooth
even function

```
u(z)=0                         if |z|<=1,
u(z)=exp(-1/(z^2-1))           if |z|>1,
f(z)=1-u(z)/10.
```

Then `.9<=f<=1`. With `t=1/(z^2-1)`, its nonzero derivative obeys
`|u'|=2 t^(3/2)sqrt(t+1)e^(-t)<4`. Thus, putting `p=E f(G)^2`,
we have `p>=.81` and `E f'(G)^2<=.16<p`. The proved scalar creation
contraction gives a unit first-chaos V with `V=U[f(V)/sqrt(p)]`.
For `H=f(V)` and `W=U H=sqrt(p)V`,

```
H=1 whenever |W|<=sqrt(p),
J=(sqrt(p)/10) E|G|u(G) in (0,.08),
a=(E H/p)J>0.
```

This is an actual mask in the marked-tree creation space, not a relaxed
Gaussian covariance model. It has positive total slack and a positive
edge coefficient, but zero slack on a whole neighborhood of `W=0`.
Its low J does not settle the remaining high-value question. The algebra
agent proposed and checked this example; the director reconstructed its
smoothness, derivative bound, and fixed-point normalization.
