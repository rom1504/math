# Wave 24 memo: joint regularity--coverage duals

## Status

The finite minimax programs, switching identities, soft cost-sensitive LP,
LP dual, aggregate converse, and finite certificates below are **proved**.
The exact `A_6,A_8,A_9` audit is checked by
`tmp/joint_regular_minimax_r24_check.py`.  No lower bound at an active
fixed density is obtained: the switching identities rewrite the joint event
exactly but do not force its lower tail.  Nothing here proves convergence.

## 1. Put the row cap inside the capped game

Let `X=Omega_{n,m}`, let `Gamma(A)` be the exact positive parent grounds,
and use ordered loss

```math
\ell(S,g)=Q(A[S])-c_A(S,g).
```

For fixed `B,t`, define one joint incidence matrix

```math
a_{Sg}^{B,t}
=1\{R_\infty(g)\le B,\ \ell(S,g)\le t\}.
```

Writing `N=|X|` and using selector probabilities `pi_S` rather than their
densities, the exact capped learner is the LP

```math
\boxed{
\begin{aligned}
\alpha_{B,t,\varepsilon}
=\min_{\pi,\alpha}\quad&\alpha\\
\text{s.t.}\quad
&\sum_S\pi_Sa_{Sg}^{B,t}\le\alpha &&(g\in\Gamma(A)),\\
&0\le\pi_S\le\frac1{N\varepsilon},
\qquad \sum_S\pi_S=1.
\end{aligned}}
\tag{J.1}
```

Thus irregular grounds remain in the action set but have an identically zero
column; regularity and coverage are not proved in separate stages.  Finite
minimax gives the exact common-prior dual.  For `nu in Delta(Gamma(A))` put

```math
p_\nu(S)=\sum_g\nu_ga_{Sg}^{B,t}.
```

Then

```math
\boxed{
\alpha_{B,t,\varepsilon}
=\max_{\nu,z,(y_S)}
\left\{z-\frac1{N\varepsilon}\sum_Sy_S\right\},
\quad
y_S\ge z-p_\nu(S),\quad y_S\ge0.}
\tag{J.2}
```

The objective is the mean of the lowest `epsilon`-fraction of the joint
coverage probabilities.  Equations (J.1)--(J.2) are the exact hard-cap
minimax/LP dual requested in (10.741), with the cap inside the game.

## 2. Exact switching form of every column

Gauge a parent ground `g` to the all-one state and write its switched signing
as `W_g=(w_{ij})`.  For a deleted shore `T=V\S`, put

```math
\begin{aligned}
r_i&=\sum_{j\ne i}w_{ij},
&R_T&=\sum_{i\in T}r_i,\\
h_T&=1_T^T W_g[T]1_T,
&b_T&=1_T^T W_g[T,S]1_S,\\
e_T&=1_S^T W_g[S]1_S,
&d_T&=q_n-Q(A[S]),\\
\gamma_T&=Q(A[S])-e_T.
\end{aligned}
```

Here `gamma_T=ell(S,g)`.  Expansion gives all of the exact
field-proportional deletion identities

```math
\boxed{
d_T+\gamma_T=R_T+b_T,
\qquad
2R_T=d_T+h_T+\gamma_T.}
\tag{J.3}
```

Ground switching by an arbitrary shore gives

```math
\boxed{
r_i\ge0,\qquad \sum_i r_i=q_n,
\qquad 0\le b_T\le q_n/2,
\qquad 0\le d_T\le q_n-q_m,
\qquad \gamma_T\ge0.}
\tag{J.4}
```

The upper cut bound uses that the shore-flipped energy
`q_n-4b_T` lies in `[-q_n,q_n]`.  Therefore the joint incidence is exactly

```math
\boxed{
a_{Sg}^{B,t}
=1\left\{
\max_i r_i(g)\le B,
\ R_T(g)+b_T(g)-d_T\le t
\right\}.}
\tag{J.5}
```

For one deletion this specializes to

```math
\boxed{\gamma_{\{i\}}(g)=2r_i(g)-d_i.}
\tag{J.6}
```

These formulas use the row cap and coverage jointly.  They also expose the
remaining obstruction.  Uniformly over `m`-sets,

```math
\mathbb E e_T=p_2q_n,
\qquad
\mathbb E\gamma_T
=\mathbb E Q(A[S])-p_2q_n.
\tag{J.7}
```

The second expression is exactly the unknown restriction excess.  Cut
stability (J.4) bounds `b_T` but gives no lower-tail mass for (J.5).  The
simultaneous edge-flip cover (10.733) and switched-block replacement
(10.702) do not repair this projection: their existential maximizers may be
positive-deficit states, whereas every column of (J.1) is an exact parent
ground.  `A_9` realizes such witness migration.  A theorem transferring
those witnesses to exact-ground columns with controlled `R_2` would be a new
input, not a consequence of the displayed switching identities.

## 3. The soft row-square theorem collapses the adversarial game

The domain-free row-square bound and global product conditioning in
`tmp/soft_row_entropy_r24.md` strictly weaken (J.1).  The calculation is
sound: if `nu` is a parent-ground prior and

```math
G=\{(S,g):\ell(S,g)\le t\},
\qquad Z=(U_m\otimes\nu)(G),
```

then conditioning `U_m x nu` on `G` has total relative entropy `-log Z`.
The chain rule gives

```math
I(S;D)+D(\pi\Vert U_m)\le-\log Z,
```

and its row-square cost is exactly `E[R_2(D)|G]`.  Allowing the selector
marginal to move is legitimate because its KL cost is already present in the
uniform-reference inequality.  There is no missing worst-selector term.

This joint target has an exact cost-sensitive LP.  For each parent ground put

```math
u_g=U_m\{S:\ell(S,g)\le t\},
\qquad c_g=R_2(g).
```

For a desired captured row-square budget `C`, adjoin a dummy action with
`u=0` and define

```math
\boxed{
\begin{aligned}
Z_*(C)=\max_{\nu\in\Delta(\Gamma\cup\{0\})}\quad
&\sum_g\nu_gu_g\\
\text{s.t.}\quad
&\sum_g\nu_gu_g(c_g-C)\le0.
\end{aligned}}
\tag{J.8}
```

The constraint is precisely
`E[R_2(D)|G]<=C` when the objective is positive.  The dummy makes the LP
feasible and gives value zero when no positive captured mass can meet the
budget.  Strong LP duality gives the exact priced weak learner

```math
\boxed{
Z_*(C)
=\min_{\theta\ge0}\max_{g\in\Gamma\cup\{0\}}
u_g[1+\theta(C-c_g)].}
\tag{J.9}
```

Thus `Z_*(C)>=z` is equivalent to: **for every row-square price
`theta>=0`, some parent ground has priced captured mass at least `z`**.
The price multiplies `u_g(c_g-C)`, so a spiky ground is charged only through
the selector mass it actually captures.

There is an equally useful captured-output form.  Any optimizer with
`Z_*(C)>0` gives zero mass to the dummy and to every zero-coverage ground:
removing such mass and renormalizing preserves the cost constraint and
strictly increases `Z`.  For such an optimizer put

```math
\mu_g=\frac{\nu_gu_g}{Z}.
```

Then

```math
\boxed{
\frac1Z=\sum_{g:u_g>0}\frac{\mu_g}{u_g},
\qquad
\sum_g\mu_gc_g\le C.}
\tag{J.10}
```

Conversely every such `mu` reconstructs
`nu_g=Z mu_g/u_g`.  Hence maximizing `Z` is the linear program of minimizing
`sum mu_g/u_g` under one moment constraint.  An optimizer uses at most two
positive-coverage grounds.  The soft theorem therefore needs no iterative
boosting: its exact joint prior can always be chosen with two captured
ground types (plus the harmless dummy at value zero).

For the power-saving route it is enough that, with

```math
t=O(n^{3/2-c}),
\qquad C=O(n^{9/4-c}),
```

one has

```math
\boxed{Z_*(C)\ge\exp\{-O(n^{3/4-c})\}.}
\tag{J.11}
```

This is strictly weaker than the capped hard-row learner: it permits a rare
spiky ground when its captured-output weight is small.

## 4. Matching converse and scoped obstruction

The same row-square mgf gives a converse for the joint LP.  For a common
`lambda` in the operator domain, Chernoff gives every parent ground

```math
-\log u_g
\ge
\lambda[q_m-t-(p_2+\epsilon_{n,m})q_n]
-\chi_{n,m}
-\lambda^2p^2c_g
-O(\lambda^2n^2).
```

Average this under the captured law `mu` in (J.10) and use
`-log Z=log E_mu(1/u_g)>=E_mu[-log u_g]`.  If `E_mu c_g<=C`,

```math
\boxed{
-\log Z
\ge
\lambda[q_m-t-(p_2+\epsilon_{n,m})q_n]
-\chi_{n,m}-\lambda^2p^2C-O(\lambda^2n^2).}
\tag{J.12}
```

If an active-ratio sequence has
`q_m-p_2q_n>=g n^{3/2}`, take
`lambda=eta n^{-3/4}`.  With
`C=O(n^{9/4-c})` and small fixed `eta`, (J.12) yields

```math
\boxed{-\log Z_*(C)=\Omega(n^{3/4}).}
\tag{J.13}
```

This rules out (J.11), whose allowed exponent is `O(n^{3/4-c})`.
Present global constants do not force such a gap for `rho>=1/2`, so this is
a precise conditional falsifier, not an active-density counterexample.

Equations (J.3)--(J.5) do not prove the opposite lower bound.  They constrain
the mean and the shore-cut range but do not force the required below-mean
tail.  A scalar relaxation makes the scope explicit: take regular fields
`r_i=q/n` and the constant complete weight
`w_{ij}=q/[n(n-1)]`.  It obeys all ground switching and deletion identities,
with

```math
e_S=p_2q,
\qquad
b_T=\frac{q|T||S|}{n(n-1)}.
```

Assign an admissible child norm `Q_S=e_S+gamma` with
`t<gamma=o(q)`.  Then every identity (J.3)--(J.4), including
`0<=b_T<=q/2`, holds while `u_g=0` although the ground is perfectly regular.
This is a relaxation witness, not a `+-1` signing: it proves only that the
listed scalar switching identities cannot establish (J.11).  Signing
realizability plus a new exact-ground overlap theorem remains available.

## 5. Exact finite joint audit

For one deletion, (J.6) makes the audit especially transparent.  At
ordered tolerance zero, the parent-ground `(u_g,c_g)` types are:

| minimizer | captured fraction `u_g` | row square `c_g` | multiplicity |
|:---:|---:|---:|---:|
| `A_6` | `5/6` | `30` | `12` |
| `A_8` | `3/8` | `64` | `8` |
| `A_9` | `1/9` | `80` | `2` |
| `A_9` | `2/9` | `96` | `12` |
| `A_9` | `1/3` | `112` | `4` |
| `A_9` | `2/9` | `112` | `3` |
| `A_9` | `1/3` | `128` | `4` |

The last two `A_9` types are dominated.  Therefore

```math
Z_*^{A_6}(C)=
\begin{cases}0&C<30,\\5/6&C\ge30,\end{cases}
\qquad
Z_*^{A_8}(C)=
\begin{cases}0&C<64,\\3/8&C\ge64.\end{cases}
```

For `A_9`, let `L(C)=1/Z_*(C)`.  Exact convexification of (J.10) gives

```math
\boxed{
L(C)=
\begin{cases}
+\infty,&C<80,\\
9-\frac9{32}(C-80),&80\le C\le96,\\
\frac92-\frac3{32}(C-96),&96\le C\le112,\\
3,&C\ge112.
\end{cases}}
\tag{J.14}
```

For example, `Z_*(80)=1/9`, `Z_*(88)=4/27`,
`Z_*(96)=2/9`, `Z_*(104)=4/15`, and `Z_*(112)=1/3`.
This is a genuinely joint audit: captured mass and row square are optimized
together.  In particular, the cap-four `A_9` grounds have zero hard capped
value at exact tolerance for a large selector tail, yet the soft joint LP
uses their `1/9` captured mass at cost `80` and pays only `log 9`.  The soft
criterion bypasses the finite separation wall exactly as intended.  These
finite constants do not imply an asymptotic active-ratio bound.

## Updated route conclusion

The hard regular capped game has the exact dual (J.1)--(J.2), and switching
turns each column into the replenishment event (J.5).  No available
switching-minimality identity forces that event on a positive lower tail.

The domain-free soft row-square theorem materially improves the target.  Its
joint overlap problem is exactly the one-price LP (J.8), with dual (J.9), and
an optimizer needs at most two captured ground types.  This removes both the
hard row cap and the worst-selector tail.  The remaining sufficient lemma is
(J.11): find one or two parent-ground types with enough uniform captured mass
and subcritical captured row-square average.  The matching conditional
falsifier is (J.13).  Proving the positive direction still requires new
minimizer-specific exact-ground overlap; the switching identities alone do
not supply it.
