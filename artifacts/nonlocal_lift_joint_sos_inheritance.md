# Seed-dependent nonlocal lifts and joint-SOS inheritance

Status: candidate theorem stated before finite computation; exact inheritance
theorem proved below.  The conclusion is a route obstruction, not a bound or
recurrence improvement.

Use the one-copy normalization

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad \operatorname{cap}(A)=\max_x|H_A(x)|.       \tag{1}
```

## 1. Minimal polynomial-state theorem

Fix a perfect-square fibre size `s` and an SOS level `d`.  For an order-`n` seed
signing `A`, let `L_s(A)` be the class of all order-`ns` signings `B`, in
`n` fibres of size `s`, satisfying only the exact quotient constraints

```math
\mathbf1^{\mathsf T}B_{ij}\mathbf1
=s^{3/2}a_{ij}\qquad(i<j).                         \tag{2}
```

Every cross block may depend jointly on the complete seed.  Thus this class
contains the nonlocal dependent-profile lifts and does not impose scalar
atoms, common block templates, separate channel payments, or independent
edge orientations.

There is no loss in using the all-one quotient channel.  If a construction
instead has Boolean fibre profiles `p_i` with
`p_i^T B_ij p_j=s^(3/2)a_ij`, switching microvertex `(i,alpha)` by
`p_(i,alpha)` converts it to (2), while preserving both `cap` and every SOS
endpoint.  Thus the audit includes an arbitrary seed-dependent choice of
one exact Boolean quotient profile.

For a signing `C`, define its degree-`2d` pseudoexpectation endpoints

```math
\begin{aligned}
U_d(C)&=\sup_{\widetilde{\mathbb E}}
             \widetilde{\mathbb E}H_C,\\
L_d(C)&=\inf_{\widetilde{\mathbb E}}
             \widetilde{\mathbb E}H_C,\\
W_d(C)&={U_d(C)-L_d(C)\over2},\\
R_d(C)&=\max\{U_d(C),-L_d(C)\},                   \tag{3}
\end{aligned}
```

where `E~(1)=1`, `E~(p^2)>=0` for `deg p<=d`, and
`E~((x_i^2-1)q)=0` whenever the degree is at most `2d`.  This is a
semidefinite program of size `n^O(d)`, hence polynomial for fixed `d`, and

```math
\operatorname{cap}(C)\le R_d(C).                  \tag{4}
```

The smallest direct joint-SDP absorption theorem strong enough for a
fixed-factor recurrence is:

> **Candidate JSA(`s,d`).**  Every centered/chiral near-minimizing seed `A`
> in the intended family admits a seed-dependent `B in L_s(A)` such that
>
> ```math
> R_d(B)\le s^{3/2}\operatorname{cap}(A)+O_s(n).
> \tag{5}
> ```

The state in (3) is demonstrably smaller than parent maximization: for fixed
`d` it is one polynomial-size moment matrix, while the parent has
`2^(ns-1)` projective spins.  Equation (5) preserves all channel
cancellations inside one SDP before applying an absolute value.  If true on
a recursively closed centered family, it would give

```math
\operatorname{cap}(B)\le
s^{3/2}\operatorname{cap}(A)+O_s(n),               \tag{6}
```

and therefore a geometrically summable `2/3`-power defect.

The rest of this note audits (5).  No finite computation was performed
before the candidate and its polynomial state were fixed above.

## 2. Exact quotient-inheritance theorem

Let the parent variables be `y_(i,alpha)`, where `i<=n` and `alpha<=s`.
The substitution

```math
\Phi(y_{i,\alpha})=x_i                              \tag{7}
```

preserves degree and takes every parent cube relation
`y_(i,alpha)^2=1` to a child cube relation.  Hence every degree-`2d`
pseudoexpectation `E~_A` pushes forward to a valid parent pseudoexpectation

```math
\widetilde{\mathbb E}_B[p(y)]
=\widetilde{\mathbb E}_A[\Phi p].                  \tag{8}
```

Let

```math
\delta_B=\sum_i\sum_{\alpha<\beta}
 b_{(i,\alpha),(i,\beta)}                          \tag{9}
```

be the total internal-fibre energy of the all-one microspin.  Equations
(2) and (7) give the polynomial identity

```math
\Phi H_B=s^{3/2}H_A+\delta_B.                       \tag{10}
```

Applying (8) to optimizing child pseudoexpectations proves, for every
`B in L_s(A)` and every `d>=1`,

```math
\boxed{\begin{aligned}
U_d(B)&\ge s^{3/2}U_d(A)+\delta_B,\\
L_d(B)&\le s^{3/2}L_d(A)+\delta_B.
\end{aligned}}                                     \tag{11}
```

Put `mu_d(A)=(U_d(A)+L_d(A))/2`.  Taking the larger absolute endpoint in
(11) yields the sharper exact form

```math
\boxed{
R_d(B)\ge s^{3/2}W_d(A)
 +|\delta_B+s^{3/2}\mu_d(A)|.}                     \tag{12}
```

This proof never separates microchannels or applies polarization.  It keeps
the complete parent moment matrix and nevertheless shows that the quotient
pseudo-solution survives every seed-dependent choice of cross blocks.

If a signed permutation sends `A` to `-A`, substitution by that signed
permutation bijects the pseudoexpectation feasible set and gives

```math
L_d(A)=-U_d(A),\qquad \mu_d(A)=0,\qquad W_d(A)=R_d(A).
\tag{13}
```

Consequently every exact compressed lift of a chiral seed obeys

```math
\boxed{R_d(B)\ge s^{3/2}R_d(A)+|\delta_B|.}        \tag{14}
```

Candidate JSA(`s,d`) can therefore hold only if the same fixed-level SOS
relaxation is already project-scale exact on every seed to which it is
applied:

```math
R_d(A)-\operatorname{cap}(A)=O_s(n).               \tag{15}
```

Thus a joint SDP does not use the lift to repair its child integrality gap;
it inherits that gap with the full factor `s^(3/2)`.

### 2.1 Scope beyond SOS

The argument uses only closure under coordinate duplication.  It applies
verbatim to any convex upper relaxation whose feasible generalized moments
on `n` Boolean variables push forward under (7).  In particular it includes
the correlation SDP, every fixed level of the cube SOS/Lasserre hierarchy,
and bounded collections of universal matrix-valued moment constraints.

It does not rule out a relaxation augmented by a new `B`-specific global
valid inequality that cuts off (8), nor a level `d=d(n)` growing far enough
to resolve the child extreme tail.  Either is an additional theorem about
the signing, not a consequence of retaining channels in the standard joint
SDP.

## 3. Degree-two floor

At `d=1`, pseudoexpectations are correlation matrices.  For `t>0`, both

```math
G_+(t)={(A+tI)^2\over n-1+t^2},\qquad
G_-(t)={(A-tI)^2\over n-1+t^2}                     \tag{16}
```

are positive semidefinite with diagonal one.  Since
`tr(A^2)=n(n-1)` and `tr A=0`, they give

```math
\begin{aligned}
U_1(A)&\ge {\operatorname{tr}A^3+2tn(n-1)
                 \over2(n-1+t^2)},\\
L_1(A)&\le {\operatorname{tr}A^3-2tn(n-1)
                 \over2(n-1+t^2)}.
\end{aligned}                                      \tag{17}
```

Therefore

```math
W_1(A)\ge {tn(n-1)\over n-1+t^2}.
\tag{18}
```

Optimizing at `t=sqrt(n-1)` proves the signing-independent floor

```math
\boxed{W_1(A)\ge {n\sqrt{n-1}\over2}.}            \tag{19}
```

Hence the ordinary joint correlation SDP cannot certify an amplification
constant below the spectral `1/2` scale.  On any chiral family with
`cap(A)<(1/2-c)n^(3/2)` for fixed `c>0`, (14) and (19) falsify
JSA(`s,1`) by a leading-order gap.

### 3.1 Exact conference benchmark

If `C` is a symmetric conference signing, so that

```math
C^2=(n-1)I,
```

then (19) is exact.  Indeed, `C <= sqrt(n-1)I` in Loewner order, and every
correlation matrix `G` has `tr G=n`; hence

```math
{1\over2}\operatorname{tr}(CG)
\le {n\sqrt{n-1}\over2}.
```

Equality is attained by the correlation matrix

```math
G=I+{C\over\sqrt{n-1}},
```

whose eigenvalues are zero and two.  Replacing `C` by `-C` gives the lower
endpoint.  Therefore

```math
\boxed{U_1(C)=-L_1(C)={n\sqrt{n-1}\over2}.}       \tag{20}
```

This identifies the obstruction exactly on conference seeds: it is not a
rounding loss and it is not caused by treating the two signs separately.

## 4. Dual pullback: the lift cannot manufacture a certificate

The primal inheritance theorem has an equally useful proof-theoretic form.
Suppose a degree-`2d` SOS certificate proves the parent upper bound `T`,
that is, in the Boolean quotient ring it represents

```math
T-H_B(y)
```

as a sum of squares of polynomials of degree at most `d`.  Apply `Phi` to
the complete identity.  Coordinate duplication preserves degree, sends
squares to squares, and gives from (10) a degree-`2d` child certificate for

```math
T-\delta_B-s^{3/2}H_A(x).                         \tag{21}
```

The same statement holds for the lower endpoint.  Thus any two-sided
parent certificate at threshold `T` pulls back, at the *same degree*, to
child endpoint certificates at thresholds

```math
{T-\delta_B\over s^{3/2}}
\quad\hbox{and}\quad
{-T-\delta_B\over s^{3/2}}.                       \tag{22}
```

In particular, neither a joint SDP dual coupling nor an operator-valued
moment proof can acquire strength merely by moving to the lift, provided
its certificate class is closed under the same-switch substitution (7).
The complete multi-channel parent certificate already contains a child
certificate of equal complexity.

There is also a precise statement for proposed extra constraints.  Let a
linear inequality in parent moments be valid for every Boolean parent
spin.  Its pullback by `Phi` is valid for every Boolean child spin.  If the
parent inequality cuts off the inherited pseudoexpectation (8), its
pullback is therefore a child-valid inequality cutting off the original
child pseudoexpectation.  A `B`-specific global cut can escape (11), but it
does so only by supplying a new integrality-gap separator for `A`; the lift
has not eliminated that obligation.

This is the relevant equivalence boundary.  Fixed-size covariance,
matrix-valued discrepancy, noncommutative moment, and same-switch
Grothendieck states are all caught whenever they:

1. contain every Boolean configuration;
2. are closed under coordinate duplication; and
3. certify the parent through a bounded-degree/complexity cone preserved
   by substitution.

### 4.1 Even all bounded-local Boolean cuts retain the conference gap

One might try to escape the bare correlation SDP by imposing every exact
Boolean marginal constraint on at most `q` coordinates, while retaining
the global PSD correlation matrix.  For fixed `q`, this still does not
remove the conference pseudo-solution.

Let `C^2=(n-1)I`, put `r=sqrt(n-1)`, and suppose

```math
{q(q-1)\over2r}\le1.                              \tag{23}
```

For each `T subset [n]` of size at most `q`, define two signed local
densities on the Boolean cube by

```math
\mu_T^\pm(z)=2^{-|T|}
\left(1\pm {1\over r}
  \sum_{\{i,j\}\subset T}c_{ij}z_i z_j\right).   \tag{24}
```

They are nonnegative by (23), have total mass one, and are consistent under
marginalization because all Fourier terms involving an eliminated variable
vanish.  Their pair correlations are respectively `c_ij/r` and
`-c_ij/r`.  Meanwhile
the global pair-correlation matrices

```math
G^\pm=I\pm C/r                                    \tag{25}
```

are PSD.  Thus the two matching global/local systems are feasible for the
hybrid SDP plus *all* `q`-local cut-polytope/marginal constraints, and they
still attain the two conference endpoints in (20).

This is an exact bounded-local no-go, not a numerical observation.  In
particular, all triangle inequalities leave (20) untouched for conference
orders `n>=10`.  More generally the same pseudo-solution survives local
windows through `q=Theta(n^(1/4))` under the explicit sufficient condition
(23).  Therefore the separator in escape item 1 below must be genuinely
global (or use growing locality beyond this regime); adding a fixed menu of
small-cut facets cannot be the missing ingredient.

### 4.2 Factor audit: signed elliptope versus correlation width

The joint finite-fibre action audit introduces the related state

```math
\Gamma(A)=\max\left|\sum_{i<j}a_{ij}(P_{ij}-Q_{ij})\right|,
\quad P,Q\succeq0,\quad\operatorname{diag}(P+Q)=\mathbf1. \tag{26}
```

It is important that `Gamma` is **not** `W_1` in general.  Put

```math
S_+(A)=U_1(A),\qquad S_-(A)=-L_1(A)=U_1(-A).       \tag{27}
```

Standard SDP duality gives

```math
\begin{aligned}
S_+(A)&=\min\{\mathbf1^Ty^+:
       \operatorname{Diag}(y^+)-A/2\succeq0\},\\
S_-(A)&=\min\{\mathbf1^Ty^-:
       \operatorname{Diag}(y^-)+A/2\succeq0\},    \tag{28}\\
\Gamma(A)&=\min\{\mathbf1^Ty:
       \operatorname{Diag}(y)-A/2\succeq0,
       \operatorname{Diag}(y)+A/2\succeq0\}.
\end{aligned}
```

Thus `Gamma` requires one diagonal budget feasible for both signs, whereas
`W_1=(S_++S_-)/2` averages two independently optimized budgets.  The exact
universal comparison is

```math
\boxed{W_1(A)\le R_1(A)\le\Gamma(A)
       \le2W_1(A)\le2R_1(A).}                     \tag{29}
```

For the first lower inequality use average at most maximum.  For the
second, `(P,Q)=(G,0)` realizes every positive endpoint and `(0,G)` every
negative endpoint in (26).  For the upper inequality, let `y^+,y^-` be
optimal in (28).  Positivity forces both vectors to be componentwise
nonnegative.  Hence `y=y^++y^-` satisfies both joint constraints and has
objective `S_++S_-=2W_1`.

Equality `Gamma=W_1` holds exactly when `S_+=S_-` and a single diagonal
vector is simultaneously optimal for both programs in (28).  Conference
matrices have this property via the constant vector
`y_i=sqrt(n-1)/2`.  Chirality supplies `S_+=S_-` but does not supply a
common optimizer: numerically, the two saved certified order-eight
minimizer orbits give respectively

```math
(W_1,\Gamma)=(11.31370849,14.42220510)
\quad\hbox{and}\quad
(12.00000000,12.00000000).                        \tag{30}
```

Consequently the exact finite-fibre inequality in the separate audit,

```math
\operatorname{cap}(A\otimes R+I\otimes D)
\le k\lVert R\rVert_{\rm op}\Gamma(A)+n\operatorname{cap}(D),
```

is not an alternative proof using `W_1`.  Substituting (29) yields only a factor-two
bound through `W_1` in general.  A mean/residual or positive/negative
operator decomposition with the same shared pointwise norm budget is
precisely the state (26); repeating it would not evade either the
conference floor or the quotient-inheritance obstruction.

## 5. Finite tests after fixing the theorem

After (5) and the state (3) were fixed, I tested the correlation (`d=1`)
and degree-four Boolean moment (`d=2`) relaxations.  The matrices were read
without alteration from the repository files named below.  For each fixed
matrix, `cap` was recomputed exactly by enumerating the `2^(n-1)`
projective spins.  SDP entries are numerical CLARABEL 0.11.1 optima through
CVXPY 1.7.5 (absolute and relative gap tolerances `10^-8`); values within
`10^-6` of an integer cap are displayed as equal and are not claimed as
independent exact SDP certificates.

| seed | role | exact fixed-matrix `cap` | `W_1` | `W_1-cap` | `W_2` | `W_2-cap` |
|---|---|---:|---:|---:|---:|---:|
| `n=6` conference | certified minimizer | 5 | 6.708204 | 1.708204 | 5.000000 | 0 |
| `n=8`, orbit 0 | certified minimizer | 10 | 11.313708 | 1.313708 | 9.999999 | 0 |
| `n=8`, orbit 1 | certified minimizer | 10 | 12.000000 | 2.000000 | 10.000000 | 0 |
| `n=10` nonconference | certified minimizer | 13 | 17.946644 | 4.946644 | 13.078252 | 0.078252 |
| `n=10` GF(9) conference | proved construction | 15 | 15.000000 | 0 | 15.000000 | 0 |
| `n=12` saved witness | heuristic construction; fixed cap exact | 18 | 20.784610 | 2.784610 | 18.365820 | 0.365820 |
| `n=14` conference | proved construction | 21 | 25.238859 | 4.238859 | 25.238859 | 4.238859 |

Input provenance:

- `computations/results/exact_m6.json`;
- `computations/results/m8_minimizer_orbits.json`;
- `computations/results/exact_m10.json`;
- `computations/results/conference_order10_gf9.json`;
- `computations/results/heuristic_m12_seed20260731.json`;
- `computations/results/conference_completion_m13.json`.

The SDP used the Boolean moment matrix indexed by square-free monomials of
degree at most `d`, with entries `Y_(S,T)=y_(S triangle T)` and `y_empty=1`.
At `d=1`, the conference values in the table are additionally proved by
(20).  The `d=2` data show both sides of the issue: level two closes the
small-order gaps at `n=6,8`, but it already has a positive gap on the exact
nonconference `n=10` minimizer, a larger one on the saved `n=12` witness,
and no numerical improvement at all on the `n=14` conference signing.
These are diagnostics, not an asymptotic lower bound on the SOS gap.

## 6. Verdict and the smallest possible escape

Candidate JSA(`s,1`) is falsified at the desired leading scale on any
chiral below-half family by the proved floor (19).  At every fixed higher
level, JSA(`s,d`) is not a new absorption mechanism: (12)--(15) reduce it
to the assertion that degree `2d` already estimates the child tail within
`O_s(n)`.  The finite level-two tests do not prove that this assertion
fails asymptotically, but they show that it is nontrivial even on saved
exact minimizers and good structured seeds.

Accordingly, the polynomial joint-moment state is simpler to *evaluate*
than full parent maximization, but the needed uniform theorem is not yet
demonstrably simpler than the child bare-tail problem.  Exact quotient
compression forces every duplicated child pseudo-solution into the parent,
and seed dependence of all cross blocks does not change that fact.

Within this exact-compression/moment framework, there are three apparent
genuine escapes from this obstruction:

1. **A new seed-specific separator.**  Find a polynomial-size family of
   globally valid inequalities whose pullbacks close the seed relaxation
   gap to `O(n)`.  The pullback itself is the concrete theorem to prove;
   merely calling the constraints `B`-specific or operator-valued does not
   help.  Section 4.1 rules out every fixed-local marginal/facet family on
   sufficiently large conference seeds.
2. **Growing resolution.**  Use `d=d(n)` (or another state whose complexity
   grows) and prove project-scale accuracy while retaining enough
   tractability for composition.  This no longer gives the requested
   bounded-complexity state for free.
3. **Break the exact quotient.**  Replace (2) by a lift whose diagonal
   fibre subspace is not exactly the child Hamiltonian, while proving that
   the resulting quotient error has a globally cancelling, summable bound.
   This is a different composition theorem, not an exact compressed lift.

The smallest logically possible escape is item 1: one explicit,
polynomially describable child-valid separator that excludes the inherited
extremal pseudoexpectations uniformly.  Until such a separator is supplied,
a seed-dependent nonlocal lift plus a joint SDP is equivalent, at the
critical point, to certifying the original seed tail.
