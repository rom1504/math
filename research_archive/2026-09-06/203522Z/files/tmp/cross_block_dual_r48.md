# Wave 48C: projective completion profiles and cross-block pressure

## Status

- **Verified:** completion profiles from one signing form an exact projective
  system.  If `R subset S`, their slacks satisfy the Bellman identity
  `(R48C.3)` below.  In particular, active supports on overlapping blocks are
  projections of one global ground relation; abstract block profiles cannot
  be assigned independently.
- **Verified overlap-cover lemma:** a fractional edge-uniform cover by maximal
  selectors proves tight decomposition whenever the exact inequality
  `(R48C.7)` holds.  It certifies the stored `A6,m=5` and `A8,m=7` cases, but
  not `A9,m=8`; the maximal `A8,m=4` family has no such uniform cover.
- **Verified universal endpoint:** tight principal decomposition holds at
  `m=2,3` for every signing, without minimality.  For `m=3`, the second-moment
  lower bound and Mantel's theorem force a positive triangle inside every
  oriented parent ground.
- **Verified jointly realizable wall:** the exact minimizer `A8` has two bad
  maximal four-selectors.  Each separately has an exact zero one-block dual,
  but the two lifted dual faces are disjoint and their common minimax value is
  exactly `1/2`.  A four-edge flip shaves every old ground by eight but raises
  the old slack-eight layer by eight, leaving the norm equal to `20`.  Thus
  ground-face overlap pressure alone cannot prove descent.
- **Reproducible exact sample, not exhaustive:** twelve distinct
  random-objective MILP solutions at the proved order-ten cap `Q_10=26` all
  have an all-size common parent.  Either `20` or all `40` parent grounds work.
- **Open:** prove a cross-block inequality which retains the full slack
  profiles, not only their zero sets.  The exact falsifier to any proposed
  ground-layer rounding is the A8 migration in section 5.

The checker and output are `tmp/cross_block_dual_r48_check.py` and
`tmp/cross_block_dual_r48.out`.

## 1. The exact projective/Bellman constraint

For an oriented parent state `d`, let `E_d` be its full energy.  For an
oriented restriction `y` on `S`, define

```math
L_S(y):=k_S(y)+c_S(y)
       =\max_{d:d[S]=y}E_d,
\qquad
s_S(y):=q_n-L_S(y).
\tag{R48C.1}
```

If `R subset S`, optimizing first over the part outside `S` and then over
extensions from `R` to `S` gives the exact composition rule

```math
\boxed{
k_R(y_R)=
\max_{y_S:y_S|_R=y_R}
\{k_S(y_S)+c_S(y_S)\}-c_R(y_R).}
\tag{R48C.2}
```

Equivalently,

```math
\boxed{
s_R(y_R)=
\min_{y_S:y_S|_R=y_R}s_S(y_S).}
\tag{R48C.3}
```

Consequently

```math
\{s_R=0\}=\pi_{S\to R}\{s_S=0\}.
\tag{R48C.4}
```

Taking `S=[n]` says that every block-active state is the restriction of a
global parent ground.  For two overlapping selectors, both active supports
are therefore projections of the same global ground relation.  This is the
first exact constraint absent from the abstract antipodal profile in Wave 47.
The checker exhausts all `6561=3^8` nested pairs `R subset S` on `A8`.

The identity constrains supports and slack values, but it does **not** choose
probability laws on those supports.  In particular, the one-block dual laws
in `(10.1215)` can still depend on `S`.  Section 5 shows this dependence is
essential even on an exact minimizer.

## 2. An exact overlap-cover inequality

Fix `m`, put

```math
\mathcal F_m=
\{S:|S|=m,\ Q(A[S])=q_*\},
```

and give its members weights `lambda_S>=0`.  Write

```math
\Lambda=\sum_{S\in\mathcal F_m}\lambda_S,
\qquad
w_e=\sum_{S\ni e}\lambda_S.
\tag{R48C.5}
```

Suppose the selector weights are edge-uniform: `w_e=r` for every parent
edge.  For every oriented parent ground `d`, with
`M_{d,e}=a_ev_d(e)`, double counting gives

```math
\sum_{S\in\mathcal F_m}\lambda_Sc_S(d[S])
=2\sum_e w_eM_{d,e}
=r q_n.
\tag{R48C.6}
```

All oriented order-`m` energies have the same residue modulo four.  Hence if
tight decomposition fails, then for every parent ground and every maximal
selector

```math
c_S(d[S])\le q_*-4.
```

Combining this with `(R48C.6)` proves the following exact sufficient
condition:

```math
\boxed{
r q_n>\Lambda(q_*-4)
\quad\Longrightarrow\quad
\exists d\in\mathcal G_n(A),\ S\in\mathcal F_m:
c_S(d[S])=q_* .}
\tag{R48C.7}
```

After normalizing `Lambda=1`, edge counting forces

```math
r=\frac{\binom m2}{\binom n2}
=\frac{m(m-1)}{n(n-1)}.
\tag{R48C.8}
```

This is a genuine cross-block theorem: it uses the simultaneous internal
energies of one parent ground across a fractional selector cover, rather than
one-block minimality or a nested chain.

Exact examples from the checker are:

| instance | `r q_n` | `q_*-4` | conclusion |
|:--|--:|--:|:--|
| `A6,m=5` | `20/3` | `4` | proves a tight pair |
| `A8,m=7` | `15` | `14` | proves a tight pair |
| `A9,m=8` | `56/3` | `20` | inequality misses |
| `A8,m=4` | -- | `8` | maximal blocks have no edge-uniform cover |

The last nonexistence is exact, not a solver-tolerance claim.  In the
lexicographic selector order printed by the checker, the ten weights reduce
under the two-selector edge equations to

```math
\lambda_0=\lambda_3=\lambda_4=\lambda_6=\lambda_8=:a,
\qquad
\lambda_1=\lambda_2=\lambda_5=\lambda_7=\lambda_9=:b.
```

The edge `01` then has load `a+b`, edge `05` has load `3a`, and edge
`16` has load `3b`.  Constant load would give
`a+b=3a=3b`, hence `b=2a` and `a=2b`, forcing `a=b=0`.
Thus no nonzero constant-load weighting exists.

At fixed density the coefficient in `(R48C.8)` is approximately `p^2`,
whereas the natural principal scale is `p^(3/2)`.  Thus `(R48C.7)` is not by
itself a scalable fixed-density proof.  Its likely use is near the diagonal
or after a stronger minimizer-specific estimate replaces the uniform-load
identity.

## 3. Tight decomposition is universal at selector size three

Every signed triangle has full-matrix norm `6`, so `q_*=6` for `m=3`.
Fix any signing `A` of order `n>=3` and any oriented parent ground `d`.
Put

```math
b_{ij}=a_{ij}v_d(ij)\in\{-1,1\},
\qquad
P=\sum_{i<j}b_{ij}=q_n/2,
\qquad
N=\binom n2.
\tag{R48C.9}
```

For uniform Boolean `x`, orthogonality of distinct quadratic characters gives

```math
\mathbb E(x^{\mathsf T}Ax)^2=4N,
```

and therefore

```math
P=q_n/2\ge\sqrt N>n/2.
\tag{R48C.10}
```

The graph formed by the positive `b`-edges has

```math
e_+=\frac{N+P}{2}>\frac{n^2}{4}.
\tag{R48C.11}
```

Mantel's theorem supplies a positive triangle.  On that selector all three
terms `b_ij` are `+1`, so its inherited oriented energy is `6=q_*`.  Thus:

> For every signing of order at least three, every oriented parent ground has
> a ground restriction on some three-selector; hence `(10.1203)` holds at
> `m=3` without using exact signing minimality.

At `m=2`, the same conclusion is immediate because `P>0` supplies a positive
edge and every two-selector has norm `2`.  The checker exhausts every
switching-normalized signing through order six (`2,8,64,1024` signings at
orders `3,4,5,6`) and every parent ground as a finite audit of the theorem.

## 4. Lifted one-block dual faces

For a selector `S`, every zero-slack child state lifts to at least one global
ground.  Thus the zero-dual face `(10.1215)` can be lifted into the common
global ground simplex:

```math
\mathcal D_S=
\left\{\mu\in\Delta(\mathcal G_n(A)):
\mathbb E_\mu M_{d,e}\le0\quad(e\in E(S))\right\}.
\tag{R48C.12}
```

The one-block fractional margin is zero exactly when `D_S` is nonempty.
Cross-block consistency would ask for an intersection of several such faces.
The Bellman identity makes all the faces live over one genuine ground
relation, but it does not imply that their intersection is nonempty.  The
next example shows that even two actual bad maximal blocks can require
mutually singular laws.

There is an exact complement-size restriction on these zero duals.  Put

```math
K_S=|E(K_n)\setminus E(S)|
=\binom n2-\binom m2.
```

For an arbitrary law `p` in the dual `(10.1215)`, lift each child state `y`
to a completion attaining `k_S(y)`.  That lift has full energy `q_n-s_y`.
Since each of its `K_S` external signed edge terms is at most one,

```math
\sum_{e\in E(S)}\mathbb E_pM_{y,e}
\ge \frac{q_n-\mathbb E_ps}{2}-K_S.
\tag{R48C.12a}
```

The sum of coordinatewise positive parts dominates the positive part of the
sum.  Substitution in `(10.1215)` and minimization over
`t=mathbb E_p s>=0` therefore gives

```math
\boxed{
\gamma_{\rm frac}(S)
\ge\min_{t\ge0}
\left\{t+4\left(\frac{q_n-t}{2}-K_S\right)_+\right\}
=(q_n-2K_S)_+.}
\tag{R48C.12b}
```

The same bound has a direct primal interpretation: take the half-flip
`f_e=1/2` on every internal edge.  Its margin at `y` is

```math
s_y+2\sum_{e\in E(S)}M_{y,e}=s_y+c_S(y)=q_n-k_S(y)
\ge q_n-2K_S.
```

Thus `(R48C.12b)` is precisely the fractional gain from zeroing the whole
internal block; it need not survive Boolean rounding.

In particular, a zero-dual antipodal obstruction requires

```math
\boxed{q_n/2\le K_S.}
\tag{R48C.12c}
```

For `m=n-r`, `K_S=r(2n-r-1)/2`.  Since exact minimizers satisfy the verified
lower bound `q_n=Omega(n^(3/2))`, `(R48C.12c)` is impossible when
`r=o(sqrt(n))`.  Hence the zero-fractional-margin obstruction from Wave 47 is
confined away from the very near diagonal.  This does **not** prove tight
decomposition there: exact minimality still has `gamma_bool=0`, and a
positive fractional margin can be lost entirely in Boolean rounding.  The
checker records the exact lower bounds `0`, `6`, and `8` for `A8,m=4`,
`A8,m=7`, and `A9,m=8`, respectively.

## 5. Exact A8 cross-block incompatibility and slack migration

The stored exact minimizer `A8` has

```math
q_8=20,
\qquad
q_*^{(4)}=12,
```

and ten maximal four-selectors.  Exactly two selectors have no parent ground
whose restriction is a child ground:

```math
S_0=\{0,3,4,5\},
\qquad
S_1=\{1,2,6,7\}=S_0^c.
\tag{R48C.13}
```

Across the eight oriented parent grounds, the inherited energy on either bad
block is always `0` or `4`, never `12`.  For each `S_i`, exactly four parent
grounds have inherited energy zero.  Uniform mass on those four grounds has

```math
\mathbb E M_{d,e}=0
\qquad(e\in E(S_i)),
\tag{R48C.14}
```

so it is an exact zero dual in `(R48C.12)`.  However, the two four-element
zero-energy ground sets are complementary.  Any zero dual for `S_i` must be
supported on its zero-energy set: every inherited energy is nonnegative, and
the sum of the nonpositive edge means must also be nonpositive.  Hence

```math
\boxed{\mathcal D_{S_0}\cap\mathcal D_{S_1}=\varnothing.}
\tag{R48C.15}
```

The incompatibility has a short exact pressure certificate.  Let

```math
F=\{05,34,16,27\}.
```

For every oriented parent ground `d`, direct integer evaluation gives

```math
\boxed{\sum_{e\in F}M_{d,e}=2.}
\tag{R48C.16}
```

Therefore every law `mu` on the parent grounds has some edge in
`E(S_0) union E(S_1)` with mean at least `1/2`.  Uniform mass on all eight
grounds has every such edge mean at most `1/2`, so the exact common-law game
is

```math
\boxed{
\min_{\mu\in\Delta(\mathcal G_8(A_8))}
\max_{e\in E(S_0)\cup E(S_1)}\mathbb E_\mu M_{d,e}
=\frac12.}
\tag{R48C.17}
```

This positive cross-block pressure still does not give an improving signing.
Flip the four edges in `F`.  Formula `(R48C.16)` lowers every old parent
ground from energy `20` to energy `12`.  But exactly eight old oriented states
with energy `12`, slack `8`, and flip sum `-2` rise to energy `20`.  Thus

```math
Q(A_8^F)=20=Q(A_8).
\tag{R48C.18}
```

This is an exact, genuinely jointly realizable obstruction to a
ground-face-only cross-block proof.  It is not a counterexample to tight
principal decomposition: the other eight maximal four-selectors contain
tight pairs.  Its implication is narrower and important: compatible or
positive dual information on the current ground face must be strengthened to
control the nearby slack layers before it can be rounded to a Boolean block
replacement.

An exact sufficient target for a proposed cross-block flip set `F` is still

```math
\boxed{
s_{[n]}(d)+4\sum_{e\in F}M_{d,e}>0
\quad\text{for every oriented state }d.}
\tag{R48C.19}
```

The A8 falsifier attains equality in `(R48C.19)` on its migrated slack-eight
layer.  Any successor theorem should therefore state how overlap information
controls positive slack levels, not just the active projections in
`(R48C.4)`.

## 6. Reproducible order-ten sample

The active ledger records the proved optimum `Q_10=26`.  The checker uses the
exact cap `|H_A(x)|<=13`, fixes the switching gauge at vertex zero, and solves
twelve deterministic random linear objectives over the resulting binary
polytope.  Exact integer reevaluation gives twelve distinct normalized
signings of norm `26`.

Every sampled minimizer has `40` oriented parent grounds and at least one
ground that works separately at every selector size.  The exact results are:

| maximal-principal profile `q_*^(m)`, `m=0,...,10` | samples | all-size parents |
|:--|--:|:--|
| `(0,0,2,6,12,16,22,22,24,24,26)` | `7` | `40/40` |
| `(0,0,2,6,12,20,22,22,24,24,26)` | `5` | `20/40` |

This is not an enumeration of order-ten minimizers.  It extends the exact
finite sample supporting the all-size-common-parent strengthening, while
leaving the conjecture asymptotically open.

## 7. Research judgment

The Bellman identity and `(R48C.7)` are the positive outputs.  They identify
two legitimate cross-block mechanisms: projective support/slack consistency
and fractional selector covers.  The universal `m=3` theorem shows that
global ground energy plus a genuinely multi-selector extremal theorem can
force a tight pair.

The A8 calculation sharply limits the route.  Separate one-block zero duals
need not align even when their profiles are jointly realizable, and forcing a
positive direction on every current ground does not control witness migration.
Thus the next plausible lemma must combine several maximal blocks at the
**full slack-profile level**, for example by deriving `(R48C.19)` from a
weighted selector cover with a quantitative penalty for low positive slack.
Without such a penalty, cross-block dual consistency is another fractional
ground-layer relaxation and cannot prove `(10.1203)`.
