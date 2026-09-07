# Audit of the proposed depth-uniform tropical/cohomological theorem

## Scope and verdict

This audit uses the metric-shell closure and response convention of Theorem
16.3, the long-depth obstruction in Theorem 16.7, and the weighted-automaton
and exact-lumpability conventions of Theorem 16.9 and the weighted-automaton
benchmark draft.

The proposed result is sound only after the following qualifications.

1. **Metric shells with gauges:** the composition formula is exact.  Its
   orientation is
   `+phi_(t-1)` at the source and `-phi_t` at the target.  All internal
   potentials cancel.  The terminal potential does not enter a directed
   row difference at all; the error in a directed row response is exactly
   `phi_0(a)-phi_0(b)`, and is therefore at most `osc(phi_0)`.  If instead
   one measures the whole endpoint kernel modulo one global scalar, both
   endpoint oscillations enter, with best sup error
   `(osc(phi_0)+osc(phi_T))/2`.

2. **Defect characterization:** zero rectangular defects of each residual
   are necessary but not sufficient.  At every shared interface one must
   additionally require

   ```math
   E_t(a,x)+E_{t+1}(x,c)
   =E_t(a,x')+E_{t+1}(x',c).                       \tag{A.cross}
   ```

   For nonempty complete finite interfaces, the rectangular conditions and
   (A.cross) are jointly necessary and sufficient for simultaneous boundary
   potentials.  For sparse supports, rectangles are not the right complete
   test; all alternating support cycles must be tested.

3. **Finite graph cocycles:** bounded path sums are equivalent to zero
   holonomy on every directed cycle only when the paths are all paths in the
   stated finite graph.  For a declared regular language, the graph must be
   the finite automaton/product graph encoding that language, and only
   reachable and co-reachable (hence pumpable) cycles matter.  A nonzero
   holonomy on such a cycle gives linear drift.  The coboundary is required
   only inside strongly connected components; labels on the acyclic edges
   between components can be arbitrary.

4. **Max-plus projective coefficient:** on the *full finite-dimensional*
   projective space `R^n/R1`, an all-finite max-plus matrix has global
   projective Lipschitz coefficient either zero or one.  Zero is equivalent
   to additive rank one, `A_(ij)=u_i+v_j`, and to constant projective image.
   The assertion is false on an arbitrary finite or restricted reachable
   subset.

5. **Syndetic resets:** if a reset word has length `m`, its quotient map has
   projective image diameter at most `rho`, and the last completed occurrence
   is at most `G` letters behind, the correct bound is

   ```math
   rho+(m+G)epsilon.                               \tag{R.bound}
   ```

   Thus `rho+2L epsilon` follows when both `m<=L` and `G<=L`, with the first
   reset also occurring in bounded time.  "Gaps at most `L`" by itself is
   ambiguous and does not imply all of these facts.

The detailed proofs and a ready-to-use theorem statement follow.

## 1. Exact metric-shell gauge cancellation

Let `(X_t,d_t)`, `0<=t<=T`, be nonempty finite metric spaces.  Let
`g_t:X_(t-1)->X_t` be a **bijective** isometry and let `lambda_t>=0`.  Define

```math
D_t(a,b)=D_{lambda_t,g_t}(a,b)
       =lambda_t d_t(b,g_t(a)).                    \tag{1.1}
```

Surjectivity of each `g_t` is needed by the endpoint-attainment proof of the
metric-shell composition law.  An isometric embedding is not enough in
general.

Suppose

```math
K_t(a,b)=D_t(a,b)+phi_{t-1}(a)-phi_t(b)+c_t,       \tag{1.2}
```

where `phi_t:X_t->R` and `c_t in R`.  With min-plus composition

```math
(K star L)(a,c)=min_b {K(a,b)+L(b,c)},             \tag{1.3}
```

put

```math
lambda_* = min_{1<=t<=T} lambda_t,
G=g_T circ ... circ g_1,
C=sum_{t=1}^T c_t.
```

Then

```math
K_1 star ... star K_T
=D_{lambda_*,G}+phi_0-phi_T+C,                    \tag{1.4}
```

where the endpoint functions on the right mean
`phi_0(a)-phi_T(z)`.

### Proof

At one internal interface,

```math
\begin{aligned}
(K_t star K_{t+1})(a,c)
={}&phi_{t-1}(a)-phi_{t+1}(c)+c_t+c_{t+1}\\
 &+min_x\{D_t(a,x)+D_{t+1}(x,c)\}.
\end{aligned}                                     \tag{1.5}
```

The two occurrences of `phi_t(x)` cancel before minimization.  Since
`g_(t+1)` is a bijective isometry, the triangle inequality and the two
endpoint choices `x=g_t(a)` and `x=g_(t+1)^(-1)(c)` give

```math
D_{lambda_t,g_t} star D_{lambda_{t+1},g_{t+1}}
=D_{min(lambda_t,lambda_{t+1}),g_{t+1} circ g_t}. \tag{1.6}
```

Induction proves (1.4).  Notice that this is exact and has no term
proportional to `T`.

The constants `c_t` are operationally useful for keeping track of absolute
baselines, but on a single acyclic time chain they are algebraically
redundant: replacing

```math
phi_t by phi_t-sum_{s=1}^t c_s
```

sets every `c_t` to zero.  They cannot all be removed by one stationary
potential on a graph with cycles unless the corresponding scalar cycle
holonomies vanish; projective responses, however, ignore such scalar
baselines.

### Directed response and endpoint constants

Use the directed response from Theorem 16.3,

```math
r(f,g)=max_z(f(z)-g(z)).                           \tag{1.7}
```

For `K=K_1 star ... star K_T`, (1.4) gives, for every ordered source pair,

```math
\boxed{
r(K(a,cdot),K(b,cdot))
=lambda_* d_0(a,b)+phi_0(a)-phi_0(b).}            \tag{1.8}
```

Indeed, `-phi_T(z)+C` cancels pointwise between the two rows, and the
metric-shell maximum is attained at `z=G(b)`.  Consequently

```math
|r(K(a,cdot),K(b,cdot))-lambda_*d_0(a,b)|
=|phi_0(a)-phi_0(b)|
<=osc(phi_0).                                     \tag{1.9}
```

The reverse directed response has the opposite potential difference.  If
rows are themselves taken projectively, with

```math
d_pr(f,g)=inf_c ||f-g-c1||_infinity
         ={1\over2}osc(f-g),                      \tag{1.10}
```

then all row constants disappear and

```math
d_pr(K(a,cdot),K(b,cdot))=lambda_*d_0(a,b).       \tag{1.11}
```

On the other hand, the residual of the entire endpoint kernel relative to
the ideal shell is

```math
R(a,z)=phi_0(a)-phi_T(z)+C.
```

Therefore

```math
osc_{a,z}R=osc(phi_0)+osc(phi_T),
\qquad
inf_c||R-c1||_infinity
={osc(phi_0)+osc(phi_T)\over2}.                   \tag{1.12}
```

Equations (1.9) and (1.12) are different claims.  A theorem should not call
both merely "endpoint oscillation" without specifying the response metric.

## 2. Exact characterization by rectangular and interface defects

Let `E_t:X_(t-1) times X_t -> R` be arbitrary residuals, with every `X_t`
nonempty.  Define their rectangular defects by

```math
\begin{aligned}
Delta_t(a,a';b,b')={}&E_t(a,b)+E_t(a',b')\\
                    &-E_t(a,b')-E_t(a',b),        \tag{2.1}
\end{aligned}
```

and their cross-interface defects by

```math
\begin{aligned}
Gamma_t(a,c;x,x')={}&E_t(a,x)+E_{t+1}(x,c)\\
                    &-E_t(a,x')-E_{t+1}(x',c),    \tag{2.2}
\end{aligned}
```

for `1<=t<T`.

### Proposition

The following are equivalent.

1. There are functions `phi_t:X_t->R` and constants `c_t` such that

   ```math
   E_t(a,b)=phi_{t-1}(a)-phi_t(b)+c_t             \tag{2.3}
   ```

   simultaneously for every `t`.

2. Every `Delta_t` and every `Gamma_t` is zero.

### Proof

Substitution of (2.3) immediately makes (2.1) and (2.2) zero.

Conversely, `Delta_t=0` is the standard additive-separability criterion.
After choosing one source and target anchor at each layer, it gives

```math
E_t(a,b)=A_t(a)+B_t(b).                            \tag{2.4}
```

Condition `Gamma_t=0` says exactly that

```math
B_t(x)+A_{t+1}(x)=kappa_t                         \tag{2.5}
```

is constant on `X_t`.  Set `phi_0=A_1`, `phi_t=-B_t` for `t>=1`, take
`c_1=0`, and for `t>=2` take `c_t=kappa_(t-1)`.  Equations (2.4)-(2.5) then
give (2.3).  This proves sufficiency.

Rectangular defects alone do not suffice.  For example, on a two-point
shared interface take `E_1(a,x)=0` and `E_2(x,c)=q(x)` for a nonconstant
`q`.  Both residuals have zero rectangular defects.  The first residual
forces the shared potential to be constant, while the second forces it to
vary like `q`; equivalently, (2.2) equals `q(x)-q(x')`.

This proposition assumes a full real-valued matrix on every consecutive
pair of layers.  On a sparse bipartite support, four-cycles may not exist or
generate the support cycle space.  Additive separability then requires zero
alternating sum around every support cycle, followed by the analogous
compatibility test on shared vertices.

For a branching family of declared words, (2.2) must hold on every allowed
adjacency in the actual finite context graph.  Compatibility along one
chosen word does not by itself produce a word-independent potential.  The
remaining obstruction around repeated contexts is precisely the cycle
holonomy in the next section.

## 3. Finite-graph cocycles and arbitrary depth

Let `G=(V,E)` be a finite directed multigraph and let
`omega:E->W`, where `W` is a real normed vector space.  For a directed path
`P=e_1...e_k`, write

```math
Hol(P)=sum_{j=1}^k omega(e_j).                     \tag{3.1}
```

No finite-dimensional hypothesis on `W` is actually needed.

### Proposition

If all finite directed paths in `G` are admissible, the following are
equivalent.

1. `sup_P ||Hol(P)||<infinity`.
2. `Hol(C)=0` for every directed closed walk (equivalently, every simple
   directed cycle) `C`.
3. For every strongly connected component `S` there is a potential
   `q_S:S->W` such that every edge internal to `S` satisfies

   ```math
   omega(u->v)=q_S(v)-q_S(u).                     \tag{3.2}
   ```

When these conditions hold, one explicit bound is

```math
||Hol(P)|| <= (|V|-1) max_{e in E}||omega(e)||.   \tag{3.3}
```

If a cycle `C` has `h=Hol(C) ne 0`, then

```math
Hol(C^n)=nh,
\qquad ||Hol(C^n)||=n||h||.                       \tag{3.4}
```

With a fixed entrance path `P` and exit path `Q`, the more operational form
is

```math
||Hol(PC^nQ)||
>=n||h||-||Hol(P)+Hol(Q)||.                       \tag{3.5}
```

### Proof

Boundedness applied to `C^n` implies (2).  To obtain (3), fix a root `r` in
an SCC and define `q_S(v)` as the sum along any directed path from `r` to
`v`.  If `P,P'` are two such paths, choose a directed path `Q` from `v`
back to `r`; both `PQ` and `P'Q` are closed, so their zero sums imply that
`P` and `P'` have equal sums.  Thus `q_S` is well-defined, and appending an
edge proves (3.2).

Alternatively, (2) directly implies boundedness: whenever a path repeats a
vertex, delete the intervening closed subpath.  Its sum is zero.  Repeating
this operation leaves a simple path with the same sum and at most
`|V|-1` edges, proving (3.3).  Formulae (3.4)-(3.5) are immediate.

Edges between distinct SCCs need not be coboundaries.  The condensation
graph is acyclic, so a path uses at most `#SCC-1` such edges.  This is why
"a coboundary on each SCC" is the correct statement and "a coboundary on
the entire graph" is unnecessarily strong.

### Declared-word version

For paths accepted from an initial set `I` to a terminal set `F`, replace
`G` by its useful subgraph: vertices and edges reachable from `I` and
co-reachable to `F`.  Then accepted path sums are uniformly bounded if and
only if every directed cycle in a useful SCC has zero holonomy.  A useful
cycle is pumpable: there are fixed paths `P,Q` for which `PC^nQ` is accepted
for every `n`, so (3.5) proves linear drift when its holonomy is nonzero.

For a regular declared word language, any extra word/context constraint
must first be incorporated into a finite automaton or product graph.  A
cycle in an underlying raw-state graph that is unreachable, not
co-reachable, or forbidden by the declared language is not an obstruction.
Conversely, for an arbitrary non-regular subset of graph paths, the simple
SCC formulation need not characterize what that subset can repeat.

If responses are projective, apply the proposition in the quotient vector
space `W/span{1}`.  Then the required cycle holonomy is zero **modulo a
common scalar**.  A scalar cycle drift is invisible projectively but is a
real linear drift for absolute weighted-automaton outputs.

## 4. The max-plus projective Lipschitz dichotomy

For `x,y in R^n`, let

```math
d_H([x],[y])=max_j(x_j-y_j)-min_j(x_j-y_j),       \tag{4.1}
```

so the metric used in Theorem 16.9 is `d_pr=d_H/2`.  Let an all-finite
matrix `A in R^(m times n)` act by

```math
(F_Ax)_i=max_j(A_(ij)+x_j).                       \tag{4.2}
```

Define the global coefficient on the full projective domain by

```math
Lip_H(F_A)=sup_{[x] ne [y]}
 {d_H(F_A[x],F_A[y])\over d_H([x],[y])}.          \tag{4.3}
```

### Proposition (zero-or-one law)

For an all-finite max-plus matrix on the full `R^n/R1`,

```math
\boxed{
Lip_H(F_A)=
\begin{cases}
0,&A_(ij)=u_i+v_j\text{ for some }u,v,\\
1,&\text{otherwise}.
\end{cases}}                                     \tag{4.4}
```

The first case is exactly constant projective image.  It is safest to call
it **additive rank one** or **max-plus factor rank one**, since "tropical
rank" has several conventions.

### Proof

For arbitrary `x,y`, put

```math
alpha=min_j(x_j-y_j),\qquad beta=max_j(x_j-y_j).
```

Monotonicity and additive homogeneity give, for every output coordinate,

```math
alpha <= (F_Ax)_i-(F_Ay)_i <= beta.               \tag{4.5}
```

Hence every max-plus matrix is projectively nonexpansive and
`Lip_H(F_A)<=1`.

If `A_(ij)=u_i+v_j`, then

```math
F_Ax=u+max_j(v_j+x_j)1,                           \tag{4.6}
```

so its projective image is the single point `[u]` and its coefficient is
zero.  Conversely, because all entries are finite, each column `j` can be
made the unique maximizer in every row by taking `x_j` sufficiently larger
than all other coordinates.  Constant projective image therefore forces
all columns of `A` to be projectively equal, which is precisely
`A_(ij)=u_i+v_j`.

Now suppose `A` is not additively rank one.  There are rows `i,k` and
columns `j,l` with

```math
A_(ij)+A_(kl) ne A_(il)+A_(kj).                  \tag{4.7}
```

Suppress all other input coordinates and choose `x_j-x_l` strictly between
the two row thresholds

```math
A_(il)-A_(ij)\quad\hbox{and}\quad
A_(kl)-A_(kj).                                    \tag{4.8}
```

On an open neighbourhood of this `x`, row `i` uniquely selects one of
`j,l` and row `k` uniquely selects the other.  Perturb `x_j` by `+s` and
`x_l` by `-s`, with `s>0` small enough to stay in that selector cell, and
leave all other perturbations between these two values.  The input Hilbert
distance is `2s`; the increments of output rows `i,k` are `+s,-s`, so the
output Hilbert distance is at least `2s`.  Together with nonexpansiveness it
is exactly `2s`, proving `Lip_H(F_A)=1`.

The all-finite and full-domain qualifications matter.  For example,

```math
A_delta=\begin{pmatrix}0&0\\-delta&0\end{pmatrix},
\qquad delta>0,                                   \tag{4.9}
```

acts on the projective coordinate `z=x_1-x_2` by

```math
z |-> max(z,0)-max(z-delta,0)
   =min(max(z,0),delta).                           \tag{4.10}
```

Its full projective image has `d_H`-diameter `delta` (and `d_pr`-diameter
`delta/2`), yet its global coefficient is one because it has slope one on
`(0,delta)`.  This explicitly shows that arbitrarily small nonzero image
diameter does not imply a coefficient below one.

On the restricted two-point domain `{z=-delta,z=2delta}`, however, the
Lipschitz ratio is

```math
{delta\over3delta}={1\over3}.                     \tag{4.11}
```

Thus the zero-or-one claim is false if "finite projective space" means a
finite set of reachable projective states, rather than the full
finite-dimensional tropical projective space.

## 5. Depth-uniform error from a syndetic reset word

The reset assertion is a metric fact and does not follow from a strict
Lipschitz coefficient.  Here is a formulation with all needed hypotheses.

Let `F_a:X->X` be the concrete transitions, let `Q_a:Y->Y` be
nonexpansive maps of a metric space `(Y,d)`, and let `pi:X->Y`.  Max-plus
projective maps with either `d_H` or `d_pr` are examples.  Suppose the local
semiconjugacy defect is measured in that same projective metric and obeys

```math
d(pi F_a(x),Q_a pi(x))<=epsilon                   \tag{5.1}
```

for every allowed letter and every relevant reachable `x`.  For a word
`u=a_1...a_t`, composition is in time order,

```math
F_u=F_(a_t)circ...circ F_(a_1),\qquad
Q_u=Q_(a_t)circ...circ Q_(a_1).                  \tag{5.2}
```

Assume a word `w` of length `m` has reset diameter

```math
diam(Q_w(Y))<=rho.                                \tag{5.3}
```

If a declared word factors as `u=pwq`, then

```math
\boxed{
sup_x d(pi F_u(x),Q_u pi(x))
<=rho+(m+|q|)epsilon.}                            \tag{5.4}
```

### Proof

Let `y_s` be the projection of the concrete orbit and `z_s` the quotient
orbit, both started with `y_0=z_0=pi(x)`.  Nonexpansiveness and (5.1) give
the usual telescoping estimate over any `r` consecutive letters:

```math
d(y_(s+r),Q_(block)y_s)<=r epsilon.               \tag{5.5}
```

At the end of the displayed occurrence of `w`,

```math
\begin{aligned}
d(y_(|p|+m),z_(|p|+m))
&<=d(y_(|p|+m),Q_wy_(|p|))\\
&\quad+d(Q_wy_(|p|),Q_wz_(|p|))\\
&<=m epsilon+rho.                                 \tag{5.6}
\end{aligned}
```

Propagating through the trailing word `q` adds at most `|q|epsilon`, which
proves (5.4).  Crucially, the diameter in (5.3) is the diameter of the
**quotient** reset map on a set containing both quotient trajectories.  If
only the concrete reset has small diameter, additional lifting and
regularity assumptions are needed.

A precise syndetic hypothesis is the following: every declared prefix is
either of length at most `m+G`, or has a factorization `pwq` with
`|q|<=G`.  (Equivalently in an infinite-run formulation, completed reset
ends have gaps at most `G` and the first completion occurs within bounded
time.)  The direct telescoping bound handles the short prefixes, and (5.4)
handles the others, yielding

```math
sup_{u,x}d(pi F_u(x),Q_u pi(x))
<=rho+(m+G)epsilon.                               \tag{5.7}
```

If `m+G<=2L`--in particular, if `m<=L` and `G<=L`--the advertised coarse form

```math
sup_{u,x}d(pi F_u(x),Q_u pi(x))
<=rho+2L epsilon                                  \tag{5.8}
```

follows from these symbols.  A bound on gaps between occurrences does not,
under every convention, bound the reset-word length or the first occurrence;
those facts must be declared.  The word family should be prefix-closed, or
all assumptions should explicitly be imposed on its prefix closure, because
(5.1) is used at intermediate states.

There is no hidden factor of two in (5.7) when `rho` and `epsilon` are both
measured in `d_H`, or both in `d_pr`.  If instead the local datum is a sup
error `||r-s||_infinity<=epsilon` between fixed representatives, then it
implies a `d_pr` error at most `epsilon` but only a `d_H` error at most
`2epsilon`; the theorem's constants must be converted accordingly.

## 6. Recommended exact theorem statement

The following is suitable as a single theorem, with the three mechanisms
kept logically distinct.

### Theorem (boundary gauges, repeatable holonomy, and tropical resets)

**(a) Gauged metric shells.**  Let `(X_t,d_t)_(t=0)^T` be nonempty finite
metric spaces, `g_t:X_(t-1)->X_t` bijective isometries, and `lambda_t>=0`.
Put

```math
D_t(a,b)=lambda_t d_t(b,g_ta).
```

If

```math
K_t(a,b)=D_t(a,b)+phi_(t-1)(a)-phi_t(b)+c_t,
```

then, for min-plus composition,

```math
K_1 star ... star K_T
=D_(min_t lambda_t,g_T circ ... circ g_1)
 +phi_0-phi_T+sum_t c_t.                          \tag{6.1}
```

For `r(f,g)=max(f-g)`, its ordered row responses satisfy

```math
r(K_(1:T)(a,cdot),K_(1:T)(b,cdot))
=(min_t lambda_t)d_0(a,b)+phi_0(a)-phi_0(b),      \tag{6.2}
```

so the uniform directed-response error is at most `osc(phi_0)` and has no
depth factor.  Projective row shape is exactly `(min_t lambda_t)d_0(a,b)`.

For residuals `E_t=K_t-D_t` on complete interfaces, such simultaneous
potentials and constants exist if and only if

```math
E_t(a,b)+E_t(a',b')-E_t(a,b')-E_t(a',b)=0         \tag{6.3}
```

for every rectangle and

```math
E_t(a,x)+E_(t+1)(x,c)
-E_t(a,x')-E_(t+1)(x',c)=0                        \tag{6.4}
```

at every shared interface.

**(b) Repeatable cocycles.**  Let a finite directed graph carry edge labels
in a normed real vector space.  Sums over all finite directed paths are
uniformly bounded if and only if every directed cycle has zero label sum.
Equivalently, on each SCC the edge label is a vertex coboundary.  A cycle
with holonomy `h ne 0` has `n`-fold sum `nh` and hence linear norm drift.
For accepted/declared words this assertion applies to the reachable,
co-reachable part of the finite automaton or product graph encoding the
declaration; only cycles pumpable in that graph are relevant.  For
projective responses, all labels and holonomies are taken modulo common
scalars.

**(c) Max-plus dichotomy and resets.**  Let `A` be an all-finite max-plus
matrix acting on the full finite-dimensional projective space.  Its global
Hilbert-projective Lipschitz coefficient is zero if
`A_(ij)=u_i+v_j`, and is one otherwise.  Thus finite nonzero projective image
diameter does not yield a strict global contraction coefficient.

Nevertheless, let `pi F_a` and `Q_a pi` have one-letter defect at most
`epsilon` in one fixed projective metric, with every `Q_a` nonexpansive.  If
a word `w` of length `m` satisfies `diam(Q_w(Y))<=rho`, and every declared
prefix either has length at most `m+G` or ends at most `G` letters after a
completed occurrence of `w`, then

```math
sup_{u,x}d(pi F_u(x),Q_u pi(x))
<=rho+(m+G)epsilon.                               \tag{6.5}
```

In particular, if `m,G<=L`, the right side is at most
`rho+2L epsilon`.

### Final assessment

Part (a) is an exact telescoping theorem, part (b) is the precise
cohomological obstruction to extending a bounded local correction through
arbitrary repetitions, and part (c) is a distinct finite-memory reset
mechanism.  They fit naturally in one depth-uniform statement, but none of
the connectivity, quotient-metric, or declared-word qualifications should
be suppressed: doing so makes the graph equivalence false and makes the
`rho+2L epsilon` constant unjustified.
