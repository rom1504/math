# Depth-uniform tropical response stability

**Status.**  Independently audited theorem draft.  The algebraic claims below
are elementary.  Their purpose is to isolate two finite
mechanisms that really prevent long-depth response drift: an exact gauge
coboundary, or a recurrent small-diameter reset.  Small one-step defect and
ordinary nonexpansiveness do not suffice.

## 1. Metric shells with interface gauges

Let `(Y,d)` be a finite metric space.  For a bijective isometry `g` and
`lambda>=0`, put

```math
D_{\lambda,g}(a,b)=\lambda d(b,g(a)).             \tag{DU.1}
```

Composition is min-plus:

```math
(K\star L)(a,c)=\min_b\{K(a,b)+L(b,c)\}.          \tag{DU.2}
```

The metric-shell identity is

```math
D_{\lambda,g}\star D_{\mu,h}
=D_{\min(\lambda,\mu),h\circ g}.                 \tag{DU.3}
```

### Theorem DU.1 (gauge-telescoping shell algebra)

For `t=1,...,T`, choose interface potentials `phi_t:Y->R`, a starting
potential `phi_0`, and scalars `c_t`.  If

```math
K_t(a,b)=D_{\lambda_t,g_t}(a,b)
          +\phi_{t-1}(a)-\phi_t(b)+c_t,           \tag{DU.4}
```

then

```math
\boxed{
K_1\star\cdots\star K_T
=D_{\lambda_*,G}+\phi_0\otimes\boldsymbol1
 -\boldsymbol1\otimes\phi_T+\sum_{t=1}^Tc_t,}   \tag{DU.5}
```

where

```math
\lambda_*=\min_t\lambda_t,
\qquad G=g_T\circ\cdots\circ g_1.               \tag{DU.6}
```

In particular the directed row response is

```math
\max_b\{K_{1:T}(a,b)-K_{1:T}(a',b)\}
=\lambda_*d(a,a')+\phi_0(a)-\phi_0(a').          \tag{DU.7}
```

Thus a potential perturbation can be present at every layer without paying
once per layer.  All internal potentials cancel before the absolute value or
directed maximum is taken; only the source endpoint potential remains in the
row table.  If `osc(phi_0)<=delta`, the directed table differs from that of
the unperturbed shell by at most `delta`, independently of `T`.
Projective distance between any two rows is exactly
`lambda_*d(a,a')`, since the source potential is a row constant.  By
contrast, the whole endpoint kernel differs projectively from its ideal
shell by

```math
{\operatorname {osc}(\phi_0)+\operatorname {osc}(\phi_T)\over2}. \tag{DU.7a}
```

These are different operational metrics and should not be conflated.

#### Proof

In the sum defining `K_t star K_(t+1)`, the terms `-phi_t(b)` and
`+phi_t(b)` cancel before minimization.  Induction and (DU.3) give (DU.5).
The terminal potential and all scalar baselines cancel between two rows.
The reverse triangle inequality gives the upper bound in (DU.7), and the
centre of the second row attains it, exactly as for the unperturbed distance
kernel. `square`

For arbitrary terminal fields, (DU.5) also gives

```math
\min_b\{K_{1:T}(a,b)+h(b)\}
=\phi_0(a)+\sum_tc_t+
 \min_b\{D_{\lambda_*,G}(a,b)+h(b)-\phi_T(b)\}.  \tag{DU.8}
```

Thus the terminal gauge reparametrizes the future field rather than creating
a new internal response channel.

## 2. A finite exact recognition test

Suppose proposed shell factors `D_t` are fixed and write

```math
E_t=K_t-D_t.                                     \tag{DU.9}
```

For one factor define its rectangular defect

```math
\square E_t(a,a';b,b')
=E_t(a,b)+E_t(a',b')-E_t(a,b')-E_t(a',b).        \tag{DU.10}
```

Across one internal interface define

```math
\begin{aligned}
\partial(E_t,E_{t+1})(a;b,b';c)
={}&E_t(a,b)-E_t(a,b')\\
 &+E_{t+1}(b,c)-E_{t+1}(b',c).                  \tag{DU.11}
\end{aligned}
```

### Proposition DU.2 (rectangle plus interface criterion)

There exist functions `phi_0,...,phi_T` and scalars `c_t` satisfying
(DU.4) if and only if

```math
\square E_t=0\quad(1\le t\le T),
\qquad
\partial(E_t,E_{t+1})=0\quad(1\le t<T).         \tag{DU.12}
```

This is a finite certificate involving individual factors and adjacent
interfaces, not all paths of the composed system.

#### Proof

Equation (DU.4) makes both alternating sums vanish.  Conversely,
`square E_t=0` implies

```math
E_t(a,b)=r_t(a)+s_t(b)                            \tag{DU.13}
```

after fixing one base row and column.  The second condition says
`s_t(b)+r_(t+1)(b)` is constant in `b`.  Take `phi_0=r_1`,
`phi_t=-s_t` for `1<=t<T`, and `phi_T=-s_T`; absorb the resulting constants
into the `c_t`.  This gives (DU.4). `square`

The conditions are exact for a reason.  A uniformly small nonzero
rectangular or interface circulation may be repeated and need not remain
small at long depth.

This test uses complete real-valued interfaces.  On a sparse bipartite
support, rectangles need not generate the support cycle space; additive
separability instead requires zero alternating sum on every support cycle.
For a branching word language, compatibility must also hold on every
allowed adjacency in its actual finite context graph.

## 3. The repeatable-cycle rigidity behind the criterion

Let `G=(V,E)` be a finite strongly connected directed graph and let every
edge carry a label `zeta_e` in a normed real vector space `Z`.  For a walk
`P`, write `zeta(P)` for the sum of its edge labels.

### Proposition DU.3 (bounded-depth stability is exact cohomology)

The following are equivalent.

1. The labels of all finite directed walks are uniformly bounded.
2. Every directed closed walk has label zero.
3. There is a vertex potential `p:V->Z` such that

   ```math
   \zeta_{u\to v}=p(v)-p(u).                     \tag{DU.14}
   ```

4. Any two directed walks with the same endpoints have equal labels.

If a closed walk `C` has `zeta(C) ne 0`, then its `k`-fold repetition has
label `k zeta(C)`.  Hence no nonzero tolerance on a repeatable cycle can by
itself imply a path-length-uniform error bound.  The same statement holds in
a response quotient such as `R^Y/R1`: scalar baselines are first set to
zero in the quotient.

#### Proof

Boundedness applied to repetitions of a closed walk proves (2).  Fix a root.
Under (2), define `p(v)` by the label of any root-to-`v` walk; strong
connectivity and a return walk show that this is path independent.  This
proves (3), and telescoping proves (4) and (2).  Conversely, delete
zero-label closed subwalks until a walk is simple; it then has at most
`|V|-1` edges, proving uniform boundedness.  Repetition gives the last
assertion. `square`

On a general finite directed graph the statement applies separately to each
strongly connected component.  An acyclic condensation can be crossed only
finitely many times, so arbitrary-depth amplification can arise only from a
repeatable component.

For a declared regular language, apply the criterion to the reachable and
co-reachable part of the finite automaton (or product graph) encoding that
language.  Only cycles that are pumpable there can force drift; a forbidden
or unreachable raw cycle is irrelevant.

This is ordinary graph cohomology.  Its role here is diagnostic: a proposed
uniform approximate-lumpability theorem must either repair recurrent
holonomy **exactly**, restrict repetition depth, or include a genuine memory
reset.  Merely bounding every local cycle by `delta` cannot work, because a
nonzero cycle can be traversed `Theta(1/delta)` times.

## 4. Why strict max-plus contraction is not the missing mechanism

For a finite real `r by s` matrix `S`, let

```math
(F_Su)_b=\max_a\{u_a+S_{ab}\}.                   \tag{DU.15}
```

Use Hilbert's projective metric

```math
d_H([u],[v])={1\over2}\operatorname {osc}(u-v).  \tag{DU.16}
```

Define the projective image diameter

```math
\Delta(S)={1\over2}\max_{b,c}
\left[
 \max_a(S_{ab}-S_{ac})-
 \min_a(S_{ab}-S_{ac})
\right].                                         \tag{DU.17}
```

### Theorem DU.4 (zero-or-one tropical contraction)

On the full projective domain `R^r/R1`, the map `F_S` is nonexpansive in
`d_H`, its projective image has diameter exactly `Delta(S)`, and its global
projective Lipschitz coefficient is

```math
\operatorname {Lip}_H(F_S)=
\begin{cases}
0,&\Delta(S)=0,\\
1,&\Delta(S)>0.
\end{cases}                                      \tag{DU.18}
```

Moreover `Delta(S)=0` exactly when

```math
S_{ab}=\alpha_a+\beta_b,                          \tag{DU.19}
```

in which case `F_Su` is projectively equal to `beta` for every input.  Thus
there is no nontrivial global contraction coefficient strictly between zero
and one for a finite max-plus linear map.  A depth-uniform theorem cannot be
obtained by postulating an ordinary contraction factor unless it has already
postulated a complete projective reset.

#### Proof

Nonexpansiveness follows from

```math
\min_a(u_a-v_a)
\le(F_Su)_b-(F_Sv)_b
\le\max_a(u_a-v_a).                               \tag{DU.20}
```

For every pair of output coordinates, their difference lies between the
minimum and maximum row difference in (DU.17).  Conversely, sending one
input coordinate far above the others makes the corresponding row dominate
all outputs.  Two such limits attain every extremal cross difference, which
proves the diameter formula.

On each full-dimensional cell of the comparison arrangement, every output
has a fixed unique maximizing row.  The derivative is therefore a coordinate
selection map.  Its projective operator norm is zero when all outputs select
the same input coordinate, and one otherwise.  If the projective map is
nonconstant, some full-dimensional cell has nonzero projective derivative:
a continuous piecewise-affine map with zero projective derivative on every
full-dimensional cell is constant.  Hence its global coefficient is one.
Finally `Delta(S)=0` says all row differences are constant across columns,
which is equivalent to (DU.19). `square`

The full-domain qualification is essential.  On a restricted reachable
subset, a nonconstant map can have an intermediate coefficient.  For
example, in one projective coordinate the all-finite matrix

```math
\begin{pmatrix}0&0\\-\delta&0\end{pmatrix}
```

acts as `z -> min(max(z,0),delta)`.  Its global coefficient is one, but on
the two-point set `{-delta,2delta}` its coefficient is `1/3`.

## 5. Small-diameter resets do give uniform stability

The preceding dichotomy does not make small `Delta(S)` useless.  It controls
absolute projective memory after the map, rather than multiplying the
previous error.

Let `F_ell:X->X` be raw updates, `G_ell:Y->Y` nonexpansive updates on a metric
space, and `pi:X->Y` an approximate state map satisfying

```math
d\bigl(\pi F_\ell(x),G_\ell\pi(x)\bigr)\le\epsilon
\quad\hbox{for every }x,\ell.                    \tag{DU.21}
```

Call a word `w` a `rho`-reset when

```math
\operatorname {diam}G_w(Y)\le\rho.              \tag{DU.22}
```

### Theorem DU.5 (recurrent-reset semiconjugacy)

Run the raw and quotient systems from matching initial states.  Suppose that
along a declared letter sequence, every sufficiently late time `t` has a
`rho`-reset block of length at most `L` which ends no more than `L` steps
before `t`.  The diameter is taken for the quotient map on a set containing
both quotient trajectories.  Then, after the first such block,

```math
\boxed{d\bigl(\pi(x_t),y_t\bigr)
\le\rho+2L\epsilon}                              \tag{DU.23}
```

at every time.  The right-hand side is independent of the total depth.
For nonuniform local errors, replace the two `L epsilon` terms by the error
sums inside the chosen reset block and after it.

#### Proof

Let the chosen reset occupy steps `s+1,...,j`.  Iterating (DU.21) from the
actual encoded state at time `s` gives

```math
d\bigl(\pi(x_j),G_w\pi(x_s)\bigr)\le L\epsilon.
```

The quotient trajectory satisfies `y_j=G_w(y_s)`, and the two images of the
reset word are within `rho`, regardless of the error at time `s`.  The tail
from `j` to `t` has length at most `L`; nonexpansiveness and (DU.21) add at
most another `L epsilon`. `square`

Equivalently, for a fixed reset word of length `m` whose last completed
occurrence is at most `G` letters behind, the exact bound is
`rho+(m+G)epsilon`.  A uniform statement over all prefixes must also bound
the first completion time (or handle the finite initial segment directly),
and the local defect must hold on the prefix closure of the declared word
family.

For max-plus quotient maps, (DU.17) computes `rho` exactly for a one-letter
reset.  In a binary zero-temperature Ising transfer with

```math
S_J(s,t)=Jst,
```

one has

```math
\Delta(S_J)=2|J|.                                \tag{DU.24}
```

Thus every nonzero bond still has global projective Lipschitz coefficient
one, but a weak bond erases prior memory down to diameter `2|J|`.  If every
gap of length `L` contains a bond with `|J|<=eta`, an `epsilon`-accurate
quotient has uniform error at most `2eta+2L epsilon` after the first weak
bond.  This is a quantitative distinction between contraction and reset
that the usual one-step nonexpansive estimate does not see.

## 6. What the combined theorem says

There are now two finite, genuinely depth-uniform mechanisms.

1. **Cohomological cancellation:** local errors are endpoint gauges, and
   internal channels cancel exactly before optimization.
2. **Memory renewal:** a recurrent continuation maps the entire old state
   into a small projective set, so only errors since the last reset matter.

They also supply a sharp obstruction.  Max-plus dynamics do not offer a
generic strict projective contraction, and a nonzero repeatable holonomy
amplifies linearly.  The transition-toll kernel

```math
K_\delta(i,j)=a|i-j|-\delta\boldsymbol1_{i\ne j}
```

has rectangular defect `2delta` already on two adjacent indices, so it is
not an endpoint gauge; its repeated drift is exactly the cycle/toll side of
the dichotomy.

This is stronger than saying that a quotient is nonexpansive.  It gives a
finite exact certificate for cancellation, an explicit no-go for generic
strict contraction, and a quantitative sufficient reset condition.  It is
still a theorem about structured tropical continuations, not a universal
repair theorem for arbitrary approximate lumpings.
