# Wave 20 memo: puncture cycles are assignment reduced costs, not landing certificates

## Status

All identities and inequalities below are **proved exactly** in the one-copy
normalization $`M=Q/2`$.  The finite claims are exhaustive integer
computations, checked by `tmp/puncture_cycles_r20.py`.  Zero directed cycles
occur in the exact minimizers $`A_8`$ and $`A_9`$ while the restriction on the
cycle vertices has positive terminal excess.  Thus cycle totals do not give a
deterministic landing theorem.  Nothing here proves convergence.

## The cross-puncture matrix is an assignment reduced-cost matrix

Retain the notation of (10.653)--(10.659).  Thus $`A`$ is an exact order-$`n`$
minimizer, $`M=M_n`$, $`q_i`$ is an optimally extended oriented ground of
$`A[-i]`$, and

```math
\langle A,q_i\rangle=M-s_i,
\qquad
r_i(q_i)=d_i-s_i,
\qquad
u_{ij}=s_i+r_j(q_i)-d_j\ge0.
```

Define the incident-field assignment matrix

```math
c_{ij}=r_j(q_i).
```

Then

```math
\boxed{
u_{ij}=c_{ij}-(-s_i)-d_j,
\qquad
c_{ii}=d_i-s_i.
}
\tag{R20.1}
```

Consequently $`\alpha_i=-s_i`$ and $`\beta_j=d_j`$ are feasible dual
potentials for the minimum-cost assignment problem with costs $`c_{ij}`$:

```math
c_{ij}\ge\alpha_i+\beta_j.
```

The diagonal assignment attains equality in every row, and its cost equals
the dual value.  It is therefore an optimal assignment, and $`u_{ij}`$ is
exactly its matrix of reduced costs.

Equivalently, for every permutation $`\pi`$,

```math
\boxed{
\sum_i u_{i,\pi(i)}
=\sum_i\left[r_{\pi(i)}(q_i)-r_i(q_i)\right]
\ge0.
}
\tag{R20.2}
```

The same identity holds after replacing a permutation matrix by any doubly
stochastic matrix.  If $`C=(i_1,\ldots,i_k)`$ is a directed cycle and all
other coordinates are fixed, (R20.2) becomes

```math
\boxed{
U(C):=\sum_{\ell=1}^k u_{i_\ell,i_{\ell+1}}
=\sum_{\ell=1}^k
\left[r_{i_{\ell+1}}(q_{i_\ell})-r_{i_\ell}(q_{i_\ell})\right],
}
\tag{R20.3}
```

where $`i_{k+1}=i_1`$.  Thus a cycle total is precisely the cost increase
from cyclically reassigning the distinguished deleted coordinates.  There
can be no negative augmenting cycle.  A zero cycle merely gives another
optimal field assignment.  This conclusion uses only (10.654); it does not
yet use that the words are cuts.

## Cut support has trivial holonomy

The cut structure gives an exact interpretation of each edge, but no new
cycle holonomy.  Gauge the child $`B_j=A[-j]`$ by its oriented ground
$`q_j[-j]`$ and put

```math
w_e=a_eq_{j,e},
\qquad
\sum_{e\in E(B_j)}w_e=M(B_j).
```

The ratio $`q_iq_j`$ is another augmented-cut word.  Write it as
$`\epsilon y_uy_v`$, with $`\epsilon\in\{\pm1\}`$, and let $`S`$ be the
negative shore of $`y`$ inside $`V\setminus\{j\}`$.  Directly from the
definition of $`u`$,

```math
\boxed{
\frac{u_{ij}}2
=
\begin{cases}
\displaystyle\sum_{e\in\delta(S)}w_e,
&\epsilon=+1,\\[3mm]
\displaystyle\sum_{e\notin\delta(S)}w_e
=M(B_j)-\sum_{e\in\delta(S)}w_e,
&\epsilon=-1.
\end{cases}
}
\tag{R20.4}
```

Hence $`u_{ij}=0`$ says exactly that $`q_i[-j]`$ is another oriented ground
of $`B_j`$: it is a zero signed cut in the same-orientation case and a
zero complement residual in the opposite-orientation case.

For every directed cycle,

```math
\boxed{
\prod_{\ell=1}^k(q_{i_\ell}q_{i_{\ell+1}})=1
}
\tag{R20.5}
```

coordinatewise.  In particular, the number of complemented transitions is
even and the symmetric difference of all transition supports is empty.
This is an identity, not a rigidity theorem: the $`q_i`$ are already global
full cut words, so their connection is pure gauge and has identically trivial
holonomy.

There is also an exact boundary localization.  Since

```math
u_{ij}=\langle A[-j],q_j[-j]-q_i[-j]\rangle,
```

the full-matrix terms telescope around a cycle and give

```math
\boxed{
\begin{aligned}
U(C)
&=-\sum_{\ell=1}^k
  \sum_{v\ne i_{\ell+1}}
  a_{i_{\ell+1}v}
  \left(q_{i_{\ell+1},i_{\ell+1}v}
       -q_{i_\ell,i_{\ell+1}v}\right)\\
&=-\sum_{e=\{a,b\}}a_e
  \sum_{\ell:i_{\ell+1}\in\{a,b\}}
  \left(q_{i_{\ell+1},e}-q_{i_\ell,e}\right).
\end{aligned}
}
\tag{R20.6}
```

All edges disjoint from the cycle vertices cancel.  Since every word
difference is in $`\{0,\pm2\}`$,

```math
\boxed{0\le U(C)\le2k(n-1).}
\tag{R20.7}
```

For a two-cycle, common-edge cancellation recovers the sharper bound
$`u_{ij}+u_{ji}\le4(n-2)`$ from (10.658).  Equations (R20.3)--(R20.7) show
that cycle totals remain boundary cocycles.  The genuinely cut-specific
consequence of a zero cycle is only a cyclic chain of alternative exact
child grounds.

## The first missing terminal data

Cycle reassignment concerns child states, while optimized restriction
concerns the norm of a common terminal matrix.  The separation is exact
already after two deletions.  Put

```math
T=V\setminus\{i,j\},
\qquad
v_{i;ij}
=M(A[T])-\langle A[T],q_i[T]\rangle\ge0.
\tag{R20.8}
```

Here $`v_{i;ij}`$ is the **terminal refresh deficit**: it measures the gain
from allowing a new maximizing cut after both vertices have been deleted.
Removing two incident stars from $`q_i`$ and adding their shared edge back
once gives

```math
\begin{aligned}
\langle A[T],q_i[T]\rangle
&=\langle A,q_i\rangle-r_i(q_i)-r_j(q_i)
  +a_{ij}q_{i,ij}\\
&=M-d_i-d_j+s_i-u_{ij}+a_{ij}q_{i,ij}.
\end{aligned}
```

Therefore

```math
\boxed{
M(A[T])
=M-d_i-d_j+s_i-u_{ij}+a_{ij}q_{i,ij}+v_{i;ij}.
}
\tag{R20.9}
```

If

```math
\varepsilon_{ij}=M(A[T])-M_{n-2},
```

then $`d_i=\delta_n-e_i`$ yields the terminal-excess form

```math
\boxed{
\varepsilon_{ij}
=e_i+e_j+(\delta_{n-1}-\delta_n)
  +s_i-u_{ij}+a_{ij}q_{i,ij}+v_{i;ij}.
}
\tag{R20.10}
```

This identity explains the sign of the cycle obstruction.  A **large**
$`u_{ij}`$, not a small one, lowers the carried terminal baseline; the fresh
deficit $`v_{i;ij}`$ can then refund that gain.  A zero cycle gives no
negative $`u`$ credit at all.  Moreover the inherited defects, scalar
curvature $`\delta_{n-1}-\delta_n`$, and the internal edge remain.

For a general deleted set $`S`$ and an anchor word $`q_i`$, define

```math
R_i(S)=\sum_{e:e\cap S\ne\varnothing}a_eq_{i,e},
\qquad
v_i(S)=M(A[V\setminus S])
       -\langle A[V\setminus S],q_i[V\setminus S]\rangle.
```

Then the exact terminal decomposition is

```math
\boxed{
M(A[V\setminus S])-M_{n-|S|}
=M-M_{n-|S|}-s_i-R_i(S)+v_i(S).
}
\tag{R20.11}
```

The baseline-removal/refresh pair $`(R_i(S),v_i(S))`$ is the first genuinely
terminal, non-cyclic datum.  Cycle totals determine neither component.
Even $`v_i(S)`$ alone is insufficient: a carried word can already be a
terminal ground while the terminal matrix remains globally nonoptimal.

## Exact finite falsifiers

The verifier includes both choices of the restored spin when its field is
zero and exhausts every projective child ground.

1. The displayed $`A_6`$ in the verifier has $`M(A_6)=M_6=5`$ and all six
   children have norm $`M_5=4`$.  It has a zero two-cycle on vertices
   $`(0,1)`$ for which both selected anchors satisfy

   ```math
   u_{01}=u_{10}=0,
   \qquad
   v_{0;01}=v_{1;01}=2.
   ```

   Thus even under exact global minimality, zero cycle cost does not control
   terminal refresh.

2. For the exact minimizer $`A_8`$ of (10.445), a zero two-cycle on
   $`(0,3)`$ lands on an order-six restriction with

   ```math
   M(A_8[-\{0,3\}])=7=M_6+2.
   ```

   A zero three-cycle on $`(0,3,4)`$ lands on an order-five restriction with

   ```math
   M(A_8[-\{0,3,4\}])=6=M_5+2.
   ```

   The three anchor refresh deficits can be $`(4,0,0)`$.

3. For the exact minimizer $`A_9`$ of (10.298), a zero two-cycle on
   $`(0,2)`$ lands on an order-seven restriction with

   ```math
   M(A_9[-\{0,2\}])=11=M_7+2.
   ```

   A zero three-cycle on $`(0,2,3)`$ lands on an order-six restriction with

   ```math
   M(A_9[-\{0,2,3\}])=9=M_6+4.
   ```

   All three selected anchor refresh deficits in the latter cycle are zero:
   the carried words are already terminal grounds, but the terminal matrix
   itself has excess four.

These examples falsify every universal implication of the form

```math
\text{terminal excess on the cycle vertices}\le f(U(C)),
\qquad f(0)=0,
```

and every deterministic landing rule whose only nonlinear certificate is
zero cut holonomy or zero cycle cost.  They are finite obstructions and do
not rule out a tail-uniform theorem with an additive constant or a new
macroscopic invariant.

## Disposition

Directed cycles do not bypass the linear cocycle from Wave 19.  They are the
reduced-cost cycles of an assignment problem whose diagonal solution is
already primal-dual optimal, while cut-ratio holonomy is identically trivial.
Low and zero cycles certify alternative child grounds but provide neither a
simultaneous multivertex exchange nor a terminal-excess drop.

The first viable continuation must couple the terminal removal/refresh data
in (R20.11) to exact global minimality.  A concrete useful target would be a
macroscopic assignment or subset for which the favorable $`u`$ and removed
internal energy cannot be refunded by $`v_i(S)`$, with a power-saving bound
on the resulting terminal excess.  Without such a refresh bound, cycle
selection is not a leading-route engine.

