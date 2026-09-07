# Asymmetric incidence capacity: an exact seed threshold at a subleading parent scale

Date: 2026-09-07. Status: proved finite combinatorial statement and scale audit;
**not** a leading-order seed-transfer theorem or a convergence argument.

Write `H_S(x)=sum_{i<j} S_ij x_i x_j` and `Q(S)=max_x |H_S(x)|`.
An incidence assignment has signs `u_ij` for ordered pairs `i != j`.
At each vertex permit at most `r` negative incidences. A mismatch is an edge
with `u_ij u_ji != S_ij`. Unlike balanced symmetric incidence types, this
one-sided capacity can retain a switching-dependent seed energy.

## 1. Exact orientation formula

Let `G` be the graph of negative edges of `S`, with `e(G)=e`. The smallest
possible number of mismatches under the capacity constraint is

```
D_r(G) = max_{U subset V} (e_G(U)-r |U|),
```

where the empty set makes the maximum nonnegative.

Proof. Every correctly realized negative edge needs exactly one negative
incidence, so the correctly realized negative edges form a graph whose edges
can be assigned to endpoints with capacity `r` per endpoint. Conversely, any
such edge assignment is realized by giving its assigned endpoint the negative
incidence and all other incidences the positive sign. All positive edges are
then correct, and precisely the unassigned negative edges are mismatches.

The endpoint assignment is the elementary integral-flow problem with a
unit-capacity source edge for each edge of `G`, infinite-capacity arcs from
that edge-vertex to its two original endpoints, and endpoint-to-sink capacity
`r`. A finite cut determined by a set `U` of original vertices puts at most
the `e_G(U)` corresponding edge-vertices on the source side; its minimum
capacity is `e-e_G(U)+r|U|`. Integral max-flow/min-cut gives the formula.

In particular, if the full graph is a densest induced subgraph, meaning

```
e_G(U)/|U| <= e/m                 (nonempty U),
```

then, for every nonnegative integer `r`,

```
D_r(G) = max(0,e-rm).
```

Indeed if `r >= e/m` all deficiencies are nonpositive. If `r < e/m`, every
deficiency is at most `(e/m-r)|U| <= e-rm`, attained at `U=V`.

## 2. Densest-full-graph criterion in a low-cap regular gauge

Fix `C`. For all sufficiently large `m`, if a signing `S` satisfies

```
Q(S) <= C m^(3/2),       min degree(G) >= .45 m,
```

then its negative graph `G` has the full graph as a densest induced subgraph.

To check this, set `a=e/m`. The cap bound gives
`a=(m-1)/4+O_C(sqrt(m))`, hence `.24m <= a <= .26m` for large `m`.
If `1 <= u=|U| <= m/8`, then `e_G(U)/u <= (u-1)/2 <= m/16 < a`.
If `m/8 <= u <= 3m/4`, principal restriction and the cap bound give

```
e_G(U)/u = (u-1)/4 - H_{S[U]}(1)/(2u)
          <= 3m/16 + 4C sqrt(m) < a.
```

Finally let `u >= 3m/4` and `s=m-u`. The number of edges incident with the
complement is at least

```
.45 m s - binom(s,2) >= (.45-.125) m s >= a s.
```

Consequently `e_G(U) <= e-as=au`, as required. Empty complements cause no
problem.

The minimum-degree hypothesis is available, uniformly over **every**
switching of a bounded-operator low-cap signing, after deleting only a
bounded number of vertices. Specifically, if `||S||op <= K sqrt(m)`, its row
fields `f=S 1` satisfy `sum_i f_i^2 <= K^2 m^2`. There are at most `400 K^2`
vertices with `|f_i| > m/20`. Delete those vertices. At each remaining vertex,
the old negative degree is at least `(m-1-m/20)/2`; deleting the bounded
exceptional set still leaves degree at least `.45` times the new order for
large `m`. Principal cap and operator bounds persist, with harmless constant
changes. The same argument applies after any diagonal sign switching because
operator norm and cap are switching invariant. Deletion changes the energy
of any fixed signing state by `O_K(m)`, not by order `m^(3/2)`.

This applies to the selectable fixed-accuracy, bounded-operator near-minimizer
normal form in `principle_synthesis_2026_09_07_global_balancing.md`. It is not
being asserted that all exact minimizers have a fixed uniform operator bound.

## 3. Actual seed sensitivity and its exact scale

For a regularized gauge as above, put

```
r = (m-1)/4 - alpha sqrt(m)/2 + O(1).
```

Because `e=(binom(m,2)-H_S(1))/2`, the orientation formula becomes

```
D_r(G) = (alpha m^(3/2)-H_S(1))_+/2 + O(m).
```

Thus zero-defect feasibility detects `H_S(1) >= alpha m^(3/2)+O(m)`.
Allowing a selectable switching genuinely exposes the seed optimization;
there is no seed-blind universality claim for this asymmetric capacity class.
No entropy count is needed for this statement.

However, consider the usual rank-one-weave coordinate normalization
`h_i(j)=sqrt(k) u_ij`. Its off-diagonal unordered energy is

```
sum_{i<j} S_ij h_i(j)h_j(i) = k [binom(m,2)-2D_r(G)].
```

A seed-energy change of order `m^(3/2)` therefore changes this parent energy
by only `O(k m^(3/2))`. At parent order `N=mk`, the normalized change is

```
O(k m^(3/2)/N^(3/2)) = O(k^(-1/2)).
```

In the proportional weave regime `k ~ m`, the physical change is
`O(N^(5/4))`, strictly below the target `N^(3/2)` scale. Increasing a
Laplace parameter to order `sqrt(m)` can make the effect visible in an
`m^2`-scale logarithmic moment, but the Markov threshold divides by that same
parameter. It does not change the physical suppression. In the standard
weave cap formula `(t+R_t)/(2t sqrt(p))`, a bounded normalized row entropy
term divided by `t -> infinity` leaves the leading baseline
`1/(2 sqrt(p)) >= 1/2`.

Accordingly this result separates asymmetric seed-sensitive compilers from
balanced seed-blind compilers, but it does **not** transfer an optimized seed
constant into a leading-order subhalf parent constant. A useful seed landing
mechanism must overcome the explicit `k^(-1/2)` loss, not merely produce a
seed-sensitive hard type or an `m^2`-scale tilted pressure.
