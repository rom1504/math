# Independent abstraction audit: code amalgamation and low-variance bridges

Date: 2026-07-31.  This note is an agent-authored research report.  It does
not modify the project's user directives or strategic files.

## 1. Exact augmented-code normalization

For a graph `G`, encode a signing by a word `a` in `F_2^{E(G)}` and put

```math
D_G=C(G)+\langle\mathbf 1\rangle,
```

where `C(G)` is the binary cut code.  If `N_G=|E(G)|`, then

```math
\operatorname{cap}_G(a)=N_G-2d(a,D_G),\qquad
M(G)=N_G-2\rho(D_G).                                      \tag{A.1}
```

Indeed, for a vertex word `z`, the signed energy is agreements minus
disagreements with the cut `delta z`; adjoining `1` takes the absolute
value.  Thus the displayed factors and directions agree with the project's
Hamiltonian convention `H_A(x)=sum_(i<j) a_ij x_i x_j`.

Let

```math
B_{m,n}=\min_{W\in\{\pm1\}^{m\mathbin\times n}}
       \max_{x\in\{\pm1\}^m,y\in\{\pm1\}^n}|x^{\mathsf T}Wy|. \tag{A.2}
```

The augmented cut code of `K_(m,n)` is used in (A.2).  Global negation of a
bipartite energy is obtained by replacing `y` by `-y`, so the usual
one-sided Gale--Berlekamp objective and the absolute objective in (A.2) are
identical.

### Proposition A (independent deep-hole amalgamation)

For all positive `m,n`,

```math
M_{m+n}\le M_m+M_n+B_{m,n}.                              \tag{A.3}
```

To prove this, take deepest words `a,b,w` for `D_(K_m)`, `D_(K_n)`, and
`D_(K_(m,n))`.  In the coordinate decomposition

```math
E(K_{m+n})=E(K_m)\sqcup E(K_n)\sqcup E(K_{m,n}),
```

every parent cut restricts to one codeword in each of the three augmented
codes.  Consequently its distance from `(a,b,w)` is at least the sum of the
three covering radii.  Hence

```math
\rho(D_{K_{m+n}})\ge
\rho(D_{K_m})+\rho(D_{K_n})+\rho(D_{K_{m,n}}),            \tag{A.4}
```

and substituting (A.1) proves (A.3).

This is a valid composition inequality, but it is not a convergence
mechanism.  The rectangular term is of leading order for comparable shores.
In particular, Bowlin's theorem below implies

```math
B_{n,n}>n^2 2^{-(n-1)}
 \binom{n-1}{\lfloor(n-1)/2\rfloor}
 =(\sqrt{2/\pi}+o(1))n^{3/2}.                            \tag{A.5}
```

Thus any proof that replaces the state-dependent bridge budget by the
standalone quantity `B_(m,n)` pays a leading, not summable, cost.

## 2. Verified primary theorem and exact response-rank obstruction

Bowlin, [*Maximum Frustration in Bipartite Signed
Graphs*](https://doi.org/10.37236/2204), Theorems 27 and 25--26, proves the
following.  With

```math
c_l=2^{-(l-1)}\binom{l-1}{\lfloor(l-1)/2\rfloor},
```

the maximum one-sided frustration of `K_(l,r)` obeys

```math
F_{\max}(K_{l,r})\le {lr\over2}(1-c_l),                  \tag{A.6}
```

with equality if and only if `r` is a positive multiple of `2^(l-1)`.
The exact mapping `B_(l,r)=lr-2F_max(K_(l,r))` therefore gives

```math
B_{l,r}\ge lr c_l,                                      \tag{A.7}
```

with equality under precisely that divisibility condition.  No absolute
versus one-sided mismatch remains because one shore can be globally flipped.

The proof gives a further structural obstruction.  A column of an `l x r`
bridge has one of `2^(l-1)` types modulo negation.  If `q` is the vector of
type counts, the complete optimized cut-response profile is a linear map
`K_l q`.  Bowlin proves that the associated matrix obtained by adjoining the
total-count row is invertible.  Therefore

```math
q\longmapsto(r,K_lq)                                    \tag{A.8}
```

is injective and has rank `2^(l-1)`.

This proves a scalable no-compression result, with a limited but exact
scope: no lower-dimensional **exact linear** quotient of the column-type
counts retains every optimized switching response.  At comparable shores,
the exact response state is exponential.  It does not rule out nonlinear or
approximate compression, but any such proposal must explicitly use an
approximation scale and cannot appeal to an unproved finite-state collapse.

## 3. A new sufficient mechanism: low-variance fractional bridge rounding

For fixed child signings `A,B`, the exact block identity is

```math
\operatorname{cap}\begin{pmatrix}A&C\\C^{\mathsf T}&B\end{pmatrix}
=\max_{x,y}\bigl(|H_A(x)+H_B(y)|+|x^{\mathsf T}Cy|\bigr). \tag{A.9}
```

For a proposed energy target `T`, define the fractional bridge body

```math
\mathcal P_T(A,B)=\left\{C_0\in[-1,1]^{m\times n}:
 |H_A(x)+H_B(y)|+|x^{\mathsf T}C_0y|\le T
 \text{ for all }x,y\right\}.                           \tag{A.10}
```

If `p=cap(A)` and `q=cap(B)`, the zero-defect energy target

```math
T_0=(p^{2/3}+q^{2/3})^{3/2}                            \tag{A.10a}
```

satisfies `T_0>=p+q`.  Hence `C_0=0` always belongs to
`P_(T_0)(A,B)`.  Fractional feasibility is therefore not the issue.  The
geometric question is whether this slab body penetrates close enough to the
vertices of the bridge cube.

### Proposition B (uniform low-variance rounding)

Suppose `C_0` belongs to (A.10), and set

```math
V(C_0)=\sum_{i,j}(1-(C_0)_{ij}^2),\qquad
L=(m+n+2)\log2.                                         \tag{A.11}
```

Then there is a sign bridge `C` for which

```math
\operatorname{cap}\begin{pmatrix}A&C\\C^{\mathsf T}&B\end{pmatrix}
\le T+\sqrt{2V(C_0)L}+{4\over3}L.                       \tag{A.12}
```

For the proof, independently round each entry to `+1` or `-1` with mean
`(C_0)_ij`.  For fixed `x,y`, the centered error

```math
Z_{x,y}=x^{\mathsf T}(C-C_0)y
```

has total variance `V(C_0)` and summands bounded in absolute value by `2`.
Bernstein's inequality and the standard inverted bound give

```math
\Pr\{|Z_{x,y}|>\sqrt{2VL}+4L/3\}\le2e^{-L}.             \tag{A.13}
```

There are at most `2^(m+n)` ordered spin pairs, so the union probability is
strictly below one with the safe choice (A.11).  Equation (A.9) then proves
(A.12).

For comparable blocks of total order `N`, a uniform construction satisfying

```math
V(C_0)=O(N^{2-\eta})                                    \tag{A.14}
```

for some fixed `eta>0` has rounding loss
`O(N^(3/2-eta/2)+N)`.  Without using a derivative estimate,

```math
(T+D)^{2/3}\le T^{2/3}+D^{2/3}                          \tag{A.15}
```

shows a power-saving `b`-defect `O(N^(1-eta/3)+N^(2/3))`, which is
geometrically summable in the campaign's criterion.

This is a genuinely weaker interface than integral bridge optimization **if
a fractional bridge is supplied independently**: the two scalar quantities
`T` and `V(C_0)` provably control the full integral parent cap.  The exact
missing lemma is now:

> Construct, for the chosen child family and near the ideal compositional
> target `T`, a member of `P_T(A,B)` with (A.14), without defining the
> children or `C_0` through the unknown optimum.

A stronger uniform form, which would compose even unstructured children, is
the following **cube-boundary penetration lemma**: for every comparable pair
with `p,q=O(N^(3/2))`, the body `P_(T_0+O(N^(3/2-delta)))(A,B)` contains a
point with `V=O(N^(2-eta))`, for fixed positive `delta,eta`.  Proposition B
would turn this directly into a power-saving `b`-recurrence.  This statement
uses only the supplied child caps, not `M_m` or `M_n`, and is therefore
noncircular.  It is also falsifiable by a family whose fractional body stays
a positive normalized variance away from every cube vertex.

The exact falsifier is a lower bound `V(C_0)=Omega(N^2)` for every
near-ideal `C_0` in (A.10) along a scalable family.  Both claims are amenable
to a new computation: optimize or bound `sum C_ij^2` over the fractional
bridge body, rather than adding further uncompressed binary separator states.
The optimization is nonconvex (maximizing a convex quadratic over a
polytope), so an LP optimum alone is not a certificate of minimum variance.
Moreover, (A.10) has `2^(m+n)` state inequalities before projective
identifications.  Its raw membership/separation problem is the same full
parent maximization that the campaign is trying to avoid.  Consequently the
definition of `P_T` is not itself a bounded-complexity state or a reduction:
progress requires an explicit formula for `C_0`, or a compressed certificate
of (A.10), in addition to the variance estimate.

An equivalent useful sufficient state is a feasible fractional bridge with
only `r=O(N^(2-eta))` nonintegral entries, since then `V(C_0)<=r`.  At a
polytope vertex, the number of fractional box coordinates is at most the rank
of its active non-box state constraints.  Thus a theorem producing a
near-ideal vertex with a power-saving active rank would establish (A.14).
Bowlin's full-rank theorem warns that this cannot hold for the complete exact
standalone response profile without using the children and their nonuniform
margins.

There is an important route-identity warning.  Obtaining (A.14) *only* by
asserting a low-rank active face reintroduces the earlier common-active-face
obligation: one must still control exponentially rare tight states and their
rank, while the established theorem controls only polynomial-scale moments.
Bowlin's rank `2^(l-1)` result makes the corresponding exact standalone claim
false.  Therefore the active-rank corollary of Proposition B is not by itself
a genuinely new ingredient.  The new content is only the variance-sensitive
rounding implication; it remains independent of the inactive route if an
analytic near-integral `C_0` and a direct proof of all its margins are found
without a common-face entropy assertion.

## 4. Audit of the current structured-landing formulation

Let `F_n` be any family defined independently of `M_n` and put
`u_n=min_(S in F_n) cap(S)^(2/3)`, `b_n=M_n^(2/3)`.  If `u_n/n` is already
known to converge, the statement

```math
0\le u_n-b_n=o(n)                                      \tag{A.16}
```

immediately proves the desired convergence.  Therefore (A.16), by itself,
is an equivalent asymptotic obligation rather than an explanation.  It is
demonstrably simpler only when accompanied by one of the following
independent mechanisms:

1. a uniform transformation of arbitrary near-minimizers into `F_n` with
   `o(n^(3/2))` cap loss;
2. an all-signings lower bound matching the explicit `F_n` upper bound; or
3. a bridge interface such as Proposition B whose state and rounding theorem
   remove the exponential integral optimization.

Pure bounded local/restriction-profile compactness is insufficient.  The
project's previously verified planted-block construction changes only
`s=n^(3/4)` vertices, and hence changes every fixed-`k` empirical principal
restriction distribution by at most `ks/n=o(1)`, while its `s^2=Theta(n^(3/2))`
edges change the normalized Boolean cap by a constant.  Thus a viable state
must include a uniform-integrability or spike-control coordinate; fixed local
statistics alone cannot establish (A.16).

The independently checked contemporaneous switching-template theorem gives
a complementary warning for bounded operational states: any template using
at most `K` switching types has cap at least `n^2/(8K)`, so `K=o(sqrt(n))`
has a superlinear landing gap.  Its proof combines a largest-cell principal
restriction with orthogonality of the quotient's degree-two characters.
Thus a fractional construction cannot claim simplicity merely because it is
generated by a bounded switching script.  A viable `C_0` must either have
effective complexity at least the critical `sqrt(n)` scale or exploit a
different algebraic certificate whose margins are checked collectively.

## 5. Other blank-slate candidates after comparison

- A scalar energy-histogram/Hoeffding bridge is exactly the already recorded
  random-bridge criterion and is far from finite optimized bridges.  Generic
  partial-coloring bounds applied to all `2^(m+n)` rank-one constraints still
  have order `sqrt(mn(m+n))`, a leading `N^(3/2)` term.  They gain a power
  only if a theorem first proves that the tight-margin constraint family has
  smaller hereditary discrepancy or active rank; that missing geometry is
  not contained in the energy histograms.
- Generic partial coloring does not yet prove the cube-boundary lemma.  With
  `Theta(N^2)` bridge variables and `2^Theta(N)` constraints, its unweighted
  discrepancy scale is `sqrt(N^2 log(2^N))=Theta(N^(3/2))`.  This consumes a
  leading amount of the slack in (A.10a), rather than a power-saving defect.
  It can gain leverage only from a proved entropy/overlap bound for the
  tight-margin rank-one rows.  Merely counting child energy levels supplies
  no such overlap geometry.
- Fixed-temperature pressure and ordinary zero-diagonal tensor products are
  existing no-go routes.  A diagonal-one tensor `R tensor S` does stay in the
  signing class after subtracting the identity and composes spectral norm
  exactly, but landing it near `M_n` would require the still-unproved matching
  lower constant.  It therefore supplies an explicit upper family, not a
  reduction of the convergence problem.

## 6. Recommendation

Retain Proposition B as the independent track's falsifiable target.  The next
test should use held-out exact child pairs and solve variance-oriented
fractional bridge models, recording `(T,V)` and active rank.  Evidence of
`V/(mn)` decreasing at a power rate would isolate the uniform lemma (A.14);
a scale-stable positive lower bound for `V/(mn)` at the required `T` would
falsify this mechanism.  The independent deep-hole inequality (A.3) and
Bowlin mapping (A.6)--(A.8) should be retained as verified normalization and
obstruction results, but not promoted as a convergence route.

## 7. Bounded primary-literature survey

The mapping used for this survey was fixed before searching: the project is
the covering radius of `D_n=C(K_n)+<1>`, via (A.1), or equivalently the
minimum over signings of a **two-sided** Boolean quadratic maximum.  A theorem
about `rho(C(K_n))` alone is one-sided frustration and does not apply after
adjoining the complement word.

1. Sole and Zaslavsky, [*A Coding Approach to Signed
   Graphs*](https://doi.org/10.1137/S0895480189174374), identify switching
   classes with cosets of the cut code and the maximum one-sided frustration
   with `rho(C(G))`.  Their general bounds concern that one-sided radius.
   For `K_n`, this is not our code: adjoining `1` turns the all-negative word
   from a candidate deep hole into a codeword.  No theorem in the paper gives
   an asymptotic limit or composition law for `rho(D_n)`.
2. Bowlin's primary theorem is imported exactly in Section 2, but it concerns
   `K_(l,r)`.  It determines the standalone rectangular bridge in a periodic
   fixed-shore regime and gives only a leading-order lower bound for balanced
   shores.  It does not impose the nonuniform child-energy margins in (A.10),
   so it removes no parent-composition obligation.
3. Girão, Hunter, Wigderson et al., [*Factorization norms and an inverse
   theorem for MaxCut*](https://doi.org/10.1007/s00208-026-03355-2), Theorem
   1.7, show that an unsigned graph with `m` edges and MaxCut at most
   `m/2+alpha sqrt(m)` contains a clique of size
   `2^(-O(alpha^9)) sqrt(m)`.  Their Lemma 9.4 also bounds ordinary MaxCut
   surplus below by a constant times adjacency trace norm.  Writing a signing
   as a negative-edge graph `G` gives the exact cut formula

   ```math
   H_A(U)=\binom n2-2e(G)-2|U|(n-|U|)+4e_G(U,U^c).       \tag{A.17}
   ```

   Thus ordinary MaxCut is only one term, coupled to the cut cardinality and
   the signing's density.  More decisively, our relevant surplus is
   `Theta(n^(3/2))` while `sqrt(m)=Theta(n)`; in Theorem 1.7 this makes
   `alpha=Theta(sqrt(n))` and its clique conclusion vacuous.  The theorem is
   verified but has the wrong scale and does not remove a present obligation.
4. Dense graphon ground-state/testability theorems normalize energies by
   `n^2`; the present signed minimizers have a vanishing first-order graphon
   signal and a second-order `n^(3/2)` objective.  Higher-order graphon
   fluctuation results concern finitely many motif counts of graphs sampled
   from a fixed graphon, not adversarial signings or the maximum over all
   Boolean spins.  The previously proved planted `n^(3/4)` block also shows
   directly that weak/action profile convergence can miss a constant in our
   normalization.  These theories therefore do not supply the required
   compactness or all-orders realization theorem.

**Survey conclusion.**  No primary source found in this bounded search states
an asymptotic limit for `rho(C(K_n)+<1>)`, a second-order graph-limit theorem
continuous for this Boolean maximum, or a composition inequality with a
summable defect.  The search returns one-sided frustration, rectangular
Gale--Berlekamp bounds, ordinary MaxCut, or first-order graph limits.  None is
promoted.  The cube-boundary penetration lemma remains the only new
falsifiable interface from this track.
