# A full-rank bridge compressed by deterministic rearrangement synchronization

Status: main-agent proof draft for independent audit.

Low rank is not the only way a quadratic bridge can have a small response
algebra.  This note gives a natural full-rank dense family in which a common
rearrangement makes all pairwise overlap features functions of scalar
magnetizations.

## 1. Model

For `a=1,...,m`, let `x^a in {-1,1}^n`, put

```math
k_a=|\{i:x_i^a=1\}|,
\qquad s_a=2k_a-n,                                             \tag{SY.1}
```

and let the internal landscape of block `a` be an arbitrary function
`h_a(k_a)`.  Couple block pairs by

```math
(x^a)^TR_{ab}x^b,
\qquad R_{ab}=\alpha_{ab}I+\beta_{ab}J,\quad \alpha_{ab}\ge0. \tag{SY.2}
```

The bridge is dense when `beta_ab!=0`, and it is full rank exactly when
`alpha_ab!=0` and `alpha_ab+n beta_ab!=0`.

### Theorem SY.1 (common nested rearrangement)

For arbitrary real `beta_ab` and nonnegative `alpha_ab`, the exact optimum is

```math
\max_{k_1,...,k_m\in\{0,...,n\}}
\left\{
 \sum_a h_a(k_a)
 +\sum_{a<b}\left[
   \alpha_{ab}(n-2|k_a-k_b|)
   +\beta_{ab}(2k_a-n)(2k_b-n)
 \right]
\right\}.                                                     \tag{SY.3}
```

For every fixed tuple `(k_1,...,k_m)`, all pairwise overlap upper bounds in
(SY.3) are attained simultaneously by one configuration family.

#### Proof

Let `P_a={i:x_i^a=1}`.  For a fixed pair of sizes `k_a,k_b`,

```math
(x^a)^Tx^b=n-2|P_a\triangle P_b|
\le n-2|k_a-k_b|.                                             \tag{SY.4}
```

Because every `alpha_ab` is nonnegative, summing (SY.4) gives an upper bound
on the complete identity-channel contribution.  Relabel coordinates so that

```math
P_a=\{1,...,k_a\}                                             \tag{SY.5}
```

for every block.  These plus sets are nested after ordering by size, and
equality holds in (SY.4) for every pair at once.  The `J`-channel is already
`s_as_b`, independent of the arrangement.  Adding the internal energies and
maximizing over the weights proves (SY.3). `square`

### Corollary SY.2 (a strict full-rank response quotient)

The microscopic configuration space has `2^(mn)` points, but its exact
optimum under the declared symmetric internal landscapes and bridges is
determined by the `(n+1)^m` magnetization tuples in (SY.3).  A maximizing
microscopic family can always be recovered from the tuple by (SY.5).

For two blocks the reusable left response to a symmetric right block is the
`n+1`-entry table

```math
p_h(\ell)=\max_{0\le k\le n}
\{h(k)+\alpha(n-2|k-\ell|)
       +\beta(2k-n)(2\ell-n)\}.                               \tag{SY.6}
```

By the bridge-query isometry, (SY.6) is the canonical exact semantic response
under arbitrary symmetric right-side weights: every table coordinate is an
available query.  The table has polynomial rather than exponential length
even when `R=alpha I+beta J` has rank `n`.  This is not a lower bound against
arbitrary symbolic encodings of the table.

## 2. What synchronized

An unrooted magnetization normally forgets every overlap.  Here the sign
condition `alpha_ab>=0` supplies an **optimizer-compatible common section**:
for every tuple of quotient labels `(k_a)`, the nested representative (SY.5)
simultaneously realizes the best value of every discarded pairwise overlap.
Thus the missing rooted features become deterministic functions

```math
q_{ab}^{\max}(k_a,k_b)=n-2|k_a-k_b|.                           \tag{SY.7}
```

This is stronger than pairwise quotient compatibility: one section realizes
all pairwise maxima globally.  It is also different from low rank: the
identity channel retains `n` algebraic directions, but symmetry plus the
common rearrangement makes only one orbit label per block response-relevant.

The sign hypothesis is structural.  With mixed positive and negative
`alpha_ab`, different pairs may demand nested and anti-nested plus sets, and
the pairwise bounds can be frustrated on a cycle.  Then magnetizations alone
need not determine the optimum by a sum of separately optimized pair
potentials; overlap holonomy reappears.

### Proposition SY.3 (mixed-sign overlap holonomy is extensive)

Let `m=3`, let `n` be even, pin every block to `k_a=n/2`, and take identity
couplings

```math
\alpha_{12}=\alpha_{23}=1,\qquad \alpha_{13}=-1,\qquad\beta=0. \tag{SY.9}
```

Optimizing each pair separately predicts `3n`, but the true joint optimum is
`n`.

#### Proof

At each coordinate the interaction is

```math
x^1x^2+x^2x^3-x^1x^3,
```

which is at most one (its other possible value is `-3`).  Hence the total is
at most `n`.  Taking all three blocks equal to the same balanced spin vector
attains `n`.  Separately, the two positive pairs can each attain dot product
`n`, while the negative pair can attain dot product `-n`, giving the false
sum `3n`. `square`

Thus summing pairwise conditioned optima does not suffice without a common
realizing section.  A three-cycle sign is a composition-created compatibility
bit for that relaxation, and its omission costs the full leading scale in
this benchmark.  A richer joint magnetization-conditioned table could still
represent the unbalanced system exactly.

### Theorem SY.4 (signed balance is the exact gauge criterion)

Let the nonzero identity-channel coefficients define a signed graph on the
blocks.  If every signed cycle has positive sign product, then there are
gauges `epsilon_a in {-1,1}` with

```math
\operatorname{sgn}(\alpha_{ab})=\epsilon_a\epsilon_b.         \tag{SY.10}
```

After replacing `x^a` by `epsilon_a x^a`, Theorem SY.1 applies with every
identity coefficient equal to `|alpha_ab|`.  Thus the exact quotient remains
one transformed magnetization label per block; negative edges do not enlarge
the state when their sign holonomy is trivial.

Conversely, take an isolated unbalanced signed cycle of length `ell`, with
unit-magnitude identity couplings and even `n`.  Pin every block on that
cycle to zero magnetization.  The sum of the individually optimal edge
responses is `ell n`, while the true joint optimum is

```math
(\ell-2)n.                                                     \tag{SY.11}
```

#### Proof

The equivalence between positive sign on every cycle and the vertex gauge
(SY.10) follows by fixing one vertex sign and transporting it along paths;
cycle balance makes the result path-independent.  A spin flip sends
`k_a` to `n-k_a` and changes the incident identity and `J` coefficients by
the same gauge, after which all identity coefficients are nonnegative and
the nested-section proof applies.

On an unbalanced cycle, the product of the realized edge products
`x_i^ax_i^b` at one coordinate is `+1`, whereas the product of the desired
edge signs is `-1`.  At least one edge is therefore unsatisfied.  With unit
weights, the coordinate reward is at most `ell-2`.  Choose a sign assignment
with exactly one unsatisfied edge on half the coordinates and its global
negative on the other half.  Every block is balanced and every coordinate
attains `ell-2`, proving (SY.11). `square`

This gives an intrinsic and checkable deterministic synchronization law:
the discarded overlaps admit one common section attaining every edgewise
optimum iff the signed interaction cocycle is gauge-trivial on cycles.
Failure produces a realizable extensive holonomy witness against the
separable pair-potential algebra, not a lower bound on every possible
response representation.

### Theorem SY.5 (thermodynamic limit for the synchronized full-rank class)

Fix the number of blocks `m` and a signed-balanced coefficient graph.  After
the gauge in Theorem SY.4, write its identity coefficients as
`a_ab=|alpha_ab|`.  Suppose the dense rank-one coefficient is scaled as
`b_ab/n`, and suppose the permutation-invariant internal landscapes satisfy,
uniformly over spins,

```math
{1\over n}h_{a,n}(x)=f_a\left({1\over n}\sum_i x_i\right)+o(1), \tag{SY.12}
```

where every `f_a` is continuous.  Then the normalized optimum converges:

```math
\lim_{n\to\infty}{1\over n}\max_{x^1,...,x^m}
\left\{
 \sum_a h_{a,n}(x^a)
 +\sum_{a<b}\left[
  a_{ab}(x^a)^Tx^b+{b_{ab}\over n}s_as_b
 \right]
\right\}

=\max_{u\in[-1,1]^m}
\left\{
 \sum_a f_a(u_a)
 +\sum_{a<b}\left[
  a_{ab}(1-|u_a-u_b|)+b_{ab}u_au_b
 \right]
\right\}.                                                     \tag{SY.13}
```

For the original signed-balanced system, apply the vertex gauge to the
arguments of `f_a` and to the `b_ab` coefficients.

#### Proof

Theorems SY.1 and SY.4 reduce the finite problem exactly to magnetizations
`u_a=s_a/n` on the parity grid of `[-1,1]`.  The identities

```math
{1\over n}(n-2|k_a-k_b|)=1-|u_a-u_b|,
\qquad {s_as_b\over n^2}=u_au_b                              \tag{SY.14}
```

turn the normalized finite objective into the right side of (SY.13), plus a
uniform `o(1)`, restricted to that grid.  The grids become dense and the
limiting objective is continuous on a compact cube, so their maxima converge
to its maximum. `square`

This is a restricted but genuine thermodynamic-limit theorem generated by
the response/synchronization theory.  It is not a consequence for arbitrary
dense signings: permutation invariance and signed-balance provide exactly the
common optimizer section absent there.

## 3. Complexity and scope

For fixed number of blocks `m`, the quotient description is polynomial in
`n`.  If `m` grows, the number of magnetization tuples is

```math
\exp(m\log(n+1)),                                              \tag{SY.8}
```

so the quotient is sub-landscape but not bounded-state.  Formula (SY.3)
still provides an exact growing factor algebra: adding a block adds one
scalar label and pair potentials given explicitly by (SY.3).  It is not a
fixed-size aggregate state when the number of blocks grows.

The theorem applies to arbitrary permutation-invariant internal functions,
not only Curie--Weiss polynomials.  It does not apply to arbitrary internal
landscapes, arbitrary dense bridges, or mixed-sign identity channels.  Its
theory-level lesson is that algebraic rank can be defeated by a common
optimizer section; without such simultaneous realizability, orbit
compression is only a formal quotient and can miss composition-created
information.
