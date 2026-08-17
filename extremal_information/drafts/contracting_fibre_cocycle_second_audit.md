# Second independent audit: contracting fibres and cycle-rank response entropy

**Verdict: PASS WITH MINOR CLARIFICATIONS.**  I find CFC.1b mathematically
correct under its stated strong-connectivity, fixed-carrier, and exact-gauge
conventions.  Its path-error accounting and cycle-rank exponent are correct.
The earlier repairs to CFC.1, CFC.1a, CFC.2, the SCC language, and the
discounted benchmark are also correct.  Before canonicalization I recommend
three small clarifications:

1. say that the centred-reward net may be chosen with centred codewords;
2. distinguish the sharp covering lower bound for codewords in the scalar
   cochain model from the same lower bound for arbitrary response encoders,
   where the packing argument changes `epsilon` to `2epsilon`;
3. keep "composable codebook" explicitly scoped to a response code on an
   already certified fixed contracting carrier, modulo exact endpoint gauge.

None changes a displayed upper bound or a claimed exponent.

## 1. Quotient norm and its dimension

For a strongly connected finite directed graph, every edge belongs to a
directed closed walk.  If a cochain `d` has zero sum on every directed simple
cycle, it has zero sum on every directed closed walk.  Fixing a root and a
directed return path then shows that the `d`-sum from the root to any vertex is
path independent.  Consequently

```math
d_e=\psi(t(e))-\psi(s(e)).
```

Thus the kernel of `chi_G` is exactly the image of the directed incidence
map.  Homogeneity and the triangle inequality are immediate from the maximum
over cycle averages, so (CFC.14h) is a normed quotient, not merely a
seminormed one.  Strong connectivity makes the underlying incidence rank
`|Q|-1`, including in the presence of loops or parallel edges, and hence

```math
\dim H_G=|E|-|Q|+1=r_G.
```

This verifies both assertions used by the volume argument.  The same statement
would be false as written on a general non-strongly-connected graph: acyclic
coterminal-path defects lie in the kernel of the recurrent cycle functional
without being global gradients.  CFC.1b correctly avoids that issue by
assuming strong connectivity; the rest of the draft correctly uses recurrent
SCC cohomology in the general case.

## 2. Upper covering bound and exact path-error accounting

Let two response specifications differ by

```math
\max_q\|\Delta u_q\|_2\le\epsilon,
\qquad
\max_e\|\Delta b_e\|_2\le(1-\rho)\epsilon,
\qquad
\chi_G(\Delta m)\le\epsilon.
```

The three contributions are as follows.

* Every Markov prefix is an `L^2` contraction, so the terminal contribution is
  at most `epsilon`.
* CFC.6 sums the centred edge errors to at most
  `(1-rho)epsilon sum_(j=0)^(t-1)rho^j <= epsilon`.
* Choose a bounded linear section `s:H_G -> R^E`.  Norm equivalence in the
  fixed finite-dimensional space gives
  `||s([d])||_infinity <= K_G chi_G(d)`.  Remove directed closed subwalks from
  a length-`t` path.  Each removed closed walk decomposes into simple directed
  cycles and contributes at most its length times `epsilon`; the remaining
  simple path has at most `|Q|-1` edges and contributes at most
  `(|Q|-1)K_G epsilon`.

The scalar error is therefore at most

```math
\epsilon t+(|Q|-1)K_G\epsilon,
```

and the total error is exactly bounded as claimed with
`C_G=2+(|Q|-1)K_G`.  No path-length-dependent centred term is missing.

The standard maximal-separated-set/volume argument in an arbitrary
`r_G`-dimensional normed space gives

```math
N(B_L,\epsilon B)\le(1+2L/\epsilon)^{r_G}.
```

Multiplication by the two dictionary covering numbers proves (CFC.14i).
If `N_B` is defined as an external cover, a centre need not initially be
centred.  It may be orthogonally projected to its edgewise mean-zero part:
because every member of `mathcal B` is centred, this does not increase any
`L^2` distance.  Adding this sentence would close a minor formal gap in the
current proof.

## 3. What "modulo potential" pays for

For two scalar cochains, write

```math
\Delta m=s([\Delta m])+\nabla\psi.
```

The gradient contributes only
`psi(q_t)-psi(q_0)` to a path.  CFC.1b deliberately treats this as exact
terminal/start calibration rather than as a codeword coordinate.  Under that
convention the stated code size is correct.  If absolute endpoint responses
are part of the declared query, the potential is not free: a cover of its
bounded `( |Q|-1 )`-dimensional class must be included.  The draft now says
this explicitly, so there is no hidden finite-state claim.

This point also limits the phrase "composable codebook."  The code is reusable
for every concatenated visible path on the *fixed, already certified* fibre
carrier.  It is not an algorithm that discovers that carrier, a finite
bisimulation of arbitrary hidden dynamics, or a claim that arbitrary endpoint
potentials have zero description length.

## 4. Lower bound and sharpness

For scalar cochains `m,m'`, repeated traversal of a cycle maximizing
`chi_G(m-m')` gives

```math
\limsup_{t\to\infty}{1\over t}
 \sup_{|p|=t}|R_m(p)-R_{m'}(p)|
=\chi_G(m-m').
```

The opposite inequality follows by loop erasure.  Hence the asymptotic
response-rate metric on scalar specifications modulo potentials is exactly
the norm of `H_G`.  Covering its radius-`L` ball by scalar-cochain codewords
therefore needs at least `(L/epsilon)^{r_G}` centres (up to the usual harmless
endpoint/range conventions).  If code centres are allowed to be completely
arbitrary response predictors rather than members of this normed model, one
should instead take a `2epsilon`-packing: two source responses assigned to the
same predictor can be at most `2epsilon` apart in rate.  This still forces

```math
\Omega_G((L/\epsilon)^{r_G})
```

centres, so the asserted exponent and scale are unchanged.  I recommend
making this factor-two distinction in the prose rather than saying every
arbitrary code literally gives an `epsilon`-cover of `H_G`.

The other two scales are distribution-free sharp as claimed.  A two-state
mean-zero eigenmode realizes terminal decay `rho^t R`; repeated aligned fresh
reward realizes `(1-rho^t)B/(1-rho)`.  Thus neither
`N_U(epsilon)` nor `N_B((1-rho)epsilon)` can be uniformly coarsened without an
additional structural assumption.

## 5. Regression audit of CFC.1, CFC.1a, and CFC.2

I reran the supplied verifier in the repository environment:

```text
contracting-fibre/cocycle checks passed: 8099
```

The checks now include rectangular transported laws, orthogonality of scalar
and centred channels, stationary-flow LP duality, the acyclic-diamond
guardrail, nonlinear stochastic secants, and the two discounted denominators.

The analytical constants also check independently.

* In CFC.6, the reward at position `s` crosses exactly `s-1` kernels, while
  the terminal residual crosses `t`, yielding
  `rho^t R+B(1-rho^t)/(1-rho)`.
* Orthogonality under the transported initial law gives (CFC.12a), so a
  centred response cannot cancel scalar cycle drift.
* Loop erasure gives the transient `(|Q|-1)M`; terminal mean and centred
  shells add `U+R+B/(1-rho)` in (CFC.9).
* The stationary-flow formula is exact by finite circulation decomposition.
* CFC.1a is correct only for paired trajectories following the same declared
  visible path and for secants sharing the transported law and strict centred
  contraction.  The revised statement contains all three qualifiers.
* For CFC.2, `h=f-Pf` is mean zero and lies in
  `[-osc(f),epsilon]`, so `E h^2<=osc(f)epsilon`; combining this with
  `||(I-P)(f-Pi f)||_2 >= (1-rho)||f-Pi f||_2` and the recovery radius gives
  exactly (CFC.17).
* The homogeneous/same-scale qualification on the Walsh application is now
  explicit.  The discounted recursion correctly separates
  `M/(1-lambda)` from `B/(1-lambda rho)`.

I found no regression in these repaired statements.

## 6. Novelty and director judgment

CFC.1b is not merely an enumeration of a finite dynamic-programming state.
The hidden fibres may be arbitrarily large, and the theorem gives a
depth-independent reward-resolution scale `(1-rho)epsilon` plus a persistent
cycle-space cost independent of the hidden path tree.  It therefore makes a
real quantitative prediction: static response entropy is sampled at the
forgetting scale, whereas invariant reward memory costs `r_G` persistent
coordinates.

The ingredients are nevertheless classical: Markov contraction, graph
cohomology/circulation duality, and volumetric covering.  The new contribution
inside this project is their exact allocation in one reusable response law,
not a new general minimization theorem for weighted automata.  In particular,
CFC.1b assumes rather than discovers the common-law contracting carrier, and
its dictionary covering numbers may still be as large as the full response
landscape.  It should be canonicalized as a rigorous unifying theorem and
benchmark, not advertised as a necessary-and-sufficient theory of dynamic
compression.

Subject to the three prose clarifications at the start, I recommend
canonicalization.
