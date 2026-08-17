# Independent audit: contracting-fibre/cocycle decomposition

**Verdict: PASS WITH REPAIRS.**  The mathematical core of CFC.1, the
stationary-flow dual, CFC.2, and the nonlinear secant lift is correct.  I did
not find a counterexample to any displayed quantitative inequality.  The
remaining repairs concern the acyclic case, the scope of the global
potential statement, and three benchmark interpretations.  In particular,
the rank-one max-plus benchmark is a consistency check only after its
compatibility graph has been made visible; it is not automatically a direct
`rho=0` instance of the common-invariant-law hypothesis.

## Checks performed

I read the proof independently and compared it with canonical Theorems
16.18, 17.1e, 17.1h, 17.1l, 17.7, the approximate residual-shell draft, and
the expander phase-refresh theorem.

The supplied verifier passes:

```text
contracting-fibre/cocycle checks passed: 7113
```

I also ran two checks not present in that verifier:

1. 300 random chains of genuinely rectangular kernels with changing fibre
   dimensions and transported laws `pi_(s+1)=pi_s P_s`; every instance
   satisfied (CFC.6) with the independently computed largest centred
   singular-value factor.
2. 80 random complete directed weighted graphs of orders two through five;
   a stationary-circulation linear program agreed to numerical tolerance
   with the maximum absolute simple-cycle mean in (CFC.14).

These additional checks matter because the bundled verifier uses one square
kernel and does not exercise the changing-fibre claim or the flow dual.

## Audit of CFC.6--CFC.12a

The orientation is correct.  With

```math
P_e:L^2(\pi_{q'})\longrightarrow L^2(\pi_q),
\qquad \pi_qP_e=\pi_{q'},
```

`P_e` maps constants to constants and target-centred functions to
source-centred functions.  Hence the expansion (CFC.12) is exact.  A reward
at position `s` passes through `s-1` kernels, giving `rho^(s-1)B`, while the
terminal residual passes through all `t` kernels, giving `rho^tR`.  Thus the
geometric series and every coefficient in (CFC.6) are correct.

The added identity (CFC.12a) is also correct and useful: the scalar and
centred channels are orthogonal in `L^2(pi_(q_0))`.  It closes the only
possible cancellation loophole in the cycle lower bound.

The sharpness examples are valid.  A single two-state loop with the
mean-zero eigenvector `(1,-1)` aligns the terminal and fresh residual terms;
an arbitrary constant edge reward supplies the scalar channel.  The revised
wording no longer asserts an invalid triangle-equality statement for a
Cartesian product.

## Audit of the cycle asymptotics and terminal means

For a graph containing a directed cycle, (CFC.8) is correct.  Loop erasure
leaves at most `|Q|-1` edges, and every erased closed walk decomposes into
simple cycles.  The upper bound is therefore

```math
t\chi_G(m)+(|Q|-1)M+U+R+{B\over1-\rho}.
```

Repeating a maximizing cycle gives the reverse limsup; (CFC.12a) prevents
the centred term from cancelling it.  The additions of `U` to (CFC.9) and
(CFC.11) repair the terminal-mean omission.

Two scope clarifications are still needed.

1. **Acyclic graph.** If `G` has no cycle, there are no paths of length `t`
   for all sufficiently large `t`, so the left side of (CFC.8) is not
   literally defined under the usual `sup emptyset=-infinity` convention.
   State that the asymptotic assertion is vacuous in this case, or explicitly
   define the displayed nonnegative supremum as zero when the path set is
   empty.
2. **SCC potential versus global potential.** Vanishing directed-cycle sums
   imply a vertex potential on each strongly connected component, but not
   necessarily on the whole directed graph.  A directed acyclic diamond can
   have two coterminal paths with unequal sums and no cycle at all.  This
   causes only a bounded transient error, so the `chi=0` boundedness
   criterion remains correct.  However, (CFC.11) with just `osc(psi)` is a
   global bound only when a global potential exists (in particular, on a
   strongly connected graph).  Otherwise add a finite transient constant or
   explicitly scope (CFC.10)--(CFC.11) to one SCC.

The phrase “after calibrating terminal means” should be interpreted as
retaining the scalar terminal value `bar u_q`.  If it is not retained, the
uniform shell includes `U`.  Section 5.1 currently quotes (CFC.20) without
`U`; it should say “after terminal means are calibrated” or display
`U+R+B/(1-rho)`.

## Audit of CFC.13--CFC.14

The flow dual is exact.  A normalized nonnegative circulation is a convex
combination of normalized directed-cycle flows, and maximizing the absolute
value of a linear functional reduces to one of the two signed extrema.  The
new empty-polytope convention repairs the acyclic case.

On a non-strongly-connected graph, however, the extensive datum is more
precisely the functional induced on stationary flows, or the recurrent-SCC
cohomology class.  Calling it the global class
`[m] in C^1(G)/B^1(G)` overcounts acyclic coterminal-path defects, which can
fail to be global coboundaries but cannot be pumped.  Section 5.1 should use
“recurrent cohomology/stationary-flow functional,” unless strong
connectivity is assumed.

## Audit of nonlinear Corollary CFC.1a

The trajectory orientation is correct if

```math
x_s=F_(e_s)(x_(s-1)),
\qquad y_s=\widehat F_(e_s)(y_(s-1)).
```

Adding and subtracting `F_(e_s)(y_(s-1))` gives exactly (CFC.14f).  Since
every secant preserves `pi`, its action on the centred part contracts by
`rho`, while (CFC.14d) increments the mean by the visible scalar `m_e`.
Induction proves (CFC.14e) with the stated constants.  Along a repeated
visible cycle, its scalar mean grows at exactly the cycle rate and the
centred term stays bounded, so the cycle asymptotics really do lift to this
nonlinear class.

Two wording repairs would make the result exact.

1. “After subtracting its initial mean” is too short: (CFC.14e) subtracts
   the propagated scalar component
   `pi(x_0-y_0)+sum_s m_(e_s)`.
2. “Optimizer switches and ties cause no extra information growth” is true
   only for paired trajectories following the same declared visible path
   and under the uniform common-law hypotheses (CFC.14b)--(CFC.14d).  Keep
   those qualifiers in the conclusion.

The max-plus scope paragraph is mathematically accurate: all-finite
max-plus maps have row-stochastic secants, but a common invariant `pi` and a
strict centred `L^2(pi)` factor are not automatic.  In fact this distinction
is essential.  An exact max-plus rank-one reset has identical-row secants,
but their row law can depend on the active selector; such a secant need not
preserve any one fixed `pi`.  The lost projective coordinate reappears as a
hidden-state-dependent scalar gain.

Consequently Section 5.2 should not call rank-one max-plus reset directly
“the `rho=0` endpoint” without qualification.  It becomes a singleton-fibre
`rho=0` instance **after** the previous/right-profile type and its directed
compatibility are promoted to the visible control graph.  At that point
CFC.7--CFC.14 recover the scalar rate, but the construction of that visible
graph is imported from ARS.2/Theorem 17.7 rather than discovered by CFC.1a.

## Audit of CFC.17

The variance tax is correct, including its constants.

* `h=f-Pf` has mean zero, upper endpoint `epsilon`, and lower endpoint
  `-osc(f)`.
* The elementary bounded-variance inequality gives
  `||h||_2^2 <= osc(f) epsilon`.
* On centred functions,
  `(1-rho)||f-Pi f||_2 <= ||f-Pf||_2`.
* Orthogonal centring and `||f-g||_infinity<=omega` give
  `||f-Pi f||_2 >= sigma_g-omega`.
* Finally `osc(f)<=osc(g)+2omega`.

These yield exactly (CFC.17).  The zero-denominator convention is harmless:
then `g` is constant and `omega=0`, so the numerator also vanishes.

The Walsh/semantic-phase paragraph needs one scope guard.  CFC.2 assumes a
single same-scale inequality `f<=Pf+epsilon`.  A general scale transfer has
`f_r<=P_rf_(r+1)+epsilon_r`; applying CFC.2 then needs either a homogeneous
response, or a separate uniform comparison of `f_r` and `f_(r+1)` (which
adds that comparison error to the toll).  The time-inhomogeneous ER.2
theorem remains the correct unconditional phase-refresh statement.  Thus
CFC.23 is a valid homogeneous/mass-sensitive complement, not by itself a
replacement for ER.2.

## Discounted benchmark correction

The discounted paragraph currently merges the scalar and centred
denominators.  For

```math
D_s=a_s+\lambda P_sD_(s+1),
```

the scalar channel has factor `lambda`, while the centred channel has factor
at most `lambda rho`.  The respective fresh-error shells are

```math
{M\over1-\lambda}
\qquad\hbox{and}\qquad
{B\over1-\lambda\rho}.
```

There is no undiscounted cycle drift when `lambda<1`.  Replacing the current
sentence with this two-channel formula would preserve the point and avoid
saying that one simply replaces every occurrence of `rho` by a product.

## Novelty and overlap assessment

The novelty claim is defensible only at the level stated below.

* Geometric damping of centred Markov rewards is classical and closely
  parallels the fresh-residual part of Theorem 16.18.
* CFC.7--CFC.14 are the scalar specialization of the cycle/coboundary and
  cycle-LP mechanisms already canonicalized in Theorems 17.1e, 17.1h, and
  17.1l.
* The rank-one max-plus example is already Theorem 17.7/ARS.2.
* CFC.17 is a new repository-level quantitative bridge between one-sided
  toll, variance, and spectral gap, but its proof uses elementary standard
  ingredients and should not be advertised as literature-level novelty
  without a literature check.

What is genuinely new and useful **inside this program** is the joint
decomposition and its nonlinear common-secant corollary: under one explicit,
checkable invariant-law hypothesis, all hidden centred reward errors are
uniformly bounded and the complete extensive obstruction is a finite visible
cycle functional.  This unifies two previously separate resources and gives
a new benchmark theorem.  It does not establish an information-theoretic
minimal state, and “strictly less information” should be qualified as an
approximation at error `B/(1-rho)+R`, not an exact simulator claim.

## Verifier gaps worth repairing

The verifier is useful but currently checks a strict subset of the draft.
Before canonicalization, add or document checks for:

1. rectangular kernels and transported changing fibre laws;
2. the orthogonality identity (CFC.12a);
3. the stationary-flow LP dual (CFC.14);
4. potential/cycle equivalence on a strongly connected graph and the
   acyclic-diamond guardrail;
5. a direct nonlinear secant example for CFC.1a;
6. the exact discounted two-denominator formula if that benchmark remains.

Subject to the scope repairs above, I recommend canonicalization.
