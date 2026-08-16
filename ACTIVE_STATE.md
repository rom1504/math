# Active research state

Evidence cutoff: the first focused AR checkpoint, ledger Section 10.143
(2026-08-16), based on `c4bc2e0`.
This is compact working context.  Use `ledger.md` and Git history only when an
assignment explicitly calls for archive comparison or proof reconstruction.

## Exact problem

For a symmetric zero-diagonal matrix `A=(a_ij)` with off-diagonal entries in
`{+1,-1}`, put

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad
Q(A)=\max_{x\in\{\pm1\}^n}|H_A(x)|,
\qquad
M_n=\min_A Q(A).
```

Determine whether `M_n/n^(3/2)` converges.  Convergence to any constant is a
solution.  Genuine nonconvergence requires two infinite subsequences separated
by a fixed positive normalized gap, or an equivalent strict
`liminf < limsup` proof.

Let `E=binom(n,2)` and let the augmented cut/coboundary code be

```math
\mathcal C_n^+
=\{(\sigma x_ix_j)_{i<j}:\sigma,x_i\in\{\pm1\}\}.
```

Under the usual sign-to-bit identification,

```math
Q(a)=E-2d(a,\mathcal C_n^+),
\qquad
M_n=E-2\rho(\mathcal C_n^+).
```

Thus the problem is equivalently the asymptotic antipodal covering-radius
deficit of the complete-graph cut code.  Be careful: one-sided frustration or
maximum-cut parameters are not this absolute quadratic maximum.

## Rigorous frontier

```math
0.336493364431\ldots
\le \liminf_{n\to\infty}{M_n\over n^{3/2}}
\le \limsup_{n\to\infty}{M_n\over n^{3/2}}
\le {1\over2}.
```

Conference matrices and dense principal restrictions give the all-order
upper limit `1/2`.  The exact values currently recorded for orders `3` through
`14` are

```math
(M_3,\ldots,M_{14})=(3,4,4,5,9,10,12,13,17,18,20,21).
```

These data are useful falsifiers, not asymptotic evidence by themselves.

## Strongest reusable obstructions

1. **Full-state bridge identity.**  For a two-block signing with internal
   matrices `A,B` and bridge `R`, the parent maximum keeps the internal
   energies and `x^T R y` coupled.  Uniform scalar composition statements
   repeatedly hide rather than remove this parent optimization.  Random or
   independently selected bridges retain a leading rectangular-norm floor.

2. **Separately paid channels.**  Scalar finite-fibre atom decompositions and
   symmetric left/right atoms bounded one at a time incur a fixed leading
   multiplier.  Ordinary polarization can be asymptotically sharp.  A valid
   escape must preserve cancellation until after channels are recombined.

3. **Canonical same-map rounding.**  Same-map Gaussian/Krivine correlations
   have nonnegative odd Hermite coefficients.  On the natural conference Gram
   representation this class cannot provide the needed lower-bound constant.
   An asymmetric construction still needs a rigorous same-spin recoupling
   theorem; calling it Grothendieck does not supply one.

4. **Bounded moments and local data.**  Fixed-level SOS, fixed replica depth,
   bounded Eulerian-cycle data, and bounded restriction profiles miss planted
   zero-entropy Boolean resonances of leading size.  At the required scale the
   degree is `Theta(n)` and naive storage recovers the complete signed coset
   histogram.

5. **Moving kernels need hidden size and a root theorem.**  A proved
   operator-valued same-switch inequality preserves joint matrix
   cancellation.  Constant roots nevertheless scalarize exactly.  Correct
   coefficient requires hidden support `exp(Omega(n log n))`; direct matching
   Fourier support and pure matching add/delete representations are
   quantitatively too weak.  A different algebraically closed orbit hierarchy
   is not ruled out.

6. **Entropy-weighted bridges are exact but unclosed.**  The Gibbs variational
   and edge-reveal Rényi identities turn a rare-event cost into relative
   entropy.  Fixed small tilt fails on conference children, while known
   isolated algebraic bridges cost `Theta(n^2)` entropy.  A useful law would
   require only `O(n)` entropy and a `Theta(n)` pressure gain without encoding
   full backward dynamic programming.

7. **High temperature is neighboring, not sufficient.**  Minimized pressure
   has exact finite-temperature formulations and strict-high-temperature
   stability results.  Convergence only below a fixed inverse-temperature
   threshold does not control the ground state.  A route must extend to every
   fixed temperature and then justify the zero-temperature limit.

8. **Arithmetic examples do not prove nonconvergence.**  Bent, Paley,
   symplectic, conference, and Hadamard families supply selected constructions,
   generally at the `1/2` scale.  Nonconvergence needs both a strict upper
   family and an all-signings lower obstruction on separated infinite epochs;
   finite residue effects and route counterexamples do not count.

9. **Action compactness is subsequential.**  Under a common `2 -> 2` bound,
   the same-spin objective is quantitatively continuous in action distance:
   `|Phi(S)-Phi(T)|<=5C sqrt(2d_M)+2d_M`.  Fixed-`C` spectral
   regularization is open, but even it only identifies the scalar and action
   cluster sets.  Convergence still needs lossless exact sign realizers at
   every sufficiently large order, the archived `AR` obligation.

10. **Terminal drift is stronger than it looks.**  Square-field Paley
    conferences are strict one-edge cap minima: every edge flip raises `Q` by
    two, so their cosets have `b=0` and normalized value tending to `1/2`.
    Any uniform terminal drift law with a unique zero therefore forces that
    zero, and the full limit, to equal `1/2`.  The statistic is compressed,
    but its proposed theorem contains the sharp lower-constant burden.

11. **Projective exchangeability is too random.**  A jointly exchangeable
    infinite sign array with tight normalized operator norm must have zero
    Aldous--Hoover mean kernel and hence iid Rademacher edges.  An online
    greedy spin then gives normalized cap at least
    `(2/3)sqrt(2/pi)>1/2` almost surely.  Any viable all-order recovery law
    must therefore be order-dependent and nonprojective.

12. **Low-entropy barycentres cannot be sign-near.**  If `W=E_mu A` and `U`
    is the uniform edge law, then
    `D(mu||U)>=(binom(n,2)-V(W))/2`, where
    `V(W)=sum_(i<j)(1-w_ij^2)`.  Thus a sign-near barycentre
    `V(W)=o(n^2)` costs `Theta(n^2)` relative entropy; an `O(n)`-cost rare
    bridge or microcanonical tilt cannot produce it.

## Campaign status

The retrieval-grounded independent panel is complete:

- 103 primary-source cards across six independently scouted domains;
- six ledger-blind specialist translations and two contrarians;
- full-archive implication audits;
- one foreign-packet experiment per specialist; and
- final independent verification of the theorem artifacts and classifications.

No class-A route survived that panel.  The user has now authorized one focused
campaign on the surviving all-order realization bottleneck, not a restart of
broad route generation.

The weakest noncircular structural target is `AR_min^->`: for each fixed
purification tolerance, choose one bounded-operator liminf action cluster
`T`; on an upward ratio-dense set of orders construct exact hollow signings
whose one-profiles are directed-close to that of `T`, with
`D sqrt(delta)=o(1)`.  Full action convergence, reverse profile inclusion,
higher joint profiles, every-order recovery, and a norm bound uniform in the
purification tolerance are unnecessary.  Principal deletion fills the
remaining `o(n)` upward gaps.  The implication and all error budgets have
been independently verified; see
`artifacts/minimal_all_order_action_recovery.md`.

The strongest new executable reduction is **sign-near weighted recovery**.
For a symmetric hollow `W in [-1,1]^(m x m)`, put

```math
V(W)=\sum_{i<j}(1-w_{ij}^2).
```

It is sufficient, on upward ratio-dense orders, to construct weighted models
of the selected cluster with the correct directed one-profile and
`V(W)=o(m^2)`.  After deleting `o(m)` exceptional rows, biased sign rounding
has normalized operator error `o(1)` and yields exact signs.  For the scalar
version, the sharper direct bound is

```math
Q(A)\le Q(W)+C\bigl(\sqrt{mV(W)}+m\bigr).
```

This rigorously removes the exact-integrality obligation.  It is not yet
known to reduce the universal Boolean-profile obligation: naive blow-ups have
`V=Theta(m^2)`, while low-entropy random laws cannot have a sign-near
barycentre.  See `artifacts/sign_near_weighted_recovery.md`.

## Proposal standard

Every candidate must state:

1. its exact native translation;
2. known theorem(s) plus one boxed new lemma implying convergence or genuine
   nonconvergence;
3. why that lemma contains demonstrably less information than full
   Boolean/coset optimization;
4. an exact finite or structural falsifier;
5. all hypotheses and normalizations of imported theorems; and
6. assumptions that would make the argument circular.

A new name, a finite cap, a solver timeout, an equivalent sufficient
condition, or another class-specific falsifier is not primary progress.

## Held targets and restart conditions

| Target | Final status | Condition for reconsideration |
|---|---|---|
| Fixed-`C` spectral regularization `SR` | B support only | A new lossless every-order realization theorem demonstrably weaker than scalar optimum comparison |
| Sign-near weighted recovery `WAR` | active, independently verified reduction | Construct `W_m` on upward ratio-dense orders with `V(W_m)=o(m^2)` and the selected cluster's scalar upper bound or directed one-profile, without solving the full target-order landscape |
| Directed extremal recovery `AR_min^->` | reduced to `WAR`, otherwise open | A direct exact-sign recovery theorem controlling the one-sided one-profile with no leading residual and no full parent optimization |
| Projective exchangeable recovery | D / rigorously obstructed | Only reconsider a genuinely order-dependent, nonprojective law; projective consistency itself is incompatible with extremality |
| Terminal coset drift | C / falsifier only | A mechanism not forced by Paley traps to prove the sharp `1/2` lower theorem |
| Adversarial pressure | C | A signed deterministic interpolation remainder surviving the outer minimum with `o(n)` defect at every fixed temperature |
| Growing arbitrary-root hierarchy | C/D | A concrete algebraically closed subexponential state proved not to determine the coset histogram |
| Genuine nonconvergence | logically open | Both a selected all-spin low tower and a universal all-signings high theorem with fixed separation |

The full convergence problem remains unsolved and the rigorous interval is
unchanged.  A new autonomous campaign should not begin without one of the
restart inputs above or another comparably strict mathematical reduction.
