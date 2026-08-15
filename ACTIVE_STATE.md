# Active research state

Evidence cutoff: `e1e298db5afcbf939b994876cf08951e479c0cfe`, 2026-08-15.
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

## Active campaign architecture

The selected *research architecture*, not a selected mathematical route, is a
retrieval-grounded independent panel:

- six literature toolkits: discrepancy; spin glasses; coding/association
  schemes; Banach/tensor theory; extremal/probabilistic combinatorics; dense
  graph limits;
- six separate ledger-blind specialists who translate the problem only after
  receiving their domain toolkit;
- two ledger-blind contrarians, assuming respectively nonconvergence and a
  short proof;
- archive verifiers who compare frozen proposals with the complete project;
- one revision and cross-domain/foreign-packet critique for surviving `A/B`
  proposals; and
- director selection of at most two candidates for independent proof and
  disproof attempts.

Do not select moving representations, entropy bridges, growing states, or
pressure continuation merely because they are named here.  They are test
cases for the panel.  Mathematical route selection remains open until archive
classification.

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

## Current live targets and stopping conditions

| Stage | Live target | Stop or advance criterion |
|---|---|---|
| Retrieval | Six toolkits of 8--20 primary sources with precise theorems and hypotheses | Reject abstract-only analogy or unverifiable citation |
| Independent research | At most three native architectures per domain, each with an exact lemma and falsifier | Freeze proposal before archive access |
| Archive audit | Actual implication-level `A/B/C/D` classification | Do not reject by shared vocabulary alone |
| Cross-domain review | Strongest failure reason and strongest foreign rescue theorem | Preserve disagreements; no consensus by imitation |
| Selection | Ranked standardized cards | Execute at most two `A/B` candidates meeting the strict selection test |
| Execution | Minimal theorem, independent proof/disproof, decisive computation only | Stop false routes immediately; after two no-progress checkpoints return to director review |

The full convergence problem is unsolved.  No mathematical architecture is
currently endorsed; the panel itself is the active, user-authorized campaign.
