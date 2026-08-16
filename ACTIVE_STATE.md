# Active research state

Evidence cutoff: focused AR checkpoint-two consolidation, ledger Section
10.144 (2026-08-16), building on checkpoint commit `77ef709`.
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

1. **Local summaries miss Boolean extremes.**  Fixed-level SOS, bounded
   moments/cycles, bounded restriction profiles, and separately paid scalar
   channels miss planted zero-entropy resonances of leading size.  Retaining
   enough generic data restores the full signed-coset response.  Ordinary
   bridge, fibre, Krivine, and fixed-gadget variants are archived, not live.

2. **Action compactness is subsequential.**  Under a common `2 -> 2` bound,
   the same-spin objective is quantitatively continuous in action distance:
   `|Phi(S)-Phi(T)|<=5C sqrt(2d_M)+2d_M`.  Fixed-`C` spectral
   regularization is open, but even it only identifies the scalar and action
   cluster sets.  Convergence still needs lossless exact sign realizers at
   every sufficiently large order, the archived `AR` obligation.

3. **Sampling and projective exchangeability are too random.**  A jointly
    exchangeable
    infinite sign array with tight normalized operator norm must have zero
    Aldous--Hoover mean kernel and hence iid Rademacher edges.  An online
    greedy spin then gives normalized cap at least
    `(2/3)sqrt(2/pi)>1/2`.  More generally, bounded-op signings have iid
    fixed-size induced limits, with a fixed `0.063846...` gap in `Phi` above
    the relevant extremal scale.  Viable recovery must be nonprojective and
    cannot use uniform mesoscopic induced sampling.

4. **Sign-near weighted recovery is equivalent and information-heavy.**  If
    `W=E_mu A`, `U` is the uniform edge law, and `N=binom(n,2)`, then
    `D(mu||U)>=N[log2-h(V(W)/(2N))]`.  Thus `V(W)=o(n^2)` forces
    `D(mu||U)>=N log2-o(N)` and `H(mu)=o(N)`.  Moreover all but `o(n^2)`
    edge signs are exposed.  Although its terminal rounding theorem is proved,
    exact recovery supplies the converse witness `W=A`; weighted existence is
    not a strict reduction.

5. **Design repair is only a last-mile module.**  An `O(n)`-edge leave is
   negligible for both the normalized cap and action profile.  Existing dense
   design theorems control fixed local statistics, not the universal outer
   profile over exponentially many colorings.

6. **Pressure and shell entropy add a no-gap obligation.**  Fixed-temperature
   pressure can oscillate through state multiplicities even when a zero state
   exists at every order.  No current interpolation theorem reaches all
   temperatures and survives the outer signing minimum.

7. **Tested near-minimizer rigidity does not narrow recovery.**  Vanishing
    fourth-moment/cycle defect misses a Boolean kernel spike even after one
    conference vertex is deleted.  A universal-vertex extension preserves
    normalized near-optimality and bounded operator scale while destroying
    square-field uniform integrability.  Exact order-six examples with the
    same positive labelled row sums have caps `11` and `7`.  Stronger forms
    either fail these tests or restore the complete outer Boolean profile.

8. **Arithmetic examples do not prove nonconvergence.**  Bent, Paley,
   conference, and Hadamard families give selected constructions, usually at
   scale `1/2`.  Nonconvergence requires a universal high obstruction and a
   low construction on separated infinite epochs; finite residue effects and
   failed AR implementations do not count.

## Campaign status

The retrieval-grounded independent panel is complete:

- 103 primary-source cards across six independently scouted domains;
- six ledger-blind specialist translations and two contrarians;
- full-archive implication audits;
- one foreign-packet experiment per specialist; and
- final independent verification of the theorem artifacts and classifications.

No class-A convergence route survived that panel.  The user then authorized
a focused two-checkpoint campaign on all-order action realization.  That
campaign is now consolidated and paused.

The weakest retained **optimizer-free** structural target is
extremal-envelope recovery with uniform integrability (`EER_UI`).  At a fixed
purification tolerance, let `K_eta` be the compact cluster set of one bounded
near-liminf sequence and put

```math
\mathcal E_\eta
=\overline{\bigcup_{T\in\mathcal K_\eta}\mathcal S_1(T)}.
```

On sufficiently covering order sets, construct exact signings whose every
one-profile law is directed-close to `E_eta` and whose energy products `xy`
are uniformly integrable.  Each source profile may match a different cluster
phase; no common target operator is needed.  If `gamma_eta` is the upward
multiplicative covering ratio, the exact tolerance condition is

```math
\gamma_\eta^{3/2}(\alpha+\eta)\longrightarrow\alpha,
\qquad
\alpha=\liminf_n {2M_n\over n^{3/2}}.
```

Weak convergence plus uniform integrability transfers the energy integrals,
and principal deletion covers the omitted orders.  A quantitative sufficient
form uses normalized operator bounds and `D sqrt(delta)->0`.  Matching only
one extremizing profile is formally weaker, but selecting it performs the
target-order Boolean maximum; scalar recovery alone is equivalent to
convergence.  See `artifacts/extremal_envelope_recovery.md`.

Sign-near weighted rounding is complete but no longer a live architecture.
If `V(W)=o(n^2)`, biased rounding gives

```math
\|T_A-T_W\|_{L^\infty\to L^1}=o(1),
\qquad d_M(T_A,T_W)=o(1),
\qquad |\Phi(T_A)-\Phi(T_W)|=o(1).
```

Exact recovery gives the converse weighted witness `W=A,V=0`; weighted
recovery is therefore an equivalent recovery obligation, not a strict
reduction.  It remains a reusable terminal rounding theorem.

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
| Extremal-envelope recovery `EER_UI` | B, weakest optimizer-free profile target; open | A nonlocal constructor or inverse finite-model theorem that controls every source profile without an almost complete sign skeleton |
| Sign-near weighted recovery `WAR` | C as architecture; A rounding module | Reconsider only with an independently proved constructor; existential weighted recovery is equivalent to exact recovery |
| Directed extremal recovery `AR_min^->` | B/C, stronger than EER | Subsumed unless phase coherence itself becomes useful in a new theorem |
| Projective exchangeable recovery | D / rigorously obstructed | Only reconsider a genuinely order-dependent, nonprojective law; projective consistency itself is incompatible with extremality |
| Uniform mesoscopic induced sampling | D / rigorously obstructed | Only a quantitatively structured diagonal before iid mixing or an optimized subset could evade the theorem |
| Terminal coset drift | C / falsifier only | A mechanism not forced by Paley traps to prove the sharp `1/2` lower theorem |
| Adversarial pressure | C | A signed deterministic interpolation remainder surviving the outer minimum with `o(n)` defect at every fixed temperature |
| Growing arbitrary-root hierarchy | C/D | A concrete algebraically closed subexponential state proved not to determine the coset histogram |
| Genuine nonconvergence | logically open | Both a selected all-spin low tower and a universal all-signings high theorem with fixed separation |

The full convergence problem remains unsolved and the rigorous interval is
unchanged.  The focused AR campaign is paused after two checkpoints.  Resume
only after an external inverse-recovery theorem, a concrete nonlocal
constructor, a continuous order-class obstruction, or another comparably
strict input.
