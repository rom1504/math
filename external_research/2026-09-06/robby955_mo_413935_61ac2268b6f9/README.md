# MathOverflow 413935: Min-max quadratic forms of signs

Research notes, exact computations, and reproducible verification for Paata
Ivanisvili's [MathOverflow question](https://mathoverflow.net/questions/413935/min-max-of-a-quadratic-form-of-plus-minus-ones).

For a symmetric zero-diagonal sign matrix `A`, define

```math
Q_A(x)=\sum_{1\le i\lt j\le n}a_{ij}x_ix_j,
\qquad
M(A)=\max_{x\in\{-1,1\}^n}|Q_A(x)|,
\qquad
F(n)=\min_A M(A).
```

The question is whether `F(n) / n^(3/2)` converges. The limit remains open.

## Film

[![Watch the MO 413935 film](figures/mo-413935-research-film-poster.png)](https://cdn.jsdelivr.net/gh/Robby955/mo-413935-research@6ef37905d92908b7c091e947ea8ce00d933cd60e/media/mo-413935-research-film.mp4)

A 1 minute 50 second visual introduction to the signed-graph game, the exact
frontier through `F(15) = 27`, and the still-undecided case
`F(16) ∈ {28, 30}`.

[Watch the full film (1080p)](https://cdn.jsdelivr.net/gh/Robby955/mo-413935-research@6ef37905d92908b7c091e947ea8ce00d933cd60e/media/mo-413935-research-film.mp4) ·
[Watch the 7-second visual loop](https://cdn.jsdelivr.net/gh/Robby955/mo-413935-research@f88019e7f8561f51d14e6b9450b0bb08c1f9cf92/media/mo-413935-hero-loop.mp4)

## The game

The definition is exactly a two-move, zero-sum game, and this reading is a
restatement, not a new formulation. A designer moves first and labels every
edge of the complete graph on `n` vertices with `+1` (the endpoints want to
agree) or `-1` (they want to disagree). A player moves second, sees the
labels, and assigns `+1` or `-1` to every vertex; each edge scores the
product of its label with its endpoints' signs. The player tries to push
the total score as far from zero as possible, in either direction; the
designer tries to make every assignment score near zero. `M(A)` is the
player's best response to the wiring `A`, and `F(n)` is the value of the
game under best play on both sides.

The two universal bounds have a game meaning. The lower bound
`F(n) >= n sqrt(n-1) / pi` is proved by an averaging argument: a player who
commits to a randomized assignment before seeing the wiring already forces
this much imbalance against every design. The upper bound near
`n^(3/2) / 2` comes from explicit conference-matrix designs against a fully
informed player. The open convergence question asks, in these terms, how
much seeing the board is worth in the limit.

## Focused outputs

The project is now split into two narrower manuscripts:

- [Finite structure in a min-max quadratic sign problem](paper/mo-413935-finite-results.pdf)
  covers the exact values through order 15, computational certificate
  boundaries, optimizer non-heredity, the weighted Bellman identity, and the
  pinned framed-Hadamard theorem proving `F(16) <= 30`.
- [Relative-gauge composition for quadratic sign discrepancy](paper/mo-413935-composition-framework.pdf)
  gives the exact labeled block-composition law, its finite-temperature
  alignment form, the audited scalar obstructions, and the precise sufficient
  cross-order theorem that remains open.

The [complete research note](paper/mo-413935-second-attempt.pdf) and its
[source](paper/second_attempt.tex) remain the broad archive. The tag
`research-frontier-2026-08-05` freezes the exploration frontier before this
publication split. Failed approaches have not been deleted.

The active research scope is in [ACTIVE_RESEARCH.md](ACTIVE_RESEARCH.md). A
self-contained prompt for a public frontier model is in
[COMPOSITION_FRONTIER_PROMPT.md](COMPOSITION_FRONTIER_PROMPT.md). The earlier
[broad prompt](FRONTIER_MODEL_PROMPT.md) is retained for provenance.

## First cases

For `n = 1` there are no edges at all, so every score is zero and
`F(1) = 0`.

For `n = 2` there is a single edge. Whatever its label, its score is `+1`
or `-1`, so every assignment has absolute score exactly 1 and `F(2) = 1`.

For `n = 3` the player defeats the designer completely: whatever the three
labels, the player can make all three edges score the same sign. The reason
is the identity

```math
(x_1x_2)(x_1x_3)(x_2x_3)=1,
```

which says the three products induced by any vertex assignment always have
product `+1`. Given labels `(a_12, a_13, a_23)`, either the labels have
product `+1`, and the player matches all three edges for a score of `+3`,
or their negation does, and the player anti-matches all three for `-3`.
Either way `F(3) = 3`, the maximum possible.

![The two triangle classes](figures/first_cases.svg)

At `n = 4` this stops: the best design holds the worst case to `F(4) = 4`,
below the maximum possible score of 6. The designer gains ground for the
first time, and the rest of the table measures exactly how much ground the
designer can keep gaining as `n` grows.

## Current rigorous frontier

The audited universal bounds are

```math
\frac1\pi
\le \liminf_{n\to\infty}\frac{F(n)}{n^{3/2}}
\le \limsup_{n\to\infty}\frac{F(n)}{n^{3/2}}
\le \frac12.
```

The exact finite sequence is

| `n` | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `F(n)` | 0 | 1 | 3 | 4 | 4 | 5 | 9 | 10 | 12 | 13 | 17 | 18 | 20 | 21 | 27 |
| `F(n)/n^1.5` | 0.000 | 0.354 | 0.577 | 0.500 | 0.358 | 0.340 | 0.486 | 0.442 | 0.444 | 0.411 | 0.466 | 0.433 | 0.427 | 0.401 | 0.465 |

The normalized row is rounded to three decimals; whether it converges
is the open question. The lower certificates at orders 11, 13, and 15 are computer-assisted. Their
completeness trusts nauty's isomorph-free graph generation and canonical
deduplication; all stream counts, digests, producer exits, surviving records,
and explicit witnesses are checked separately. The order-15 value comes from a
threshold tower: no order-15 signing has maximum at most 25, and extending the
order-14 Paley conference matrix by its best sign row gives 27 exactly. The
final tower level is replayed in full from a committed, digest-pinned
catalogue of all 1,313,164 order-14 signings with maximum at most 25. The
order-16 result is only

```math
F(16)\le 30.
```

One rung higher, the order-18 Paley conference matrix has maximum 33, so
puncturing and even order-17 parity give `F(17) <= 32` unconditionally.

Monotonicity from `F(15) = 27` and even order-16 parity leave only
`F(16) = 28` or `F(16) = 30`. Neither value is claimed here. An exact
reduction narrows the first case: any order-16 signing with maximum at most
28 has all 120 of its order-14 pair-deletions at maximum exactly 27. This
follows from an exhaustive completion search over the committed order-14
catalogue: exactly 163 of its 1,313,164 signings admit extension rows of
order-15 maximum 27 (316 rows in all), and no completion of any such pair
reaches maximum 28; the completion minimum is 30, independently matching
the framed-Hadamard ceiling. The 316 extensions collapse to exactly 30
rooted canonical classes of order-15 minimizers, the complete list of
minimizers having a deletion of maximum at most 25
(`verification/order15_minimizers_rooted.g6`). Deciding `F(16)` is thereby
reduced to order-15 minimizers all of whose deletions have maximum 27; none
is currently known. Replay with
`python3 verification/research_order16_reduction.py --jobs 4`.

Other banked results include:

- `F(n) <= F(n+1) <= F(n)+n` and the parity-refined puncturing bound;
- the exact augmented cut-code identity
  `F(n) = choose(n,2) - 2 rho(D_n)`;
- the Gaussian lower bound `F(n) >= n sqrt(n-1) / pi`;
- an exact fixed-density cut-discrepancy equivalence up to `O(n)`;
- an exact energy-weighted covering-radius formula for one-vertex extension;
- finite optimizer non-heredity, including the complete `2+8` obstruction;
- the exact relative-gauge max-plus convolution for arbitrary two-block
  composition;
- an exact finite-temperature conditional-alignment chain whose
  zero-temperature slope is the optimizer-compatible composition gain;
- the framed order-16 construction and matching lower bound 30 inside its
  pinned oriented-Hadamard family.

The broad archive also records the cube and elliptope relaxations, linear cut
code and signed MacWilliams identity, graphon and spectral losses, negative
replica transport obstruction, scalar microcanonical profile, labeled Fourier
occupancy hierarchy, Paley alignment, Hadamard lifts, and all failed
amplification attempts.

## The active wall

Set

```math
H(n)=F(n)^{2/3}.
```

A uniform estimate

```math
H(n+k)\le H(n)+H(k)+O\!\left((n+k)^{1-\delta}\right)
```

for some `delta > 0` would force convergence. The missing theorem is not a
separate bound on the internal blocks and the rectangular cross block. Such a
bound loses the full leading scale because both contributions can cancel on
the same spin assignment.

The surviving exact state is the labeled relative-switching fiber. For a
deficit threshold `s`, let `b_s(g)` count subthreshold product states in gauge
fiber `g`. Then

```math
b_s(g)=0
```

is exactly the assertion that gauge `g` achieves the desired composition
gain. The next useful result must prove an empty fiber, preferably an
abundance of good fibers, at the power-saving threshold. The complete Fourier
factorization of `b_s` and the mixed four-cycle alignment Hamiltonian are the
current starting points.

Scalar energy profiles, low moments, the weighted Hamming union bound, scalar
negative-replica transport, and graph-plus-rectangular transport have audited
leading-order obstructions. The Paley route still lacks a minimax rigidity
lower bound; on the known dense aligned sequence that lower bound is
equivalent to proving the full limit is `1/2`.

## Repository map

- [Active research program](ACTIVE_RESEARCH.md)
- [Composition-only frontier prompt](COMPOSITION_FRONTIER_PROMPT.md)
- [Finite-results source](paper/finite_results.tex)
- [Composition-framework source](paper/composition_framework.tex)
- [Broad research note](paper/mo-413935-second-attempt.pdf)
- [Independent audit](AUDIT.md)
- [Literature and concept map](LITERATURE.md)
- [Continued proof search](RESEARCH_CONTINUATION.md)
- [Proof-search ledger](RESEARCH_LEDGER.md)
- [Original ledger](STATUS.md)
- [Verification guide](verification/README.md)

## Reproduce the focused checks

Use Python 3.10 or newer. The following checks use exact integer or rational
arithmetic for their pass/fail decisions:

```bash
python3 verification/verify_attempt.py
python3 verification/research_exact_small_n.py --max-n 10
python3 verification/research_order11_certify.py
python3 verification/research_order13_certify.py
python3 verification/research_order15_certify.py
python3 verification/research_cross_block_composition.py
python3 verification/verify_nonlinear_bellman.py
python3 verification/verify_relative_profile_composition.py
python3 verification/verify_labeled_shell_parseval.py
python3 verification/verify_negative_replica_alignment.py
python3 verification/verify_framed_hadamard_lift_30.py
```

Some exhaustive checks require nauty, NetworkX, NumPy, or `z3-solver`; the
order-13 full-stream scan is expensive. The order-15 certificate replays its
deciding tower level in full with
`python3 verification/research_order15_certify.py --full --jobs 4`; the
lower tower levels can be re-derived from the committed driver in
`verification/order15_tower/` in about a day of machine time. The
verification guide gives exact dependencies, expected output, deterministic
seeds, stream digests, corruption controls, and trust boundaries.

The framed order-16 construction also has an independent strict-C verifier:

```bash
cc -std=c11 -O3 -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror \
  verification/verify_framed_hadamard_lift_30.c \
  -o /tmp/verify_framed_hadamard_lift_30
/tmp/verify_framed_hadamard_lift_30
```

Build the two focused manuscripts with:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error \
  -jobname=mo-413935-finite-results \
  -output-directory=paper paper/finite_results.tex
latexmk -pdf -interaction=nonstopmode -halt-on-error \
  -jobname=mo-413935-composition-framework \
  -output-directory=paper paper/composition_framework.tex
```

## Provenance

Rob Sneiderman directs and preserves the project. Proofs, computational
certificates, nonclaims, and failed routes are retained for independent audit.
The search programs were developed with AI assistance and independently
verified as described above; the author is solely responsible for all results
and claims.
