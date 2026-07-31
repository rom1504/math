# Phase 2 falsification/computation report

This report consolidates the second sustained phase's structured-landing
computations.  It distinguishes exhaustive certificates from deterministic
samples and heuristic nonlinear optimization.  None of the finite results
below proves or disproves convergence.

## Reproduction

All temporary binaries were placed under the repository-local temporary
directory, never the system `/tmp` directory.

```bash
g++ -O3 -std=c++17 -Wall -Wextra -pedantic \
  computations/phase2_subset_caps_gray.cpp \
  -o /home/math/quadra/tmp/phase2_subset_caps_gray

.venv/bin/python computations/phase2_paley_orbit_landing.py \
  --p 5 \
  --evaluator /home/math/quadra/tmp/phase2_subset_caps_gray \
  --output computations/results/phase2_paley_orbits_p5.json

.venv/bin/python computations/phase2_restriction_state_audit.py \
  computations/results/phase2_paley_orbits_p5.json \
  --p 5 --max-restriction-size 7 \
  --output computations/results/phase2_pc26_restriction_state_audit.json
```

The smaller landing samples were produced with commands of the form

```bash
.venv/bin/python computations/phase2_structured_landing_audit.py \
  SOURCE.json --matrix-key MATRIX_KEY --label LABEL \
  --child-orders K --samples 10000 --seed 20260731 \
  --evaluator /home/math/quadra/tmp/phase2_subset_caps_gray \
  --output computations/results/phase2_landing_LABEL.json
```

The fractional diagnostics were produced with commands of the form

```bash
.venv/bin/python computations/phase2_fractional_bridge_variance.py \
  computations/results/exact_mN.json computations/results/exact_mN.json \
  --target-slack S --restarts 6 --iterations 4 --seed SEED \
  --output computations/results/phase2_fractional_variance_nN_slackS.json
```

## Orbit-complete balanced restrictions of `PC(26)`

**Classification: exhaustive and certified computationally.**

The script constructs `PGammaL(2,25)` from translations, nonsquare scaling,
inversion, and Frobenius.  It verifies directly on the saved integer
conference matrix that every generator acts by switching and possible global
negation.  These operations preserve Boolean cap.  The generated permutation
group has the checked order

```math
2q(q^2-1)=2\cdot25\cdot(25^2-1)=31,200.
```

Its action partitions all
`binom(26,13)=10,400,600` balanced subsets into 391 orbits.  One exact
Gray-code evaluation per orbit gives the complete distribution

| child cap | subset orbits | labeled subsets |
|---:|---:|---:|
| 24 | 92 | 2,610,400 |
| 26 | 177 | 4,913,300 |
| 28 | 113 | 2,750,800 |
| 30 | 9 | 126,100 |

Consequently every balanced principal restriction of `PC(26)` has cap at
least 24.  Since `M_13=20`, this particular structured family has the exact
finite landing gap

```math
24^{2/3}-20^{2/3}=0.9522722949\ldots,
\qquad
\frac{24^{2/3}-20^{2/3}}{13}=0.0732517150\ldots.
```

This falsifies exact balanced landing for `PC(26)`.  It does **not** prove a
linear gap, falsify little-oh landing, or rule out another structured family.

There is also a scalable limitation on this particular symmetry compression.
The number of balanced-subset orbits is at least

```math
\frac{\binom{q+1}{(q+1)/2}}{|PGammaL(2,q)|}.
```

For `q=25`, the denominator used here is the actual checked order 31,200, so
the lower bound is 334 orbits and the exact count is 391.  Along square fields
`q=p^2`, the group order is `2q(q^2-1)`, polynomial in `q`, whereas the central
binomial coefficient is exponential.  Thus an automorphism-orbit label alone
is not a bounded-complexity state.  This statement does not rule out a much
coarser invariant that controls the cap across many orbits.

## Spectral-state falsifier

**Classification: exhaustive and certified computationally.**

For every one of the 391 representatives, exact power traces and the monic
characteristic polynomial were computed with Newton identities and checked by
exact Cayley--Hamilton substitution.  Modulo global negation, there are 187
complete spectral states.  Eleven spectral states each contain both a cap-24
orbit and a cap-28 orbit.  Hence even the complete child spectrum does not
determine Boolean cap inside this single balanced Paley restriction family.

This directly falsifies a structured state consisting only of conference
completion data, spectrum, and Paley-order information.  Any useful state
needs genuinely Boolean information, not just additional spectral moments.

## Bounded Boolean restriction profiles

**Classification: exhaustive and certified computationally at order 13;
uniform generalization open.**

For each orbit representative, the follow-up script classifies every
principal restriction under switching, permutation, and global negation and
records the complete class histogram.  The numbers of universal signing
classes are

| restriction order | classes |
|---:|---:|
| 4 | 2 |
| 5 | 4 |
| 6 | 10 |
| 7 | 27 |

The exact collision audit gives:

| retained state | distinct states on 391 orbits | ambiguous cap states |
|---|---:|---:|
| full spectrum | 187 | 11 (`24,28`) |
| spectrum + complete profiles through 4 | 187 | 11 (`24,28`) |
| spectrum + complete profiles through 5 | 334 | 1 (`24,28`) |
| spectrum + complete profiles through 6 | 376 | 0 |
| spectrum + complete profiles through 7 | 389 | 0 |
| profiles through 4, no spectrum | 20 | 13 |
| profiles through 5, no spectrum | 258 | 14 (`24,28`) |
| profiles through 6, no spectrum | 376 | 0 |
| profiles through 7, no spectrum | 389 | 0 |

Thus profiles through five vertices are rigorously insufficient at this
order.  Profiles through six happen to determine the cap on the complete
`PC(26)` balanced family, but they distinguish 376 of 391 orbits and therefore
give almost no finite compression.  This is a reproducible candidate, not a
theorem: the next falsifiable question is whether the 16-component vector of
class counts through order six controls cap, even approximately, on held-out
Paley orders or nonconference signings.  A single pair with the same vector
and different cap falsifies exact control.

A bounded held-out computation exhausts the complete nonconference universe
at order 8.  The program
`computations/phase2_profile_collision_n8.cpp` checks all `2^20=1,048,576`
root-gauged signings modulo global negation, computes every cap exactly, and
finds only 131 size-4-through-6 profile states.  No state contains two
different caps.  This is **exhaustive finite positive evidence**, saved in
`computations/results/phase2_profile_collision_n8.json`; it is not a uniform
theorem.

One bounded held-out search then tested order 9.  A deterministic affine
permutation sampled 1,000,000 of the `2^27` global-negation representatives,
finding 1,000 distinct profiles and no cap collision.  This result is
**sampled, not exhaustive** and is saved as
`computations/results/phase2_profile_collision_n9_sample1m.json`.  The search
was stopped rather than escalated to order 10.  Absence of a collision in this
batch is positive finite evidence only.

The vector has fixed dimension 16, but naive state enumeration is still not a
practical composition method.  A coarse counting bound is polynomial with a
very high exponent:

```math
\#\text{profiles}
\le\prod_{t=4}^6\left(\binom nt+1\right)^{c_t}
=n^{O(88)},
\qquad(c_4,c_5,c_6)=(2,4,10).
```

At the only larger complete test it distinguishes 376 of 391 available
orbits, and no rule is known for selecting a low-cap state or updating the
state under a bridge.  The observed state counts are 131 on the entire order-8
universe, 1,000 in the million-signing order-9 sample, and 376 on the special
391-orbit order-13 family; they do not yet define a stable scaling law.
Moreover, a proof using only the spectral moments
encoded by restrictions through six cannot reach the target scale: the sixth
moment estimate `||A||<=tr(A^6)^(1/6)` is typically only of order `n^(2/3)`,
giving the Boolean spectral bound `cap(A)<=n^(5/3)/2`, rather than the needed
`O(n^(3/2))`.  Any positive use of the profile therefore needs a genuinely
Boolean inequality or a growing-order statistic, not a fixed-moment spectral
argument.

The main track's separate switching-rank-one cell theorem also rules out a
bounded number of coarse rank-one bridge cells at the required scale.  The
restriction-profile candidate is not such a cell decomposition, so it is not
automatically killed by that theorem.

## Sampled landing data

**Classification: every reported child cap is exact; a minimum is certified
only in rows marked exhaustive.**

| structured parent | child order | subsets | best cap | `M_k` | best cap gap |
|---|---:|---:|---:|---:|---:|
| `PC(10)` | 5 | all 252 | 4 | 4 | 0 |
| `PC(14)` | 7 | all 3,432 | 9 | 9 | 0 |
| `PC(18)` | 9 | deterministic 10,000 | 12 | 12 | 0 |
| `PC(26)` | 13 | deterministic 10,000, held-out seed | 24 | 20 | 4 |
| Hadamard-derived double, order 12 | 6 | all 924 | 5 | 5 | 0 |
| Hadamard-derived double, order 28 | 14 | deterministic 10,000 | 25 | 21 | 4 |

The sampled rows only upper-bound the unknown structured minimum.  The
`PC(26)` row is superseded by the orbit-complete certificate above.  The table
shows that exact landing at orders 5, 7, and 9 does not extrapolate to order
13; it is not by itself an asymptotic trend.

## Fractional-bridge variance diagnostic

**Classification: exact linear constraints and numerically checked feasible
points; heuristic, non-global minimization of variance.**

For fixed children `A,B`, the fractional bridge polytope is exactly

```math
-T\le H_A(x)+H_B(y)+x^{\mathsf T}Cy\le T
\quad\text{for all projective }(x,y),\qquad -1\le C_{ij}\le1.
```

At the ideal power target
`T=(M_A^(2/3)+M_B^(2/3))^(3/2)`, randomized independent rounding has variance
parameter `V=sum(1-C_ij^2)`.  The script searches exposed LP vertices and uses
monotone tangent ascent to heuristically reduce `V`.  For equal saved children
at orders 5, 6, and 7, the best observed `V/n^2` at cap slacks `0,1,2,4` was:

| child order | slack 0 | slack 1 | slack 2 | slack 4 |
|---:|---:|---:|---:|---:|
| 5 | 0.4251 | 0.1970 | 0.1835 | 0 |
| 6 | 0.4440 | 0.4018 | 0.1661 | 0.1648 |
| 7 | 0.08895 | 0.05698 | 0.03658 | 0.03658 |

These values are upper bounds on the minimum possible variance, not lower
bounds.  The method sometimes misses known integral-feasible vertices, and
deeper runs remained erratic.  Order-8 class-pair runs were stopped after the
LP cost increased without a checkable trend.  The diagnostic therefore does
not support a scaling claim.

It does isolate a precise obligation for this route.  With a union-bound
factor `L=Theta(n)`, the rounding addition is of order `sqrt(VL)+L`; obtaining
an `o(n^(3/2))` cap error requires at least `V=o(n^2)`, and a summable
`b=M^(2/3)` defect requires correspondingly quantitative decay.  No such
uniform variance estimate is currently known.

## Track judgment

The primary result of this track is negative but scalable in architectural
scope: automorphism orbits remain exponentially numerous and complete spectra
do not control Boolean cap.  The certified `PC(26)` landing gap prevents the
earlier sample from being dismissed as unlucky, but does not establish a
linear obstruction.

The only surviving finite compression candidate found here is the complete
histogram of switching classes on restrictions through six vertices.  It
determines cap both on all order-8 signings and on the balanced `PC(26)`
family, but is nearly injective on the latter orbit set and has no bridge
update theorem.  It should be retained only as a sharply falsifiable invariant,
not as evidence that structured landing has been proved or made substantially
easier.
