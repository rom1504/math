# Archive verification: coding / association-scheme proposals

Date: 2026-08-16

## Bottom line

The native augmented-cut-code normalization and all three conditional
implications are correct.  At implication level, however, the first two
architectures do not give new reductions:

| Proposal | Class | Decisive reason |
|---|---:|---|
| `L_star` | **C** | Its used content is exactly the archived sharp one-vertex insertion obligation with summable normalized error.  The every-terminal-root quantifier is an unused strengthening, not new leverage. |
| `L_mom` plus growing Terwilliger closure | **C** | The numerical moment assertion is equivalent to the desired sharp universal lower limit by the elementary finite-space `L_p`/`L_infinity` comparison.  The advertised subexponential exact closure is stipulated, not constructed.  Conventional fixed-degree/packing implementations are separately **D**. |
| `L_drift` | **B** | The asymptotic terminal-lumpability assertion is a genuinely different, strict quotient: `(r,b)` provably does not recover the outer histogram.  Its one-edge formula is archived, but no archived theorem contains the uniform large-order drift hypothesis. |

There is no `A` proposal.  `L_drift` merits at most **one tightly bounded,
disproof-first execution checkpoint**: try to amplify the exact terminal
dead ends below into an infinite family, or find two terminal families with
the same limiting `z` and separated `b/N`.  It does not merit an open-ended
attempt to prove (10) before surviving that test.

## 1. Native normalization: verified exactly

Let `N=binom(n,2)`, put `u_e=(1-a_e)/2`, and let

```math
C_n^*=\{(s_i+s_j)_{i<j}:s\in\mathbb F_2^n\},
\qquad C_n^+=C_n^*+\langle\mathbf1\rangle.
```

For `n>=3`, the incidence-map kernel consists of the two constant vertex
words, so `dim C_n^*=n-1`.  The all-one edge word is not a cut (sum the
three equations on a triangle), hence `dim C_n^+=n` and `|C_n^+|=2^n`.

For a switch `x` and augmentation sign `sigma`, the Hamming disagreement is

```math
{1\over2}\left(N-\sigma\sum_{i<j}a_{ij}x_ix_j\right).
```

Therefore, in the report's one-copy convention,

```math
d(u,C_n^+)={N-Q(A)\over2},\qquad
Q(A)=\max_x|H_A(x)|,
```

and consequently

```math
M_n=N-2\rho(C_n^+).
```

The dual is also exactly as required:

```math
(C_n^+)^\perp
=\{F:\deg_F(v)\equiv0\pmod2\ \forall v, |F|\equiv0\pmod2\},
```

the even-cardinality Eulerian sector.  Thus arbitrary signings are arbitrary
ambient roots/cosets.  Codeword-root pair or triple data are packing data and
remain blind to the covering radius.  This agrees with ledger 10.138.1 and
`artifacts/moving_projection_cut_code_coset_audit.md`.

## 2. `L_star`: Class C

### Exact proposed mechanism and implication

For a core `B` and star `b`,

```math
Q(E_b(B))
=\max_y\bigl(|H_B(y)|+|b\cdot y|\bigr).
```

This follows from `max_{t=+-1}|r+ts|=|r|+|s|`; it is a genuinely joint
two-orientation identity.  No positive and negative channel is separately
paid.

The proposed lemma says that for fixed `eta in (0,1/2)`, uniformly for every
large `n` and every terminal-band core `Q(B)<=n^(3/2)`, some star satisfies

```math
Q(E_b(B))
\le\left(1+{3\over2n}\right)Q(B)+C n^{1/2-\eta}.       \tag{S}
```

Applying this only to an order-`n` minimizer gives

```math
M_{n+1}\le\left(1+{3\over2n}\right)M_n+C n^{1/2-\eta}.
```

Since

```math
(n+1)^{3/2}\ge n^{3/2}\left(1+{3\over2n}\right),
```

the normalized positive increment is at most `C n^(-1-eta)`, which is
summable.  A nonnegative sequence with summable positive variation converges.
The implication is correct.

The exact derivative coefficient in the archive is

```math
M_n\left[(1+1/n)^{3/2}-1\right].
```

Its difference from `(3/(2n))M_n` is `O(n^(-1/2))` on the project scale and
is absorbed by `C n^(1/2-eta)` because `eta<1/2`.  Thus the scales align.

### Archive alignment

This is the same implication as:

- ledger 6.2, sharp one-vertex insertion;
- ledger 10.44, equations (10.254)--(10.256), exact insertion discrepancy and
  the summable-error criterion;
- ledger 10.47, the complete variable-radius insertion hierarchy; and
- `artifacts/cap_discrepancy_insertion.md`, Sections 1--2.

The archive only needs the estimate for an exact minimizer.  `L_star` asks it
for every terminal-band core, but that extra quantifier is discarded in the
displayed implication.  At implication level it is therefore an overstrong
restatement of the archived obligation, which is Class C rather than B.

### What is and is not obstructed

The thick-cap example in ledger 10.44 / `cap_discrepancy_insertion.md` proves
`Delta_A(b)>=n-4` for one endpoint-balancing row on a noncompetitive
quadratic-cap core.  It does **not** prove that every row fails, and it does
not give a competitive infinite family with

```math
\min_b\Delta_A(b)-{3M(A)\over2n}\ge c\sqrt n.
```

Accordingly, it is not a D-level falsifier of (S).  Likewise, the order-six
insertion jump is only finite normalization evidence.

The path-versus-matching example in the proposal correctly proves that the
parent endpoint `Q(B)` does not determine the outer histogram.  But (S)
still takes the complete core `B` as input and contains a maximum over all
`2^(n-1)` antipodal spins.  It supplies no closed compressed state.  The
rooted-transfer audit shows that exact iteration through all future stars
generates the multi-overlap hierarchy; the natural external-field support
function is injective and recovers the full energy word.  This does not
falsify a one-step lossy proof of (S), but it prevents treating the `n`-bit
witness row as an established state compression.

### Required audit flags

- Arbitrary roots: yes, stronger than needed; the implication uses an exact
  minimizer at every large order.
- Packing versus covering: correctly covering-rooted; no packing import is
  used as a substitute.
- Channels: joint until the final maximum; no separate-channel loss.
- Leading coefficient/error: correct `3/2`, power-saving and summable.
- Conference-only: no; the conference frontier is used only to place
  minimizers inside the terminal band.
- State size: the constraint family is exponential and no compression theorem
  is supplied.  Exact repeated closure is full-information.
- Surviving issue: prove or scalably falsify the archived weighted thick-cap
  discrepancy estimate for optimal cores.  The uniform-all-root version is
  not the minimal theorem.

## 3. `L_mom`: Class C overall; familiar implementations D

### Exact implication, and converse

The proposal defines

```math
S_{n,k}(A)=2^{-n}\sum_x|H_A(x)|^{2k}
```

and asks, for every `epsilon>0`, for some `alpha(epsilon)>0`,
`k=ceil(alpha n)`, and every signing `A`,

```math
S_{n,k}(A)^{1/(2k)}
\ge(1/2-\epsilon)n^{3/2}.                            \tag{M}
```

The forward implication is correct: `Q(A)>=||H_A||_(2k)`, so (M), minimized
over `A`, and the all-order conference upper frontier give convergence to
`1/2`.

But (M) is also a consequence of that target.  On a space of `2^n` points,

```math
S_{n,k}(A)^{1/(2k)}
\ge2^{-n/(2k)}Q(A)
\ge2^{-n/(2k)}M_n.                                  \tag{M-conv}
```

If `M_n/n^(3/2)->1/2`, choose `alpha` sufficiently large that
`2^{-1/(2alpha)}` loses less than the prescribed `epsilon`, then take `n`
large.  This proves (M).  Hence the numerical content of `L_mom` is
implication-level equivalent to the sharp universal lower limit, the archived
ledger 6.3 obligation.  Abstract noninjectivity of one power sum does not
establish a strict reduction when the requested lower bound is equivalent to
the answer by (M-conv).

### Growing-state clause and archive collision

The moment expansion as a signed ordered even-degree Eulerian multigraph sum
is correct.  The required scale is also correct: the radius displacement is
`Theta(n^(3/2))` in ambient length `Theta(n^2)`, so generic moment/Krawtchouk
degree must be `Theta(n)`, not fixed.

What is not supplied is the proposed sign-uniform PSD closure with
`exp(o(n))` total block state and `o(n^(3/2))` root error.  It is a requirement
that a proof exist, not a construction of such a state.  The aligned archive
results are:

- ledger 3.15 and 9.2: signed Eulerian/Krawtchouk data and exact moment
  convolution; linear moments retain cancellation and entropy losses;
- ledger 10.138.1--10.138.3: packing is coset-blind, the moving-Gram root is a
  weighted full coset enumerator, and correct scale needs degree `Theta(n)`;
- ledger 10.139: operator rank does not remove the root obligation; direct
  matching support and pure add/delete matching transitions have scalable
  coefficient no-gos;
- ledger 10.140.2 and `artifacts/eulerian_free_energy_identity.md`: fixed
  cycle/replica depth misses zero-entropy resonances and vertex deletion opens
  exponentially many boundary sectors; growing to the required degree can
  restore the complete signed Eulerian/coset histogram; and
- `artifacts/rooted_switching_transfer_operator.md`: exact growing rooted
  closure reaches the full labeled energy word and is projectively
  noncontractive.

Thus:

- a Schrijver/codeword-root or moving-projection packing implementation is D
  by coset blindness;
- fixed-degree dual-cycle/Terwilliger truncation is D by the explicit
  resonance/boundary-sector obstruction;
- scalar partial-transversal support `exp(O(n))` is D at the correct scale,
  which requires hidden support `exp(Omega(n log n))`; and
- the generic degree-`Theta(n)` all-root proposal is C, because its numerical
  lemma is the desired lower theorem in `L_p` clothing and its claimed
  subexponential closure has not been exhibited.

No proved no-go covers every conceivable noncommuting graph-orbit algebra.
If an actual algebraically closed, subexponential, arbitrary-root hierarchy is
constructed and is proved nonrecovering, it should be reclassified on that
specific construction; merely including that outcome in the lemma does not
make the present proposal B.

### Required audit flags

- Arbitrary roots: correctly every signing in (M); imported packing theorems
  do not meet it.
- Full histogram: one numerical moment is noninjective, but the sharp bound is
  target-equivalent; exact degree-`Theta(n)` generic closure can recover the
  signed coset histogram.
- Channels/cancellation: the final even moment is nonnegative, but its
  Eulerian expansion is signed.  Dropping negative diagrams or paying graph
  types separately is invalid.
- Leading loss: none in the statement; an `exp(Theta(n))` multiplicative
  uncertainty becomes a fixed loss after the `Theta(n)`-th root.
- Conference-only: conference signings supply only the legitimate all-order
  upper frontier, not the universal lower theorem.
- State size: `exp(o(n))` is stipulated.  Known exact states are exponential
  at the required accuracy; the moving-kernel support lower bound is stronger,
  `exp(Omega(n log n))`, for the scalar class.
- Surviving issue: none as a strict reduction.  A concrete new closed algebra,
  rather than the meta-level closure clause, would be needed for reassessment.

## 4. `L_drift`: Class B

### Exact implication

For the quotient Cayley graph generated by edge coordinates, put

```math
r(U)=d(U,C_n^+),\qquad
z(U)={N-2r(U)\over n^{3/2}},
```

```math
b(U)=|\{e:r(U+e)=r(U)+1\}|.
```

The proposed `L_drift` asks that, uniformly for every sufficiently large `n`
and every coset with `z in I superset [0.33,0.51]`,

```math
\left|{b(U)\over N}-\beta(z(U))\right|\le\epsilon_n,
\qquad\epsilon_n\to0,
```

where continuous `beta` has one zero `c`.

For a deepest coset, `b=0`: no coordinate neighbor can have distance larger
than the covering radius.  The known frontier puts its `z` in the mandatory
interval for all sufficiently large orders.  Compactness and the unique-zero
condition then force every subsequential limit to be `c`.  The conditional
convergence implication is correct.  Same-layer edges in the quotient do not
affect this argument.

### Exact archived one-edge formula

Represent `U` by a signing `a`, put `Q=Q(a)`, and for augmented cuts `v` set

```math
g_v={Q-a\cdot v\over2},\qquad
N_v=\{e:a_ev_e=-1\}.
```

Ledger 3.21 proves, for an edge set `S`,

```math
Q(a^S)=Q-2\min_v\{g_v+|S|-2|S\cap N_v|\}.
```

For `S={e}`, this gives the exact identity

```math
\boxed{
b(U)=N-\left|\bigcup_{g_v\le1}N_v\right|.}          \tag{D1}
```

Indeed, an edge is outward precisely when it is uncovered by the top two
energy layers.  Thus `b` is not a generic association-scheme intersection
number: it is a signing-rooted near-ground incidence statistic already
implicit in the archive.  What is new is the assertion that this statistic,
uniformly over all terminal roots, asymptotically becomes a continuous
function of the single scalar `z`.

### It is genuinely a strict quotient

Unlike `L_mom`, strictness can be proved on actual augmented-cut-code cosets.
Gauge the order-six root so all edges incident with vertex `0` are positive.
Consider internal negative-edge sets

```math
F_1=\{12,13,23\},
\qquad
F_2=\{13,14,23,24\}.
```

Exact enumeration gives both roots

```math
r=3,\qquad Q=9,\qquad b=0,
```

but their full outer distributions differ:

```text
F1: A_3=2, A_4=6, A_7=24, A_8=24, A_11=6, A_12=2.
F2: A_3=1, A_4=3, A_5=6, A_6=10, A_7=12, A_8=12,
    A_9=10, A_10=6, A_11=3, A_12=1.
```

Each list sums to `2^6`.  Hence `(r,b)` is noninjective even on the actual
coset family and cannot recover the full outer histogram.  This establishes
the proposal's positive strictness claim, not merely a description-length
claim.

The caveat is computational rather than informational: evaluating `r` and
the top-two-layer union in (D1) still contains the original Boolean maximum
at the given root.  A scalar output does not by itself give a tractable proof
of its uniform law.

### Exact finite coset-graph census through order eight

The temporary exact checker gauge-fixes the vertex-0 star, enumerates the
`2^(N-n)` residual patterns modulo global complement, evaluates every
one-copy cap over `2^(n-1)` spins, and tests all `N` edge generators.  All
arithmetic is integral.  The source is
`tmp/retrieval_panel_2026_08/verifications/coset_drift_enum.cpp`.

| `n` | cosets | `rho`, `M_n` | exact nondeep dead ends | same-layer outward-degree variation |
|---:|---:|---:|---:|---|
| 4 | 4 | `1,4` | 0 | none |
| 5 | 32 | `3,4` | 0 | none |
| 6 | 512 | `5,5` | 25 | at `Q=9`, `b in {0,4}`; also variation at `Q=11` |
| 7 | 16,384 | `6,9` | 1,260 | at `Q=11`, `b in {0,2,4,5,10}` |
| 8 | 1,048,576 | `9,10` | 130,655 | at `Q=12`, `b in {0,1,3}`; large variation on higher layers |

The nearest nondeep-dead-end layers have normalized defects

```text
n=6: Q=9,  z=0.612372...
n=7: Q=11, z=0.593944...
n=8: Q=12, z=0.530330...
```

so they lie just above, not inside, the proposal's mandatory upper endpoint
`0.51`.  Within `[0.33,0.51]`, the exhaustive results are:

```text
n=4: Q=4, b=0 (deepest);
n=5: Q=4, b=0 (deepest);
n=6: Q=7, b=1 on all 90 cosets; Q=5, b=0 (deepest);
n=7: Q=9, b=0 (deepest);
n=8: Q=10, b=0 (deepest).
```

Therefore the `n<=8` census refutes exact complete regularity and displays
many local traps, but it does **not** falsify the asymptotic terminal-band
statement.

There is, however, a sharper exact order-ten warning.  Ledger 10.1401
certifies `M_10=13`.  The signing whose negative edges are

```text
01 02 03 04 05 06 14 15 16 17 18 19 24 25 27 38
47 49 56 58 59 78 89
```

has exactly `Q=15`, and every one of its 45 single-edge flips has `Q=17`.
Thus

```math
r=15<\rho=16,\qquad b=0,
```

while both `15/(10 sqrt(10))=0.474341...` and the deepest value
`13/(10 sqrt(10))=0.411096...` lie in `[0.33,0.51]`.  This is an exact
finite terminal-band nondeep dead end.  It is not a logical falsifier because
`L_drift` starts only at an unspecified sufficiently large order, but it
shows that uniform terminal lumpability cannot be justified by finite
complete-regularity intuition.  The fixed-seed checker is
`tmp/retrieval_panel_2026_08/verifications/coset_drift_sample.cpp`.

### Archive/no-go boundary and minimum revision

- Complete regularity would imply the desired dependence on distance, but no
  complete-regularity hypothesis is known.  The proposal correctly proves
  that `S_n` symmetry cannot imply complete transitivity; this does not rule
  out approximate terminal lumpability.
- Ledger 3.21 and `artifacts/eulerian_free_energy_identity.md` show exact
  nonglobal one-edge traps.  They do not supply an infinite terminal-band
  family, so they do not contain all hypotheses needed for Class D.
- The state does not pay orientation channels separately, has no recurrence
  coefficient loss, uses every ambient coset rather than a packing root, and
  is not conference-only.
- The output state `(z,b/N)` is constant-size and strictly nonrecovering, but
  its evaluation is not known easier than `N+1` rooted optimizations.  No
  cross-order closure theorem is supplied.

Minimum revision for this B proposal:

1. replace the complete-regularity analogy by the exact top-two-layer formula
   (D1);
2. state explicitly that order 10 already has a mandatory-band nondeep dead
   end; and
3. make the first checkpoint a scalable-trap test.  Either amplify such dead
   ends to infinitely many orders (which kills `L_drift`) or prove a theorem
   excluding them uniformly in a fixed terminal neighborhood.  Only after
   that should one attempt existence/continuity/unique-zero of `beta`.

## 5. Director-facing recommendation

Reject `L_star` and `L_mom` as execution candidates in their current form:
their implications respectively recover the archived insertion obligation
and the desired sharp lower theorem.  Do not label either D wholesale:
`L_star` has no scalable all-row falsifier, and no theorem rules out every
possible noncommuting growing algebra for `L_mom`.

Retain `L_drift` as B only because its strict quotient is now proved on actual
cosets and no archive result contains its asymptotic lumpability hypothesis.
If the panel has a stronger A/B candidate, prefer it.  Otherwise authorize
one bounded disproof-first campaign, with the exact stop condition:

```math
\exists U_n,V_n:\ z(U_n)-z(V_n)\to0,
\quad |b(U_n)-b(V_n)|/N_n\not\to0,
```

or an infinite family of terminal nondeep dead ends separated from the
deepest value.  Either result closes the route.  Failure to obtain a scalable
construction after one substantive checkpoint should return it to hold,
rather than opening adjacent complete-regularity variants.
