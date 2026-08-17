# Independent audit: mesoscopic-core iteration frontier

**Verdict: PASS after four statement/scope repairs.**

This audit freezes `mesoscopic_core_iteration_frontier.md` at SHA-256
`5852ee6605ab0424bd88aca59e19ac60a8e50f4bc0c1221fc563592c20e0ec2c`.
The multi-anchor inequality is correct, and the proposed iteration really
does avoid cumulative shell-width loss because every witness is certified by
a fresh perturbation of the same exact minimizer.  The draft does not prove
reservoir persistence, a growing packing, or a fixed-edge-scale packing.

Two hypotheses are missing from the general theorem, the PP.4 initialization
must exclude `beta=1`, and the finite-ball example should be stated more
precisely.  The three-cut construction is valid as a cut-realizable Helly
falsifier and its `O(n^(3/2))` cap follows by a standard conditioned-sign
probabilistic argument, but it is not an exact-minimizer or thin-shell
counterexample.

## 1. Multi-anchor conventions and exact-flip algebra: PASS

On the common-correct reservoir

```math
R=\{e:a_e(z_i)_e=+1\ \forall i\},
```

every old anchor satisfies `(z_i)_e=a_e`.  If a new response has
`a_ez_e=-1`, it differs from **every** old anchor at that coordinate.  Thus
one hypergeometric exclusion on `R`, rather than one separately paid
exclusion per anchor, genuinely gives simultaneous actual-distance control.

For an `r`-set `F subseteq R`, let `q` be the number of negative
correlations `a_ez_e=-1` on `F`.  Exact minimality gives

```math
M\le\langle a^F,z\rangle
=M-d-2r+4q,
```

so

```math
q\ge r/2+d/4,\qquad 0\le d\le2r.
```

This agrees exactly with the archived finite-flip certificate underlying
FB.3.  The hypergeometric bound

```math
\Pr\{q_F(w)\ge r/2\}
\le\exp[-2(1/2-\theta)^2r]
```

is correct for every response whose negative-correlation population on `R`
is at most `theta p_k`.  Union over at most `2^n` augmented cuts handles the
adaptive maximizer.  There is no extra factor depending on the number of
anchors.

## 2. Projective-distance constants: PASS

Let `h_i=d_E(z,z_i)`.  The common-reservoir conclusion gives
`h_i>theta p_k` for every `i`.  Coordinatewise energy addition gives

```math
\langle a,z\rangle+\langle a,z_i\rangle
\le2(E-h_i),
```

and hence

```math
E-h_i\ge M-s-d/2\ge M-s-r.
```

Therefore

```math
d_P(z,z_i)\ge\min\{\theta p_k,M-s-r\}.
```

All factors of two are correct.  The complementary estimate is precisely
AO.20 with the two unequal shell deficits retained.  It is not a new way
around the absolute-overlap ceiling.

## 3. Required repair 1: positivity is missing in MI.1

MI.1 calls the selected response a positive augmented cut but assumes only
`d<=2r`.  As in the first MB.2 draft,

```math
\langle a,z\rangle=M-d\ge M-2r
```

need not be positive when `2r>=M`.  Add

```math
2r<M.
```

Alternatively, distinguish the selected response orientation from the
positive representative of its projective line; the latter remains in the
deficit-`2r` shell, but the signed common-reservoir conclusion belongs to the
former orientation.  For iterative use the clean `2r<M` hypothesis is
preferable.  It holds eventually in MI.2 because `r=O(n)=o(M)`.

## 4. Required repair 2: “adds one state” needs a positive gap

The displayed inequality remains formally true if `M-s-r<=0`, but then it
does not prove that `z_(k+1)` is projectively distinct from the old anchors.
Thus either retitle MI.1 as only a distance inequality or assume

```math
s+r<M.
```

With `p_k>=1` this makes the right side of (MI.4) positive and justifies
“adds one jointly separated state.”  Again MI.2 satisfies this eventually.

## 5. No cumulative loss and the iteration dichotomy: PASS

Put `s_*=max{s,r}`.  Every old anchor has deficit at most `2s_*`, while
each newly certified anchor has deficit at most `2r<=2s_*`.  At the next
step MI.1 is applied to the same base signing `A`, not to the previously
flipped signing.  Therefore neither the base cap, the deficit, nor previous
pairwise distances are modified.  There is genuinely no additive `kr` loss.

If `p_k>=beta M`, then with `theta=1/4`,
`r=ceil(9n log2)`, and `s_*=o(M)`,

```math
d_P(z_(k+1),z_i)
\ge\min\{\beta M/4,M-s_*-r\}
=(\beta/4-o(1))M
```

for fixed `0<beta<=1`.  Induction preserves the existing packing.  Either it
reaches `K` members or the exact reservoir test fails first.  This dichotomy
is logically complete and remains valid even for growing `K`, conditional on
reservoir persistence.

The phrase “exact state variable” should mean “the exact state used by this
proof,” not “a proved minimal or necessary state.”  MI.2 supplies a useful
conditional iteration theorem, but its packing-or-failure alternative is an
induction once MI.1 is known.

## 6. Required repair 3: PP.4 does not initialize `beta=1`

For the PP.4 pair the available estimate is

```math
p_2=|Z|\ge M-2s.
```

This exceeds `beta M` eventually for every fixed `beta<1`, but need not do
so for `beta=1`.  Therefore the sentence saying MI.2 recovers MB.2 from the
PP.4 pair must either assume `0<beta<1`, use a threshold such as
`beta(M-2s_*)`, or invoke MB.2 directly.  The endpoint `beta=1` is valid in
the abstract dichotomy only when its stronger reservoir premise is actually
given.

With any fixed `beta<1`, MI.2 recovers an `Omega(M)`-separated third state.
It does not recover MB.2's sharper `(1/2-o(1))M` constant unless one also
uses its varying `theta` and `O(n log^2 n)` choice of `r`.

## 7. The three-cut Helly obstruction: PASS with exact scope

For disjoint `S,T` of size `k=Theta(sqrt n)`, the four cells have sign
patterns

```text
A_0:(+,+,+), B_0:(+,+,-), C_0:(+,-,+), D_0:(+,-,-),
```

for `(z_0,z_1,z_2)`.  Cell sums `(0,m,m,-m)`, up to parity-one errors,
give all three displayed energies `m+O(1)`.  They are feasible after choosing
the constant in `m=Theta(n^(3/2))` below the sizes of `B_0` and `C_0`.

The common-correct reservoirs have the claimed geometry:

* `R(z_0,z_1)` contains the positive part of `B_0`, of size
  `(|B_0|+m)/2=Theta(n^(3/2))`;
* `R(z_0,z_2)` similarly contains the positive part of `C_0`;
* `R(z_1,z_2)` contains the negative part of `D_0` and is in fact
  `Theta(n^2)`;
* the triple reservoir is exactly the positive part of `A_0`, hence has
  size at most `|A_0|=k^2=Theta(n)`.

Thus there is no purely set-theoretic pairwise-to-joint Helly implication,
even inside the augmented-cut response family.

The bounded-cap existence claim is also valid, but the draft should spell
out its probabilistic status.  Independently choose the signs in each cell
uniformly subject to the prescribed cell sum.  For a fixed augmented cut,
the conditional expectation of its contribution from cell `C` has magnitude
at most the absolute prescribed cell sum.  Hoeffding for sampling without
replacement makes the sum of the four centred contributions subgaussian
with variance proxy `O(E)`.  A union bound over at most `2^n` cuts gives a
simultaneous fluctuation `O(sqrt(En))=O(n^(3/2))`.  Its mean is at most
`3m+O(1)`, so some exact `+-1` signing has `Q(A)=O(n^(3/2))` while retaining
all four cell sums.

This construction is **not** an exact-minimizer example, and the three
energies are not proved to be within `o(n^(3/2))` of `Q(A)`.  It therefore
does not falsify a higher-order linkage theorem specifically for thin-shell
states of exact minimizers, nor does it falsify adaptive witness selection
in (MI.14).  It proves only the narrower and useful statement in the draft:
cut algebra, bounded cap, and large pairwise reservoirs alone do not force a
large joint reservoir.  Replace “exact signing” by “exact `+-1` signing” to
avoid confusion with “exact minimizer.”

## 8. Required repair 4: make the finite-ball example literal

The conclusion of Section 3 is right: an `o(E)` projective ball can contain
a large packing at edge scale `M=Theta(n^(3/2))`.  Its current
`n^(3/4)`-mask wording does not literally make all constant-rate-code
distances `Theta(n^(3/2))`; constant relative mask distance there produces
`Theta(n^(7/4))` cut distance.

A cleaner exact construction is to reserve `L=Theta(sqrt n)` vertices and
take a constant-rate binary code of length `L` and fixed relative distance.
Every code mask lies in one projective cut ball of radius

```math
L(n-L)=Theta(n^(3/2))=o(E),
```

and every pair has projective cut distance `Theta(n^(3/2))`.  The ball
contains `exp(Theta(sqrt n))` such points.  This proves the intended claim
more strongly and without an ambiguity.  If the original larger ball is
retained, say “distance at least `Theta(n^(3/2))`,” not “distance
`Theta(n^(3/2))`.”

The subsequent observation about radius `o(M)` is correct by the triangle
inequality.  However, distance from all prior representatives is the output
of MI.1 **under** its reservoir hypothesis, not logically equivalent to that
hypothesis; other mechanisms could produce it.

## 9. FB.3/AO.20 comparison and theorem-level judgment

The projective volume estimate is correct.  A radius-`Theta(M)` cut ball
corresponds to vertex masks of size `Theta(sqrt n)` and has

```math
\exp(Theta(\sqrt n\log n))
```

words.  FB.3's guaranteed logarithmic shell cardinality is `Omega(r/n)`;
for every `r=o(n^(3/2))` this lower bound is `o(sqrt n)` and cannot beat that
ball volume.  This says only that **FB.3's bound** is insufficient, not that
the actual shell lacks entropy.

The mathematical classification is:

| Component | Classification |
|---|---|
| MI.1 | theorem-level conditional multi-anchor consequence of the archived exact-flip certificate |
| MI.2 | correct induction/no-accumulation corollary, conditional on a new higher-order reservoir premise |
| (MI.6) | exactly AO.20, not a new projective-scale mechanism |
| three-cut construction | scalable cut-realizable falsifier of pairwise Helly, not an exact-minimizer falsifier |
| fixed-scale Cut-DH progress | none |
| new unconditional near-minimizer structure | none |

Thus the draft is more than a vocabulary-only reformulation: MI.1 proves a
simultaneous multi-anchor separation statement, and the three-cut example
decisively rejects a naive pairwise-Helly closure.  But it introduces no new
optimality engine beyond FB.3's finite-flip certificate and AO.20, and MI.2
does not by itself reduce the exact-minimizer frontier.  The genuinely new
open content remains a thin-shell-specific higher-order reservoir or
non-recycling theorem.
