# Cap and local-field comparison of all cyclic order-38 orbits

Date: 2026-07-31. This is an agent-authored computational report. Exact
finite certificates, heuristic search observations, and scaling hypotheses
are kept separate.

## 1. Exact structural input

The exhaustive cyclic search at `k=3`, `s=19` contains 627 oriented
self-indexed ASDS pairs and exactly three orbits under simultaneous unit
multipliers, bridge translations, and bridge reversal. Representatives have
normalized conference-graph clique invariants

```text
natural orbit size   K4    K5    class
171                  555    0    Paley
342                  615   65    non-Paley
114                  610   48    non-Paley
```

For every representative, the comparison program independently rechecks
the sequence sums, complementary autocorrelation identity, self-indexed ASDS
identity, signing alphabet, and

```math
S^2=37I_{38}.                                             \tag{LC1}
```

Thus the two non-Paley examples are genuinely distinct structured signings,
not alternate gauges of the prime Paley conference. This is an exact finite
classification under the stated natural operations. It has no asymptotic
content by itself.

## 2. Reproducible cap search

For each of the three representatives, the program runs, for both energy
signs:

- 200,000 independent random starts followed by exact steepest single-spin
  ascent;
- 5,000 independent random starts followed by exact best-improvement ascent
  over every one- and two-spin flip.

This gives 400,000 single-flip terminal states and 10,000 pair-flip terminal
states per orbit. Every retained witness energy and local-field penalty
identity is recomputed with exact integers. The result is the same in all
three graph classes:

```math
\operatorname{cap}(S)\ge109.                            \tag{LC2}
```

The universal conference arithmetic theorem gives

```math
\operatorname{cap}(S)\le113,                            \tag{LC3}
```

so the current rigorous interval for each fixed representative is
`[109,113]`. Neither the large heuristic search nor separate CP-SAT decisions
found energy 111 or 113, but the decision models returned `UNKNOWN`.
Therefore **109 is not claimed as the exact cap**.

The sampled single-flip terminal-energy distributions are:

```text
Paley (555,0):
  85: 3,220; 93: 106,311; 101: 205,523; 109: 84,946

non-Paley (615,65):
  77: 4; 85: 1,727; 93: 62,814; 99: 1,596;
  101: 291,528; 107: 37,620; 109: 4,711

non-Paley (610,48):
  77: 12; 85: 1,803; 93: 68,313; 99: 1,637;
  101: 269,671; 107: 42,039; 109: 16,525.
```

These frequencies distinguish the sampled landscapes sharply: both
non-Paley orbits exhibit terminal levels 77, 99, and 107 absent from the
Paley sample, while their frequencies at 109 differ. This is reproducible
heuristic evidence, not a certified invariant; absence from a sample does
not prove absence from the Paley landscape.

## 3. Local-field deficit comparison

For an oriented witness, put `u_i=sign(H) w_i(Sw)_i`. At `k=3`,

```math
114-|H|={1\over24}\sum_i(u_i-5)(u_i-7).               \tag{LC4}
```

Every orbit has an explicit exact energy-109 witness with local-field
histogram

```text
u=1: 2 coordinates, u=5: 27 coordinates, u=9: 9 coordinates.
```

Both non-Paley orbits additionally yielded exact energy-109 witnesses with

```text
u=1: 1, u=5: 30, u=9: 6, u=13: 1.
```

Both patterns have penalty 120 and deficit five. The second pattern was not
seen in the Paley search, but no exact Paley exclusion is known. The central
conclusion is negative: graph inequivalence and very different local-search
landscapes do **not** currently produce a cap gap or a different leading
local-field deficit. All three classes reach the same certified value and
share an extremal field signature.

## 4. Precise scaling hypothesis isolated by the data

The exact small cases and the order-38 certificates line up as

```math
\operatorname{cap}(S_k)=2ks-(2k-1),\qquad s=2k^2+1,   \tag{LC5}
```

for `k=1,2`, and as a certified lower bound with no better heuristic witness
for all three natural orbits at `k=3`. An exact uniform lower lemma sufficient
for the negative/saturation direction is:

> For every member of a proposed infinite self-indexed ASDS family, there is
> a Boolean spin whose oriented local-field penalty is at most
> `8k(2k-1)`.

By (LC4)'s general form, this gives

```math
\operatorname{cap}(S_k)\ge2ks-(2k-1)
=\left(\frac12-o(1)\right)(2s)^{3/2}.                \tag{LC6}
```

Together with an infinite **non-Paley** self-indexed ASDS construction,
(LC6) would be a scalable obstruction to obtaining a uniform sub-`1/2`
Boolean cap merely from this two-fiber structure. This is a concrete uniform
lemma, but the present evidence covers only `k<=3`; it is a research
hypothesis, not primary scaling progress.

Conversely, a positive landing route would need a non-Paley family whose
minimum local-field penalty grows by a positive fraction of `k^3`, rather
than the observed `O(k^2)` penalty corresponding to an `O(k)` energy deficit.
No such signal appears at order 38.

## 5. k=4 search status

For `k=4`, `s=33`, the CP-SAT formulation uses 16 independent symmetric
entries of `a`, 33 entries of `c`, their exact row sums, and all 16 independent
complementary-autocorrelation equations. Two safe symmetries are removed:

1. translate `c` so `c_0=+1`;
2. reverse `c` so `c_1>=c_{-1}`.

The normalized model reproduces a `k=3` certificate, but the five-minute
`k=4` run returned `UNKNOWN`. This is neither existence nor nonexistence and
is not counted as progress. The saved status is
`computations/results/phase2h_k4_self_indexed_asds_search.json`.

## Reproduction

```bash
PYTHONPATH=computations .venv/bin/python \
  computations/phase2h_order38_cap_local_fields.py \
  --output computations/results/phase2h_order38_cap_local_fields.json

.venv/bin/python computations/two_fiber_cyclic_conference.py \
  --search 4 --time-limit 300 \
  --output computations/results/phase2h_k4_self_indexed_asds_search.json
```

The first result is classified as heuristic search with exact witnesses. The
second is explicitly classified as a solver outcome only.

## 6. Bounded exact local-field-pattern decision

A final bounded exact audit asked whether the `(615,65)` non-Paley orbit can
have energy 113. The conference moment identities and penalty budget compress
the oriented local fields to exactly six possible histograms:

```text
5^22 7^15 11^1
5^23 7^12 9^3
3^1 5^20 7^15 9^2
3^2 5^17 7^18 9^1
3^3 5^14 7^21
1^1 5^17 7^20.
```

This enumeration uses all four exact constraints: 38 coordinates, first
moment 226, second moment 1406, and total penalty 24. Moreover, `-S` is
switching/permutation equivalent to `S`: switch the second fiber, swap the
two fibers, and reverse both cyclic indices. Thus positive energy 113 is the
only sign that must be considered.

The exact model views a spin as a Seidel switching and asks whether the
switched row-sum histogram equals one of the six lists. Global negation is
fixed and lexicographic maximality removes all 19 simultaneous cyclic shifts.
This is a genuine arithmetic and symmetry reduction from unrestricted
threshold search.

It did not, however, become a new proof mechanism. The remaining model still
contains the same 38 spin variables and 703 edge XORs, with row-field
indicators added. All six bounded trials returned `UNKNOWN`. Energy 111 has
51 exact moment/penalty histograms and is less compressed. Following the
stopping rule, no longer threshold runs were made and no cap conclusion is
drawn.

The reproducible model is
`computations/phase2i_order38_local_field_patterns.py`; the honest partial
status is
`computations/results/phase2i_nonpaley615_energy113_pattern_search.json`.

## 7. Exact `k=4` postscript

The `UNKNOWN` CP status in Section 5 records the chronology, but is no longer
the mathematical frontier. A separate exhaustive integer audit proves that
the strong cyclic ansatz has no `k=4` solution. It rejects 7,998 of 8,008
internal supports with exact integer PSD separators, reduces the ten
survivors to one multiplier orbit, and exhausts the remaining bridge profile
in 43,268,109 backtrack nodes. See
`computations/audit_k4_strong_cyclic_nonexistence.cpp` and
`computations/results/k4_strong_cyclic_nonexistence.json`. This does not
exclude general order-66 conference matrices or noncyclic two-fiber systems.
