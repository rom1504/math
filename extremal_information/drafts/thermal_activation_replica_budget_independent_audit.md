# Independent audit of `thermal_activation_replica_budget.md`

**Verdict:** the core two-point theorem is correct, including its exponential
rate, adaptive-sampling coupling, sharp strict-exponent converse, and binary
mutual-information consequence.  The nonflat-background identities and hard
conjunctive product are also correct after making positivity explicit.
Canonical promotion should nevertheless wait for two substantive repairs:
the pressure-resolution/minimax sentence is not presently a well-defined
observation theorem, and the prose at the critical replica exponent is
stronger than the proved hypotheses.  Three scope/normalization repairs are
also needed.

This audit checked the draft against repository Theorem 7.1, Theorem 27.1,
Proposition 27.2, Theorem 27.3, and
`artifacts/two_temperature_bridge_audit.md`.

## 1. Exact marked-phase and background identities

### Flat two-level carrier: PASS

For `H^(1)=E 1_P`, direct splitting over `P` and `P^c` gives

```math
Z_1=1-q+qe^{\beta E},
\qquad
p={qe^{\beta E}\over1-q+qe^{\beta E}}.
```

Thus `(q,E)` exactly answers the declared scalar pressure and marked-mass
queries.  The hard conjunction in (TA.6) has marked fraction `qq'` and reward
`E+E'`, so (TA.7) is exact.  The draft correctly says that this is **not**
ordinary additive Hamiltonian composition.

One normalization phrase should be repaired.  At finite size it is the
extensive rarity cost `J=-log q`, not the intensive rate `I=J/n`, that adds:

```math
J_{\wedge}=J+J',\qquad E_{\wedge}=E+E'.
```

For systems of sizes `n,m`, the combined intensive rates are weighted
averages, not sums.  Replace “rarity exponents and rewards add” by this exact
statement or define “rarity exponent” explicitly as `-log q`.

### Lemma TA.0: PASS after a positivity repair

For a nonflat background `K`, (TA.7b)--(TA.7c) follow by changing measure.
Conditional on `P` and `P^c`, the old and tilted Gibbs laws agree.  Therefore

```math
\|\mu_{K^+}-\mu_K\|_{TV}=|p_K-q_K|.
```

The displayed equality `p_K-q_K`, monotonic language (“boost”), and the
one-sided bound in (TA.7g) require `E>=0`.  The draft never states this in the
nonflat paragraph.  **Required repair:** state `E,F>=0` in Lemma TA.0, or put
absolute values in (TA.7d), (TA.7g), and the associated prose.  Positivity is
already present asymptotically in TA.1 through `delta>0`, but it does not
logically repair the standalone lemma.

The factorization

```math
q_{K\oplus L,P\times Q}(\beta)=q_{K,P}(\beta)q_{L,Q}(\beta)
```

is exact, so (TA.7f) passes.  Notice that the state here contains an entire
function `beta -> q_K(beta)`, not two real numbers.  It is a valid presented
state, but no finite-description or computational-compression claim should
be inferred for arbitrary `K`.

Sequential maximal coupling proves (TA.7g) for a deterministic sample budget
or a stopping rule bounded by `K_samp`.  The phase-counting converse is also
correct under the two stated conditions.  If an unbounded random stopping
time is intended, an expectation or tail condition must replace the current
deterministic bound.

## 2. Thermal exponent and pressure asymptotics

### (TA.10)--(TA.11): PASS

Write

```math
u_n(\beta)=q_n(e^{\beta E_n}-1).
```

The hypotheses imply, uniformly for `0<=beta<=B`,

```math
u_n(\beta)\le e^{-n(I-B\delta)+o(n)}.
```

At `B`, positivity of `B delta` makes the first term dominate the subtracted
`q_n`, and hence

```math
u_n(B)=e^{-na_B+o(n)}=o(1).
```

Consequently `n^(-1)log(1+u_n(B))=n^(-1)e^{-na_B+o(n)}`.  The factor `1/n`
and all exponential rates in TA.1 are correct.

The phrase “exponential precision needed” is valid only as an
**exponential-rate** statement.  Exact finite pressure values distinguish
the candidates immediately; a precision lower bound needs a specified
noisy or quantized observation model, addressed in Section 5 below.

## 3. Adaptive Gibbs experiments

### One-sample total variation: PASS

For the flat model, both laws are uniform inside the two cells.  Since
`E_n>0` asymptotically,

```math
d_n(\beta)=p_n(\beta)-q_n
```

is exact and increases with `beta`.  In the subcritical regime,

```math
p_n(B)=e^{-na_B+o(n)},
\qquad
d_n(B)=e^{-na_B+o(n)},
\qquad
{q_n\over p_n(B)}=e^{-nB\delta+o(n)}.
```

All three exponents are correct.

### Adaptive coupling: PASS

Couple the observer's internal random seed.  While the histories agree, its
chosen temperatures agree; maximal coupling then fails at the next round
with probability at most `d_n(B)`.  Thus

```math
d_{TV}(\mathcal L_1,\mathcal L_0)\le K_nd_n(B).
```

No independence across the adaptively chosen temperatures is being silently
assumed beyond the stated fresh Gibbs sample at each round.  Condition
(TA.12) gives an exponentially vanishing right side.

The “energies are post-processings” sentence needs a small semantic repair.
Energy features computed under each of the two **known candidate
Hamiltonians** are post-processings of a sampled configuration.  An oracle
that separately reveals the energy under the unknown candidate is an extra
observation channel and is not literally a post-processing of the
configuration alone.  Clarify which meaning is intended.  (For this special
two-level model that extra label has the same exponential detection
threshold, but (TA.18) as written audits configuration-only transcripts.)

### Mutual information: PASS

With a uniform binary prior, the mutual information is the
Jensen--Shannon divergence of the two transcript laws.  Vanishing total
variation forces this divergence to vanish (the likelihood ratios against
the mixture are bounded by two).  Thus the information assertion is valid;
it is not using a false reverse Pinsker inequality.

### Sharp converse: PASS in the stated strict-exponent regime

Under (TA.13), `K_np_n(B)->infinity`, while `q_n/p_n(B)->0`.  The threshold
`K_np_n(B)/2`, Chernoff under `s=1`, and Markov under `s=0` give precisely
(TA.20).  This remains valid even when `K_nq_n` itself diverges.  Hence the
strict lower/upper exponential thresholds in (TA.12)--(TA.13) are correct.

There is a sharper exact formulation already implicit in the proof:

```math
K_nd_n(B)\to0 \quad\Longrightarrow\quad\text{blindness},
```

whereas

```math
K_np_n(B)\to\infty,
\quad q_n/p_n(B)\to0
\quad\Longrightarrow\quad\text{consistent distinction}.
```

Since `d_n(B)~p_n(B)`, these identify the threshold up to its critical
window.  At `log K_n/n -> a_B`, the `o(n)` terms decide the answer.

**Required repair:** replace the imprecise claims “fewer than
`exp(n(a_B-o(1)))`” in line 229 and “more than
`exp(n(a_B+o(1)))`” in the falsifier with the exact product conditions above,
or retain only the strict limsup/liminf conditions (TA.12)--(TA.13).  As
written, those phrases can include, for example, `K_n=e^{na_B-1}`, for which
`K_np_n(B)` need not tend to zero, so the claimed `delta/2-o(1)` risk does
not follow.

## 4. Zero-temperature gap

(TA.14) is immediate and correctly normalized:

```math
{\max H_n^{(1)}-\max H_n^{(0)}\over n}={E_n\over n}\to\delta.
```

For transcript-only observations satisfying `K_nd_n(B)->0`, the standard
common-part form of Le Cam's two-point argument indeed gives minimax absolute
risk

```math
{\delta+o(1)\over2}(1-d_{TV}(\mathcal L_0,\mathcal L_1))
=\delta/2-o(1).
```

There is no missing factor of two here.

## 5. Pressure-resolution/minimax clause: NOT YET A THEOREM

Lines 227--231 and 306--308 refer to “uniform resolution” and adversarial
pressure perturbations without defining the observation model.  Literal
exact access to the pressure curve distinguishes the two finite landscapes
from any nonzero difference.  A deterministic rounding grid can also put
arbitrarily close values in different bins.  Thus a scale comparison alone
does not imply the asserted minimax result.

**Required repair:** define a robust `L^infinity` pressure oracle.  For
example, set

```math
F_{n,s}(\beta)={1\over n}\log Z_{n,s}(\beta),
\qquad
\Delta_n^P=\|F_{n,1}-F_{n,0}\|_\infty=F_{n,1}(B).
```

Suppose the recorded curve may be chosen adversarially within `eta_n` of the
true curve.  If

```math
2\eta_n\ge\Delta_n^P,
\qquad K_nd_n(B)\to0,
```

the adversary can return the midpoint curve under both candidates, and Le
Cam gives the advertised `delta/2-o(1)` lower bound.  This exact statement
has the exponent

```math
\Delta_n^P={1\over n}e^{-na_B+o(n)}.
```

Alternatively, delete the pressure-plus-sample minimax sentence and retain
the rigorously proved pressure separation and transcript no-go as separate
claims.

## 6. Code/CSP benchmark and conjunctive composition

The code calculation passes.  For a binary linear `[n,k_n]` code,

```math
q_n=2^{k_n-n},
\qquad
I=(1-R)\log2,
```

and the replica threshold has the stated exponential rate.  Direct product
adds codimensions and rewards exactly.  The CSP statement likewise passes
for a hard all-constraints-satisfied bonus and disjoint conjunction.

Two scope qualifications should remain explicit in any canonical version:

1. This is a **membership-bonus** code landscape, not a covering-radius,
   coset-distance, or decoding-energy theorem.
2. `(q,E)` answers the declared permutation-invariant pressure/phase-mass
   queries and hard conjunction.  It does not reconstruct labelled Gibbs
   sample laws, arbitrary geometric futures, or ordinary additive
   composition.  The draft's (TA.25) is a correct exact falsifier for such an
   extension.

The phrase in (TA.23) that the experiment “needs” that many replicas should
be read as an exponential-rate threshold.  For a literal finite-size sample
complexity, use `1/p_n(B)` and retain its subexponential factor.

## 7. Archive collision and novelty classification

### Theorem 7.1

There is no collision with the posterior-width theorem.  Theorem 7.1 prices
mutual information needed to encode one of many geometrically separated
response vectors.  TA.1 uses a binary family and shows that a particular
sampling query acquires vanishing information because it assigns
exponentially small mass to the distinguishing event.  The two theorems can
be combined, but neither is a renaming of the other.

### Theorems 27.1--27.3

Proposition 27.2 already contains the same rare-state pressure obstruction as
the flat special case of TA.1.  Theorem 27.1 explains why retaining the
microcanonical branch recovers it, and Theorem 27.3 already supplies the
finite phase-count multiplication algebra.  Thus the two-level carrier,
entropy--reward balance, and hard product factorization are elementary
specializations of existing count/pressure structure.  The genuine
repository increment is narrower: the exact adaptive transcript-TV bound,
the sharp replica exponent, and the resulting operational minimax statement
once its noise model is made formal.

### Entropy-tilted bridge audit

Equations (2.4), (3.5), and (7.2) of
`artifacts/two_temperature_bridge_audit.md` already express the same
entropy-versus-reward principle for a rare bridge basin: a set of probability
`p` costs `log(1/p)` under exponential tilting.  That artifact studies
change-of-measure/optimization over disorder bridges, whereas TA.1 studies
hypothesis testing from Gibbs samples over configurations.  Its exact sample
threshold is not proved in the bridge audit, but the basic activation balance
in (TA.24) should not be advertised as a new principle.  Add this comparison
to the scope section.

**Classification:** after repairs, this is a rigorous, scoped **Level-3
operational theorem/no-go**, not a new adversarial-statistical-mechanics
architecture.  It gives a useful theorem outside the signing problem, but
only for a deliberately hard-conjunctive marked phase.  The draft's current
portfolio judgment is otherwise appropriately conservative.

## 8. Required repairs before promotion

1. State `E,F>=0` in Lemma TA.0 (or use absolute TV differences).
2. Formalize the pressure observation/noise model as above, or remove the
   pressure-plus-replica minimax clause.
3. Replace both `exp(n(a_B plus/minus o(1)))` sample-budget slogans by exact
   `K_nd_n(B)` / `K_np_n(B)` conditions or strict exponent inequalities; add
   that the critical sample window depends on subexponential factors.
4. Say explicitly that `-log q` adds under unequal-size conjunction, while
   intensive rarity rates are weighted averages.
5. Clarify the “energy statistic is a post-processing” sentence and add the
   collision with Proposition 27.2/Theorem 27.3 and the entropy-tilted bridge
   audit.

After these repairs, all displayed asymptotic exponents, the adaptive-TV
argument, the sharp strict-exponent converse, mutual-information statement,
hard conjunction, and code/CSP normalization pass independent audit.

## 9. Post-repair disposition

**PASS (2026-08-17).**  The revised draft implements all five repairs.  In
particular, positivity now licenses the one-sided TV formulas; extensive and
intensive rarity normalizations are separated; (TA.14a)--(TA.14b) define an
actual adversarial `L^infinity` pressure-oracle experiment; the replica
claims use strict exponents or the exact `K_nd_n(B)` / `K_np_n(B)` criteria;
and the code/CSP and archive comparisons have the required scope.  The
midpoint-oracle argument and the remaining Le Cam bound are correct.  No
further mathematical repair is required for promotion as a scoped Level-3
theorem/no-go.
