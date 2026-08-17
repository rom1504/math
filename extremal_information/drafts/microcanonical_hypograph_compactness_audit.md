# Adversarial audit: microcanonical hypograph compactness

**Verdict: REPAIR.**  The compact hypograph sup-convolution theorem is
correct, including its treatment of empty fibres, and the finite
count-convolution estimate is correct after a uniformity/zero-count repair.
The draft is not ready to promote as written, however.  Exact-fibre
hypographs do **not** by themselves control total exponential mass or the
uniform probability of a recovered branch when the descriptor image may
have exponentially many points.  The BEG paragraph also mixes the raw-count
normalization of (MH.1) with the uniform-prior probability normalization of
the cited paper.  Finally, this is best classified as an adjacent,
speed-sensitive instance of contextual response theory, not as a genuinely
orthogonal algebra.

## 1. Theorem MH.1

### 1.1 Sup-convolution: pass

Interpret every profile as taking values in `[-infinity,M]`, with no
`+infinity` values.  Let

```math
h_n(z)=\max_{m(x,y)=z}(f_n(x)+g_n(y)),
\qquad
h(z)=\max_{m(x,y)=z}(f(x)+g(y)).
```

The stated proof gives both defining halves of hypograph convergence.

- **Upper bound.**  If `z_n -> z` and the relevant limsup is finite, choose
  exact maximizing decompositions `(x_n,y_n)`.  Compactness gives a
  convergent subsequence.  Continuity gives `m(x,y)=z`, while
  hypograph upper semicontinuity gives

  ```math
  \limsup_n h_n(z_n)
  \le f(x)+g(y)
  \le h(z).
  ```

  Hypograph convergence on a compact domain also supplies the eventual
  uniform upper bounds needed when splitting the limsup: otherwise points
  at which a component tends to `+infinity` would have a convergent
  subsequence contradicting its hypograph upper bound.

- **Recovery.**  If `h(z)>-infinity`, upper semicontinuity on the compact
  fibre makes the maximum attainable.  Separate full-sequence recovery
  points `x_n -> x` and `y_n -> y` give

  ```math
  z_n=m(x_n,y_n)\to z,
  \qquad
  \liminf_n h_n(z_n)\ge f(x)+g(y)=h(z).
  ```

  If `h(z)=-infinity`, the recovery inequality is vacuous, as it should be.

No convexity is being used.  The result is the compact, sign-dual analogue
of the epi-addition stability result cited from Rockafellar--Wets.  The
reference to Proposition 7.56(a) is directionally accurate, although that
proposition is stated for addition on Euclidean spaces under total
epi-convergence and a horizon transversality condition; the elementary
compact proof here is what actually establishes the more general continuous
map `m` version.

### 1.2 Empty fibres: pass

Because `m(K_1 times K_2)` is compact and hence closed, if
`m^{-1}(z)` is empty then it remains empty in a neighbourhood of `z`.
Thus both the upper and recovery conditions hold with value `-infinity`.
No separate selection argument is required.

### 1.3 Tilted maxima and optimizers: pass with one degenerate-case qualifier

Adding a continuous `V` preserves hypograph convergence.  Compactness then
gives convergence of maximum values and the outer-limit assertion for exact
maximizers.  If the limiting maximum is finite, recovery at any limiting
maximizer is automatically asymptotically maximizing, because the recovery
lower bound and convergence of global maxima squeeze its value.

The phrase “every limit maximizer has an asymptotically maximizing recovery
sequence” should either assume

```math
\max_z\{h(z)+V(z)\}>-\infty
```

or define “asymptotically maximizing” in the topology of the extended real
line.  When the entire limiting objective equals `-infinity`, recovery is
vacuous and an additive optimality-gap formulation is not meaningful.

### 1.4 Compact-hyperspace sentence: repair

“Closed hypographs in the compact cylinder form a compact hyperspace” is not
literally correct with ordinate space `R`: a full hypograph is unbounded
below and `K times R` is not compact.  One may instead:

1. normalize to a common upper bound `M` and use the compactified ordinate
   `[-infinity,M]`; or
2. truncate hypographs to `K times [-R,M]` for every `R` and use the local
   Attouch--Wets topology.

With the first convention, the nonempty compact subsets of the cylinder do
form a compact hyperspace, and the downward-closed property is closed under
limits.  Individual bounded-above profiles without a common normalization
do not form one compact family.

## 2. Corollary MH.2

### 2.1 Exact count-convolution comparison: pass after precise hypotheses

For positive `C_n(z)`, let

```math
D_n(z)=
\#\{(x,y):m(x,y)=z,\ A_n(x)B_n(y)>0\}.
```

Then “largest summand at most sum at most number of summands times largest
summand” proves (MH.8).  If `C_n(z)=0`, both the log-count profile and the
max-profile are `-infinity`; their difference should not be written as an
ordinary real number.  A clean statement is

```math
0\le h_n^{\rm count}(z)-h_n^{\rm max}(z)
\le {\log\max(1,D_n(z))\over a_n}
```

on the common effective domain, with both profiles `-infinity` off it.
To transfer a hypograph limit, require the uniform condition

```math
\sup_z\log\max(1,D_n(z))=o(a_n)
```

and explicitly assume that the two input log-count profiles hypographically
converge.  Pointwise `o(a_n)` depending on `z` is not sufficient for arbitrary
moving sequences `z_n`.

### 2.2 Branch probability claim: false without descriptor-complexity control

Hypograph convergence of **exact-fibre** log multiplicities controls the
largest fibre near a descriptor.  It does not control the sum over all
descriptor values.  Here is a decisive counterexample.

For each `n`, put `ceil(exp(a_n))` distinct descriptor points in a compact
interval, give every point multiplicity one, and let those points become
dense.  The exact-fibre log profile hypographically converges to the constant
zero profile.  Nevertheless the total number of states has exponent one,
and the uniform probability of each fibre is `exp(-a_n+o(a_n))`, not
`exp(o(a_n))` as “below the maximum by `Delta=0`” would predict.

Therefore the sentence following (MH.8) is valid only under an additional
Laplace-principle hypothesis, for example

```math
\log |Q_n(X_n)|=o(a_n),
```

or a uniform subexponential covering-number/coarse-bin condition ensuring

```math
{1\over a_n}\log\sum_q A_n(q)
=\max_q {1\over a_n}\log A_n(q)+o(1).
```

Under that condition, if a recovered fibre has entropy exactly `Delta`
below the global maximum, its **fibre event** (not each state in it) has
uniform probability `exp(-a_n Delta+o(a_n))`.  Without the condition, a
proper large-deviation state should use local-ball masses before sending the
ball radius to zero, rather than exact fibres.

This is a scope issue for the proposed rare-event state, not merely a wording
issue.  It must be fixed anywhere the hypograph is claimed to answer
log-sum/pressure or probability queries.

### 2.3 Abstract realizability: pass, but it is deliberately weak

For a bounded nonnegative usc `s` on compact `K`, compactify/truncate the
ordinate and take finite Hausdorff nets of its hypograph.  Retain the top
sampled height at each projected grid point and assign

```math
N_n(q)=\left\lceil e^{a_n r_n(q)}\right\rceil.
```

The height error is at most `log 2/a_n`.  By diagonalizing the mesh speed one
can ensure both mesh size tending to zero and logarithmic grid cardinality
`o(a_n)`, because every fixed-scale cover of a compact metric space is
finite.  This constructs the claimed hypographic recovery sequence.

The statement realizes `s` only as an **unstructured abstract landscape with
chosen descriptor multiplicities**.  It gives no recovery theorem inside a
specified Ising, code, graph, or quadratic-form class.  This limitation
should be explicit; otherwise “realizability” sounds much stronger than what
is proved.  Profiles taking `-infinity` or profiles unbounded below require a
separate convention and are not covered by the stated nonnegative version.

## 3. BEG benchmark and bounded-temperature falsifier

### 3.1 BEG mapping: correct modulo normalization

The cited BEG paper uses

```math
U(x^n)=\sum_i x_i^2-{K\over n}\left(\sum_i x_i\right)^2,
\qquad x_i\in\{-1,0,1\},
```

and its empirical occupation vector gives exactly

```math
u(L)=L_++L_- -K(L_+-L_-)^2.
```

It takes the a priori law to be uniform with mass `3^{-n}` per spin word, so
its large-deviation log-**probability** profile is

```math
s_{\rm prob}(L)=-\sum_jL_j\log L_j-\log3.
```

By contrast, definition (MH.1) is a log-**count** profile.  In that
normalization the occupation entropy is

```math
s_{\rm count}(L)=-\sum_jL_j\log L_j.
```

The two differ only by the harmless global constant `log 3`, but the draft
must say which convention it is using.  The constrained supremum over
`u(L)=u` and the claim about a nonconcave interval at ensemble-nonequivalent
parameters agree with the primary source.  The canonical family retains only
the Legendre--Fenchel/concave-envelope information, so it does not recover
that nonconcave part.

### 3.2 Bounded-beta example: pass

Uniformly for `|beta|<=B`,

```math
0\le {1\over n}\log\left(1+{e^{\beta n\delta}\over\lceil e^n\rceil}
                      \right)
\le {1\over n}\log(1+e^{-n(1-B\delta)+o(n)})\to0
```

when `delta<1/B`.  Yet the maximum energy density is `delta` in the planted
landscape and zero in the bulk-only landscape.  This is a valid uniform
falsifier for **bounded** inverse-temperature queries.  It does not show that
the complete family of all temperatures, including `beta -> infinity`,
misses the maximum, and the draft does not claim otherwise.

The final warning is also correct: speed-`n` hypographs identify one maximal
state with `exp(sqrt(n))` maximal states and cannot retain extremal-process
decorations or genealogy.

## 4. Is this genuinely orthogonal?

Not at the level of its algebra.  If the declared queries are all continuous
tilts `V`, then

```math
V\longmapsto\sup_q\{s(q)+V(q)\}
```

is exactly a contextual response functional, and (MH.4) is its ordinary
sup-convolution composition law.  The proof is standard hypographical/Gamma
compactness.  Thus MH.1 should not be advertised as a new algebra independent
of response roofs.

There is nevertheless a real scoped distinction.  A finite-dimensional
upper response roof under **linear** tilts convexifies its lifted data and
cannot retain a nonconcave microcanonical branch.  The usc hypograph, queried
by arbitrary continuous tilts and scaled log multiplicities, refuses that
convexification and records zero-density branches at their exponential
scale.  It is therefore a useful rare-event compactness branch inside the
broader contextual-response program, and a good diagnostic for which query
class is being declared.  It is not yet an orthogonal route to structured
extremal realization.

## 5. Required repairs before promotion

1. Compactify the ordinate or state a truncation/local-hypograph topology and
   impose a common upper normalization.
2. Qualify optimizer recovery when the limiting optimum is identically
   `-infinity`.
3. Define `D_n(z)`, handle zero counts in the extended-real convention, and
   require the decomposition-count bound uniformly in `z`.
4. Add a subexponential descriptor-image/coarse-bin hypothesis before any
   total-mass or uniform branch-probability conclusion; say “fibre event,”
   not “states.”
5. Scope realizability to arbitrary abstract landscapes, not a structured
   model class.
6. Reconcile the BEG raw-count and uniform-prior probability normalizations.
7. Classify the result as a scoped speed-sensitive response/hypograph branch,
   not as a fully orthogonal theory.

After these repairs, the mathematical core merits promotion as a rigorous
benchmark theorem and bounded-temperature falsifier.
