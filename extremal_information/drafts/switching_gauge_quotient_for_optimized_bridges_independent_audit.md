# Independent audit: switching gauge quotients for optimized bridges

**Verdict: PASS.**

All covariance, optimization, fixed-label, and orbit-count statements are
correct.  In particular, the two bridge counts in the source answer
different quotient problems:

```math
\#\{\text{bare bridges modulo independent row/column gauges}\}
=2^{(m-1)(n-1)},
```

whereas

```math
\#\{\text{parent-switching classes over two fixed complete child
switching classes}\}
=2^{mn-1}.
```

The director conclusion is also correct: Theorem 36.28 measures a relative
switching gauge against one labelled fixed bridge, while independently
minimizing over the full exact-sign bridge fibre quotients both child
switches completely.  The remaining `2^(mn-1)` is only an anchored
coefficient-fibre count, not a response-packing or hardness theorem.

The source audited and frozen for this verdict is

```text
extremal_information/drafts/switching_gauge_quotient_for_optimized_bridges.md
sha256 6b24e3448bc96986c0b3a9673890914efe311dbeea15e00cf938ceebb364282c
```

No repair is required.  The phrase “minimal covariance hypothesis” should
be read as the clean natural structural guarantee used by GQ.2, not as a
logical iff: accidental equality of optimized values can hold without
covariance, as the source itself explicitly notes.

## 1. Block covariance and energy normalization: PASS

Let `D=diag(S,T)`.  Direct block multiplication gives

```math
D\begin{pmatrix}A&B\\B^T&C\end{pmatrix}D
=\begin{pmatrix}SAS&SBT\\TB^TS&TCT\end{pmatrix}
=P(A^s,C^t;SBT).
```

Since the block Hamiltonian is

```math
H_{P(A,C;B)}(x,y)=H_A(x)+x^TBy+H_C(y),
```

with no missing factor two, conjugation gives exactly

```math
H_{P(A^s,C^t;SBT)}(x,y)
=H_{P(A,C;B)}(s\odot x,t\odot y).
```

The Boolean change of variables is bijective, so the complete energy
multiset, both one-sided extrema, and the absolute cap are preserved.  GQ.1
is exact, not merely an inequality.

## 2. Exact bridge-family hypothesis and all optimized consequences: PASS

The hypothesis needed by the stated proof is the pushforward equality

```math
\mathcal B(A^s,C^t)=S\mathcal B(A,C)T.
```

Because diagonal sign matrices are involutions, `B -> SBT` is then a
bijection with inverse `B' -> SB'T`.  GQ.1 preserves each parent cap under
this bijection, proving equality of the attainable cap multisets.  It
follows immediately that their infimum, minimum, supremum, maximum, and
quantiles agree.

For a probability statement, the exact additional hypothesis is likewise
the one in the source: the law for `(A^s,C^t)` must be the pushforward of
the law for `(A,C)` under `B -> SBT`.  Then the full cap distribution and
every integrable distributional statistic agree.  If an auxiliary cost is
included, it must obey `cost(A^s,C^t;SBT)=cost(A,C;B)`; the same argument
then applies to the cap-plus-cost multiset.

The unrestricted exact-sign bridge family is covariant.  Operator norm,
rank, and a prescribed zero/support pattern are also invariant under left
and right diagonal signs.  A sign-sensitive constraint which is not closed
under this transport would not meet GQ.4 and receives no invariance
conclusion.

In particular,

```math
\min_{B\in\{+-1\}^{m\times n}}Q(P(A^s,C^t;B))
=\min_{B\in\{+-1\}^{m\times n}}Q(P(A,C;B)).
```

This remains true with a fixed public second child `C` when the bridge is
reoptimized separately for each first-child switch, using `B -> SB`.
Covariance is sufficient rather than logically necessary for an isolated
numerical equality, exactly as Section 2 states.

## 3. Abstract labelled versus unlabelled language: PASS

From

```math
V(gX,\tau_gk)=V(X,k)
```

one obtains

```math
r_{gX}(k')=r_X(\tau_g^{-1}k').
```

Thus a transported response table is generally a permutation of labels,
not the same pointwise table.  Any relabelling-invariant readout—minimum,
unordered multiset, or integral under the transported measure—factors
through `X/G`.  A laboratory-fixed context language does so precisely when
the table happens to be pointwise invariant on its common labelled fibre.
If one fixed label exposes a gap `delta`, a common decoded quotient state
cannot approximate both systems below error `delta/2`.  These are exactly
the contextual/Myhill--Nerode quantifiers used elsewhere in the repository.

For block systems, switching the first child sends the natural label
`(C,B)` to `(C,SB)`.  Consequently the table in Section 3 correctly
distinguishes optimized/transported languages from a pinned bridge or
noncovariant bank.  The “common switch” entry describes the structural
guarantee; special accidental cap equalities are not excluded.

## 4. Four-vertex fixed-bridge falsifier: PASS

With `A=J_3-I_3` and all three bridge entries negative, switching the fourth
vertex turns the parent into the all-positive `K_4` signing, whose cap is
six.

After switching only the internal `K_3` block by `s=(-1,-1,1)` while
holding the bridge fixed, set the fourth spin to `+1`.  The energy is

```math
x_1x_2-x_1x_3-x_2x_3-x_1-x_2-x_3.
```

Over the eight triples in lexicographic order it takes the values

```text
2, 4, 0, -2, 0, -2, 2, -4,
```

so the absolute cap is four.  Global spin inversion justifies fixing the
fourth spin.  Transporting the bridge gives
`SB_0=(1,1,-1)^T` and restores cap six by GQ.1.  The example therefore
falsifies fixed-label invariance and nothing stronger.

## 5. Bare bridge quotient: PASS

The row/column switching group has `2^(m+n)` elements.  If `SBT=B` for a
complete exact-sign bridge, then

```math
s_it_j=1\qquad\text{for all }i,j.
```

All `s_i` and all `t_j` must be one common sign, so the kernel has exactly
two elements.  Every orbit consequently has size `2^(m+n-1)`, and dividing
the `2^(mn)` bridges gives

```math
2^{mn-m-n+1}=2^{(m-1)(n-1)}
```

orbits.  Gauging the first row and first column positive is unique, and the
remaining entries are exactly the rectangle products

```math
B_{ij}B_{i1}B_{1j}B_{11},\qquad i,j\ge2.
```

They form a complete set of `(m-1)(n-1)` binary invariants.

## 6. Fibre over anchored complete child classes: PASS

For a complete internal signing of order `m>=2`, diagonal switching is free
modulo the global shore sign: `SAS=A` forces `s_is_j=1` for all `i ne j`,
hence `s=+-1`.  Its coordinate-labelled switching orbit therefore has
`2^(m-1)` members, and similarly the order-`n` orbit has `2^(n-1)`.

There are two equivalent independent checks of GQ.5.

First, the triple space has size

```math
2^{m-1}\,2^{n-1}\,2^{mn}=2^{mn+m+n-2}.
```

Complete-parent switching has effective orbit size `2^(m+n-1)` and acts
freely modulo its one global sign.  The quotient therefore has size

```math
2^{mn+m+n-2-(m+n-1)}=2^{mn-1}.
```

Second, gauge `A'` and `C'` back to fixed representatives `A,C`.  Each
returning shore switch is unique only up to its global sign.  Reversing the
chosen sign-vector representative on exactly one shore changes `B` to
`-B`; reversing both leaves it fixed.  Completeness gives no further
stabilizer.  Thus the residual relation is precisely `B~-B`, again leaving
`2^(mn-1)` classes.

This is not the bare-bridge quotient.  Fixing the internal representatives
uses the `m+n-2` nontrivial row/column gauge bits and makes their alignment
relative to the child coefficients observable.  The decomposition

```math
mn-1=(m-1)(n-1)+(m-1)+(n-1)
```

correctly separates rectangle holonomy from the two alignment labels.
Neither calculation says that distinct anchored coefficients have distinct
optimized caps.

## 7. Consequence for Theorem 36.28: PASS

Theorem 36.28 freezes one bridge `B` while its old child and query fill
range through named switchings.  Under an old switch `S_s` and fill switch
`S_q`, covariance would transport that bridge to `S_sBS_q`; the theorem
deliberately retains `B` instead.  Its scalar gaps therefore measure a real
relative gauge against a fixed public frame.

If the bridge is instead optimized independently over every exact-sign
matrix, GQ.2 gives the stronger collapse

```math
\min_BQ(P(A^{s_i},D^{q_j};B))
=\min_BQ(P(A,D;B))
```

for every pair `(i,j)`.  Thus the whole switching response table used in
Theorem 36.28 becomes constant under that optimized future language.  Its
`Theta(n)` labelled response lower bound cannot be transferred to the
original bridge minimum.

This does not solve optimized composition.  Once representatives are
anchored, GQ.5 still leaves `2^(mn-1)` bridge coefficients modulo the sole
antipodal redundancy, and the present theorem proves no smaller sufficient
state for minimizing over them.  Accordingly the source's director
conclusion is exact:

- switch labels themselves are pure gauge for independently optimized
  covariant bridges;
- the dense anchored bridge fibre is the remaining search object;
- Theorem 36.28 remains relevant only to fixed/noncovariant/shared-bridge
  architectures or pinned laboratory-frame queries.

It proves no optimized-bridge information lower bound, recurrence, or
convergence result.

## 8. Archive comparison: PASS

The quadratic dense-bridge audit already observes that switching a child
transports a generic fixed bridge, so GQ.1's fixed-label warning is a
formalized archive collision rather than a new phenomenon.  The polynomial
fractional-bridge report proves covariance for a special universal-double
family, and Theorem 21.41 gives signed-graph gauge/holonomy for a structured
common-Hadamard language.  None states the present general bridge-family
factorization together with both fibre counts.

The exponent `2^(mn-1)` also appears in the microcanonical disorder-counting
composition theorem, where it comes from bridge/orientation parameter
multiplicity after an antipodal representative choice.  That numerical
collision is compatible but is not GQ.5's quotient classification.

The source accurately labels its increment as an elementary scope theorem,
not a new lower bound.  Its principal value is the correction to how the
new fixed-reference response theorem may—and may not—be used in director
reasoning about the original optimization.
