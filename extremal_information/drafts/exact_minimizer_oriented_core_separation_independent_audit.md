# Independent audit: exact-minimizer oriented core separation

**Verdict on the frozen source:** **PASS WITH SCOPE REPAIRS.**  OC.1,
OC.2, and OC.3 are correct, including the near-order estimate, orientation
forcing, cut factors, and minimal-size near-clique conclusion.  The proposed
`L_balance` consequence needs a sequential quantifier and a threshold below
the known lower constant.  The archive comparison must credit one-sided
partition superadditivity as pre-existing, and the result refines the shape
of the obstruction without making the missing lemma demonstrably weaker.

## 1. Frozen source

```text
extremal_information/drafts/exact_minimizer_oriented_core_separation.md
sha256 830bb63273683ce216257ff4977853a6b288f47f0d5d0128584a6a967e458bb3
```

The audit uses

```math
H_A(x)={1\over2}x^TAx,
\quad P(A)=\max H_A,
\quad N(A)=-\min H_A,
\quad Q(A)=\max\{P(A),N(A)\}.
```

For every hollow matrix, `E_X H_A(X)=0`, so `P,N>=0`; this justifies the
nonnegativity used in OC.10 and below.

## 2. OC.1: one-sided partition superadditivity

Let `u,v` maximize the positive energies of `A[T]` and `A[R]`.  The two
full spins `(u,v)` and `(-u,v)` have common internal contribution

```math
P(A[T])+P(A[R])
```

and opposite cross contributions.  At least one cross contribution is
nonnegative.  Hence

```math
P(A)>=P(A[T])+P(A[R]).
```

Applying this statement to `-A` gives the `N` inequality with no change of
normalization.  OC.2 is therefore exact.

This lemma is not new to the archive.  It appears explicitly in
`artifacts/nested_restriction_paving.md`, equation (1), in
`artifacts/concentration_compactness_boolean_profiles.md`, equation (8),
and in ledger equation (10.13).  The oriented-core synthesis below can be
new even though its first ingredient is not.

## 3. OC.2: near-order heredity and orientation forcing

Write `k=o(n)` and `m=n-k`.  Theorem 36.15 gives

```math
M_n\le M_m+M_k+\sqrt{2(\log2)mkn}.
```

Its random-sign bound gives `M_k=O(k^(3/2))`.  Since `k=o(n)`, both the
last term and `M_k` are `o(n^(3/2))`.  Conversely, taking any `m`-vertex
principal restriction of an exact order-`n` minimizer and averaging the
omitted spins gives `M_m<=M_n`.  Therefore

```math
0\le M_n-M_m=o(n^(3/2)),
```

which verifies OC.6 for every sublinear deletion sequence; no density or
regularity assumption on `T` is hidden.

For `R=[n]\setminus T`, principal monotonicity and the definition of `M_m`
give

```math
M_m\le Q(A[R])\le M_n,
```

so `Q(A[R])=M_n-o(n^(3/2))`.  Positive superadditivity and the orientation
`P(A)=M_n` give

```math
P(A[R])\le M_n-P(A[T])
          \le M_n-(t-o(1))n^(3/2).
```

The fixed `t>0` gap is much larger than the `o(n^(3/2))` uncertainty in
`Q(A[R])`.  Hence, for all sufficiently large orders, the complement's
absolute cap cannot be its positive cap:

```math
N(A[R])=Q(A[R])=M_n-o(n^(3/2)).
```

Negative superadditivity and `N(A)<=Q(A)=M_n` then imply

```math
0\le N(A[T])\le M_n-N(A[R])=o(n^(3/2)),
```

and also

```math
M_n-o(n^(3/2))\le N(A[R])\le N(A)\le M_n.
```

This proves all five lines in OC.5.  The conclusion does not identify a
full optimizer or control the cross block; it is genuinely a one-sided
principal-cap statement.

## 4. OC.3: signed-cut identities

After switching a positive core maximizer to `1`, let `b` be the switched
core and let `1^S` flip exactly the coordinates in `S`.  Only crossing
edges change sign, and each changes its energy contribution by `-2b_(ij)`.
Thus

```math
H_b(1^S)=P-2w(S).
```

Every Boolean spin is represented up to the irrelevant global sign this
way.  Maximality of `1` gives `w(S)>=0`, while minimizing over `S` gives

```math
\max_S w(S)={P+N\over2}.
```

A uniform random cut crosses each edge with probability one half, so

```math
E_Sw(S)={1\over2}\sum_(i<j)b_(ij)={P\over2}.
```

Consequently

```math
\max_Sw(S)-E_Sw(S)={N\over2}=o(n^(3/2)),
```

and all factors in OC.11--OC.13 are correct.  Singleton cuts are exactly
the switched signed row sums; their nonnegativity and total `2P` follow.

### Minimal-size near-clique clause

Let `k=(sqrt(2t)+o(1))n^(3/4)`.  Then

```math
\binom{k}{2}=(t+o(1))n^(3/2).
```

The assumed lower bound on `P(A[T])`, together with the trivial upper bound
`P(A[T])<=binom(k,2)`, yields

```math
\binom{k}{2}-P(A[T])=o(n^(3/2))=o(k^2).
```

In the maximizing gauge, if `r` core edges are negative, then
`P=binom(k,2)-2r`; hence `r=o(k^2)`.  The near-clique statement is correct.
When `k>>n^(3/4)`, OC.2 only gives `P=O(n^(3/2))`, so the average row bias
as a fraction of `k` is `O(n^(3/2)/k^2)=o(1)`; the stated limitation is
also correct.

## 5. Required repairs to `L_balance`

The displayed phrase "every `T=o(n)`" is not a finite-order quantifier.
A correct version is:

> There exist `0<t<c_-` and `c>0` such that, for every sequence of globally
> oriented exact minimizers `A_(n_j)` and every sequence
> `T_j subseteq[n_j]` with `|T_j|/n_j->0`,
> ```math
> \liminf_j {P(A_(n_j)[T_j])\over n_j^(3/2)}\ge t
> \quad\Longrightarrow\quad
> \liminf_j {N(A_(n_j)[T_j])\over n_j^(3/2)}\ge c.
> ```

The restriction `t<c_-` is needed for the advertised near-top implication:
then one may take `d_0=c_--t>0`.  Merely asserting the existence of some
`t>0` could choose a threshold above the guaranteed oriented ground scale
and would not establish the intended `L_tail` statement.

With this repair, failure of a fixed-rate upper-tail bound at threshold
`t` gives a PC.3 core with the left-hand property, while OC.2 forces its
negative cap to be `o(n^(3/2))`, contradicting the right-hand property.
The implication chain is sound.

There is an important classification point.  On exact minimizers OC.2
automatically turns every fixed-positive-level sublinear positive core into
a strongly one-sided core.  Therefore `L_balance` is, on this class, an
equivalent way (up to threshold slack) to exclude the PC core.  It refines
the geometry of a possible counterexample, but is not yet a demonstrably
weaker or easier missing lemma than `L_core`/`L_tail`.

## 6. Archive and scope audit

Archived ingredients:

- one-sided partition superadditivity (not merely absolute principal
  monotonicity);
- principal restriction monotonicity;
- random-bridge near-order padding and Theorem 36.15's all-sublinear
  principal heredity;
- switched cut identities and row-field nonnegativity at a maximizer.

I found no archived theorem combining a **global-scale sublinear positive
core**, near-order cap heredity, and both one-sided superadditivity laws to
force

```math
N(A[T])=o(n^(3/2)),
\qquad
N(A[T^c])=M_n-o(n^(3/2)).
```

OC.2--OC.3 are thus a new and useful synthesis.  They narrow the allowed
shape of tail failure to a cut-positive almost-one-sided core paired with
an oppositely oriented near-minimal complement.  They do not prove that
such a decomposition is impossible, do not constrain its cross block, do
not prove scalar packing without the separate selector theorem, and do not
give cross-order recurrence or convergence.

The final classification should therefore say **PROVES A STRUCTURAL
REFINEMENT OF THE OBSTRUCTION**, not "weakens the missing arrow."  No
assumption-distance/frontier reset follows until the one-sided core is
excluded by an exact-minimality theorem.

## 7. Disposition

```text
OC.1: PASS, but archived rather than new.
OC.2: PASS.
OC.3: PASS.
Minimal-size near-clique statement: PASS.
L_balance implication: PASS after sequential and t<c_- repairs.
Novelty: new synthesis, not a new one-sided superadditivity law.
Strict-reduction/frontier claim: REJECT.
```

## 8. Frozen repaired-source recheck

The final source now credits the archived one-sided superadditivity law,
uses the sequential `0<t<c_-` formulation, and classifies OC.2--OC.3 as a
structural refinement rather than a weaker missing arrow.  I rechecked the
repaired file in full.

```text
extremal_information/drafts/exact_minimizer_oriented_core_separation.md
sha256 ef9a8c9125f56919897601ba67a0099705976261b5f1e4d80c674d2c54cc8575
final verdict PASS
```
