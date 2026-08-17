# Proof audit: sign-switching visibility synchronization

**Verdict.** PASS.  The result is an exact theorem about hollow sign
matrices.  Its useful new datum is the orbit--character visibility `nu`,
not a renaming of the equivariance norm.

## 1. Constant and normalization checks

For `k` mismatched unordered edges, `A-P^TAP` has two entries of magnitude
two per edge.  Therefore

```math
\|A-P^TAP\|_F^2=8k.
```

The average word length in an elementary abelian group on `s` independent
involutions is `s/2`.  These facts give

```math
\|C-\bar C\|_{op}
\le {1\over2}\sum_i{\sqrt{8k_i}\over r}
=\sum_i{\sqrt{2k_i}\over r},
```

so the coefficient in (SV.10) is correct.

The orbit identity also has the right normalization.  With
`u=z/sqrt(n)`, the pole deficit is `<u,Cu>`, and the character weights sum
to one.  Positivity of `C` makes the twirled diagonal entries nonnegative,
so no cancellation is being used in the division by `p_(jchi)`.

## 2. Logical boundary

Three hypotheses do separate jobs.

* Projective orbit closure lets the proof bound an orbit average using the
  declared individual pole deficits.
* Multiplicity freeness turns the twirl into scalar character defects.
* Visibility prevents a bad scalar defect from hiding in a character that
  every pole orbit barely sees.

Deleting any one invalidates the displayed inference.  In particular,
multiplicity freeness alone does not imply `nu=1`, and exact automorphism
alone does not turn a poorly conditioned orbit frame into a dimension-free
certificate.

The presentation may be singular.  The proof never inverts `G`; it first
proves an operator inequality on `U=range(V)` and then pulls it back.  Thus
the generalized relative inequality remains valid for the four PC.3 poles,
whose Gram matrix has rank three.

## 3. Benchmark audit

For the hollow PC.3 signing the exact characteristic polynomial gives
spectral roof `r=5`, not the Hadamard-completion roof `4`.  Each Boolean
pole has Rayleigh value `4/5`.  The two projective orbit rows see the three
characters with weights `(0,1/2,1/2)` and `(1/2,1/2,0)`, so (SV.11) gives
`2/5`.  Direct character compression gives `(2/5,0,2/5)`, proving
sharpness.

The Cayley benchmark has three spectral values and hence is genuinely
non-Hadamard.  Translation diagonalizes it in Boolean Walsh characters;
there `nu=1` and the theorem reduces exactly to its diagonal Rayleigh
table.

## 4. What is not proved

The Frobenius edit estimate is intentionally severe: at `r~sqrt(n)`, a
vanishing error requires much less than a generic `o(n^2)` edit budget.
No claim is made that near-minimizing signings possess a visible
multiplicity-free switching action.  Establishing such an action, or a
stronger discrepancy estimate for approximate actions, remains the new
structural obligation.
