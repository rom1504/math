# Wave 31 idea ranking

Read before ranking: `README.md`, Wave 30 `Updated frontier`, and the Wave 30
version of `STEERING.md`.  The required scale is

```math
L=n^{3/4-c},\qquad k=\Theta(L/\log n),\qquad
r=\Theta(n^{1/4+c}\log n).
```

The random-label wall rules out arguments based only on one favorable local
completion per selector.  Pair/triple gluing and control of every refinement
are not open targets.

## Ranked ideas

1. **Exact-minimizer exchange control of batch signatures and `C_sig`
   (selected).** Compare favorable principal completions on selector
   exchanges, sum the exact ground-state inequalities around paths/cycles,
   and seek a batch theorem bounding both coordinate-signature entropy and
   the coarsest hybrid row field.  The deliverable is either the Wave 30
   completion lemma or a precise unavoidable missing term/counterexample.

2. **Bare collision moment/dual contradiction (selected).** Work directly
   with `J_r`, the hit incidence matrix, and an arbitrary selector law.  Seek
   a lower bound from higher codegrees, fractional moments, Hölder/Finner,
   or the KKT structure of a bad cover dual, without preselecting child
   completions.  State the weakest extra incidence inequality that closes
   `Pr(J_r)>=e^{-O(rL)}` and test whether pointwise planting can imply it.

3. **Sequential harmonic composition at fixed density (selected).** Derive
   how the exact matched harmonic cost `Phi_beta` composes when outside spins
   are integrated one at a time.  Test Doob/Efron--Stein, conditional
   variance, and entropy-chain bounds for an `O(n^{1/2-2c})` all-`k` result;
   otherwise produce a rigorous accumulation wall stronger than the existing
   statement that naive one-deletion summation is too costly.

4. **Adversarial-law entropy dichotomy.** Heavy atoms already give coverage
   by planting; for high-min-entropy laws, seek a Johnson-slice
   isoperimetric or dependent-random-choice theorem forcing a mesoscopic
   subbatch with usable overlap.  Explicitly test sparse constant-distance
   supports, which may defeat naive expansion.

5. **Aggregate multi-selector SDP and low-signature rounding.** Formulate a
   common coset as `r` favorable cut columns whose coordinate rows have at
   most `O(k)` projective patterns.  Relax row clustering through a Gram/SDP
   model and try Grothendieck or vector-discrepancy rounding while retaining
   the coarsest `A^2` cap.

6. **Energy-sensitive rate-distortion clustering.** Use the positive
   `B_(n,m)+t` deficit allowance to alter exact child completions before
   grouping coordinate signatures.  Replace worst-case Hamming Lipschitz
   loss by signed local-field costs and seek `O(k)` prototypes at total
   admissible energy loss.

7. **Parent-Gibbs batch completions.** Sample favorable completions from a
   common parent Gibbs prior and bound the entropy of their coordinate
   signature process plus `C_sig` by a matrix MGF.  This is distinct from the
   falsified pointwise low-information argument only if it proves a joint
   mesoscopic collision estimate.

8. **Regularize a bad dual law.** Show that a selector law witnessing small
   coset coverage can be smoothed, with controlled loss, into a low-degree or
   bounded-marginal-density law.  Slice harmonic analysis might then attack
   the structured adversary.  A valid result must not replace the bad law by
   the uniform slice, whose mean is circular.

9. **Direct hybrid-row inequality.** Express `C_sig` as the maximum of
   `x^T A^2 x` over hybrid spins made from completion signature cells.  Seek
   an exact-minimality inequality improving the generic
   `2(n-1)q_n` cap by `n^{1/4+c}` for such hybrids, perhaps via signed
   perturbations of `A` or the completion deficit inequalities.

10. **Multiscale harmonic deletion.** Control conditional partition sums for
    mesoscopic outside blocks where a log-determinant or quadratic-Rademacher
    MGF is usable, then compose blocks rather than individual deletions.
    Audit whether accumulated costs can be summable over a fixed-ratio
    landing before attempting a cluster expansion.

## Assignments

- Route 1: exact-minimizer exchange/signature theorem.
- Route 2: bare collision moment and adversarial dual.
- Route 3: fixed-density harmonic composition.
- Main agent: aggregate multi-selector low-signature relaxation (Idea 5),
  falling back to an exact codegree criterion if the relaxation is circular.
