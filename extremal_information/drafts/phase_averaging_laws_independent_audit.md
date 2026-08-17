# Independent audit: logarithmic statistics of an automatic extremal phase

Audited file: `phase_averaging_laws.md`

Verdict: **FAIL pending three local repairs; the theorem statements and proof
formulas are otherwise correct.**

## Repair resolution

The canonicalized version makes `h` an integer at least two, uses a
version-independent explicit trapezoidal sum, and tests one fixed sequence
`q_n=L(t_n)+0.02/sqrt(h^(r(n)))`.  It also restricts the phase-circle language
to the endpoint-compatible case and strengthens the Walsh consequence as
recommended.  The repaired verifier passes; the mathematical audit below
therefore applies without a remaining substantive objection.

## Required repairs

1. **Integer block base.**  PA.1 starts with a real `h>1`, but the proof and
   displayed sums use integer endpoints `n=h^r`.  The intended Hadamard
   application has integer `h`.  State `h` to be an integer at least two, or
   replace every block endpoint by the appropriate ceiling/floor.  With
   integer `h`, the endpoint discrepancies at `s=1` and `s=h` are one term and
   are `o(h^R)` (and `o(h^(alpha R))`), so no formula changes.

2. **The verifier does not run in the repository environment.**  NumPy here
   has `np.trapz` but not `np.trapezoid`; the current script raises
   `AttributeError` before testing anything.  Replace `np.trapezoid` by
   `np.trapz`.

3. **The verifier uses a triangular array, not one sequence.**  Inside the
   outer truncation loop it defines

   ```python
   q = phase(t) + 1.0 / (r + 2) / np.sqrt(base)
   ```

   so the value assigned to a fixed `q_n` changes when the truncation index
   changes.  A direct fixed-sequence test is, for example,

   ```python
   q = phase(t) + 1.0 / np.sqrt(base)
   ```

   which satisfies the blockwise hypothesis with error `h^(-k/2)` on block
   `k`.

## Proof checks

- On block `h^r <= n < h^(r+1)`, the harmonic sum is the Riemann sum for
  `int_1^h psi(L(t))dt/t`.  Uniform response error gives `o(1)` per late
  block; Cesaro averaging over the `log_h N+O(1)` blocks proves PA.2.  A final
  partial block has bounded harmonic mass and therefore vanishes after
  division by `log N`.
- The complete unweighted blocks have masses `h^r`; summing them gives the
  coefficient `1/(h-1)` in PA.4.  The extra endpoint when `s=h` is `O(1)` and
  is harmless.
- For `alpha>0`, the complete weighted blocks have masses `h^(alpha r)`;
  geometric summation gives `1/(h^alpha-1)`, and
  `sum_(n<=N)n^(alpha-1)=N^alpha/alpha+o(N^alpha)`.  Differentiation gives
  `C_alpha'=(alpha/s)(L-C_alpha)`, so constancy of `C_alpha` is equivalent to
  constancy of `L`.
- Under `u=log_h t`, `dt/(t log h)=du`; hence `cos(2 pi u)` and
  `cos(4 pi u)` have the same pushforward arcsine law.  PA.11 correctly shows
  that their positive-power response phases remain distinct.
- The Walsh consequence follows from Theorem 30.1 and Corollary 30.2.  In
  fact PA.3 makes the ordinary and every positive-power subsequential phase
  nonconstant, so “may retain” can be strengthened to “retains.”

## One interpretation caveat

For a general profile in PA.1, `L(1)=L(h)` was not assumed.  Thus
`dt/(t log h)` is always Haar measure on a fundamental logarithmic interval,
but `L` defines a continuous function on the quotient phase circle only when
its endpoints agree.  The Walsh profile does have equal endpoints, so the
circle language is valid in the intended application.
