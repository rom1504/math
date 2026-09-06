# Independent audit of the original-class fractional-rounding gain

Date: 2026-09-06, 07:05 UTC. Outcome: PASS for
`resumed_response_original_fractional_rounding_gain_2026_09_06.md`.
The conclusion is a conditional improvement of a feasible spin mean
for a FIXED arbitrary bounded-op hollow signing. It is not a matrix
coefficient edit or a new universal numerical lower certificate.

1. Condition on B, the old randomness, and the mean u BEFORE the fresh
   independent rounding v. With `B=A/sqrt(n-1)`, the row variance is
   exactly `[nd-(1-u_i^2)]/(n-1)`, uniformly `d+O(1/n)`. Hollowness
   makes `v_i` independent of `(Bv)_i`, conditionally on that data.

2. For `f_epsilon(x)=sqrt(x^2+epsilon^2)`, the uniform smoothing error
   is at most epsilon and the third derivative is `O(epsilon^-2)`.
   The sum of third absolute moments of the original row summands and
   matched Gaussian summands is `O(n^-1/2)`: use
   `E|v_j-u_j|^3<=2(1-u_j^2)` and the common coefficient magnitude
   `1/sqrt(n-1)`. Second-order Lindeberg replacement is therefore
   uniform in both u and the shift `(Bu)_i`. Choosing
   `epsilon=n^-1/6` gives the stated `O(n^-1/6)` error. The final
   variance replacement uses the square-root coupling bound and costs
   `O(n^-1/2)`, including when d is zero. No variance-dependent
   Berry--Esseen hypothesis has been smuggled in.

3. If `d>=d_0>0`, then `||Bu||^2/n<=L^2(1-d)`. The already verified
   convexity of `x -> Gamma_sigma(sqrt(x))`, together with monotonicity
   in sigma and in the nonnegative mean argument, proves the uniform
   gap `g_0=Gamma_sqrt(d0)(L sqrt(1-d0))>0`. The normalized signing
   has Frobenius squared n, so any valid operator bound has `L>=1`.
   In particular the fixed damping `eta=g_0/(4L)` is feasible.

4. Exact quadratic expansion gives the retained expected energy gain
   `g_0^2/(8L)-o(1)`. Both the initial rounding and a final independent
   rounding preserve signed quadratic energy exactly in expectation,
   because the matrix is hollow. The assertion remains conditional
   on a pointwise slack bound; positive slack only in expectation
   would not justify the same deterministic g_0.

5. The banked pure terminal endpoints are Boolean, so they have d=0.
   Deliberately scaling a positive-energy Boolean endpoint to create
   slack d costs exactly d times its normalized energy. This particular
   guaranteed increment is at most `d/(4 pi L)`, below that cost at
   the current frontier. This correctly rules out an improvement
   justified by THIS estimate alone, not actual damping improvements
   or a more detailed mixed-field argument.

The source's hard-versus-soft endpoint discussion and its refusal to
deduce uniform slack from vanishing smoothing slack are also correct.
