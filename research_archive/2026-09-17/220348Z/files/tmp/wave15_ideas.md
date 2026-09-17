# Wave 15: ten ranked continuations after ledger 10.67

1. **Adaptive balanced partition with small total internal norm.**  Prove that
   every global minimizer `A_n`, for a useful sublinear order `s`, has an
   equipartition `V_1,...,V_k` such that
   `sum_i Q(A[V_i])=O(n sqrt(s))+o(n^(3/2))`.  Since
   `C=A-direct_sum_i A[V_i]`, triangle then gives
   `[Q(C)-q_n]_+<=sum_i Q(A[V_i])=o(n^(3/2))`, exactly the new target.  Test
   random balanced partitions and identify the sharp concentration parameter
   if the statement is false as written.

2. **Potts-mask Hanson--Wright theorem.**  Encode a random `k`-coloring by
   centered simplex vectors, so the within-block error is a vector-valued
   quadratic chaos.  Apply or reconstruct Hanson--Wright/decoupling bounds to
   `Q(A hadamard (M-E M))`.  Determine whether the Frobenius term gives the
   desired `n sqrt(s)` and whether the operator-norm term can be removed using
   global minimality rather than the current `||A||_op=O(n^(5/6))` bootstrap.

3. **Bicriteria partition: cross control plus internal-excess capture.**  A
   small cross-only overshoot is useful only if the blocks retain nontrivial
   `sum_i[Q(A[V_i])-q_s]` in the common-mosaic identity.  Form a Lagrange or
   minimax problem over balanced partitions that simultaneously minimizes
   `Q(C)` and maximizes captured internal excess.  Derive the exact scalar
   recurrence this would give for `q_n/n^(3/2)` and seek a probabilistic or
   deterministic rounding theorem.

4. **Recursive dyadic bisection rather than a one-shot partition.**  At each
   node choose a cut whose two induced matrices and cross matrix satisfy a
   centered width inequality, then recurse to order `s`.  Analyze whether the
   errors form a square-summable martingale across levels instead of the
   leading `n s` one-shot triangle bound.

5. **Spectral paving tailored to Boolean norm.**  Use a paving theorem on the
   high-singular subspace of a global minimizer to choose blocks with small
   `||A[V_i]||_op`; then `Q(A[V_i])<=m_i||A[V_i]||_op`.  Audit exactly what
   Marcus--Spielman--Srivastava/Kyng--Luh--Song type results control, whether
   their dependence on `||A||_op` is fatal, and whether the small number of
   high singular values can be handled separately.

6. **Cut-norm/vector-discrepancy coloring.**  Treat every local Boolean
   quadratic state as a discrepancy constraint for vertex coloring.  Replace
   the exponential family by a separation oracle or by the near-ground
   down-set structure already developed in the ledger.  Seek a Banaszczyk or
   partial-coloring argument with total discrepancy `O(n sqrt(s))`.

7. **Amplify the A9 cross-mosaic wall.**  Test lexicographic substitution,
   tensor products, regular blow-ups, and conference products of the `3+6`
   obstruction.  Track global minimality, `Q(C)-q_n`, block order, and the
   hybrid excess exactly.  A leading `Theta(n^(3/2))` amplification would
   falsify the adaptive cross route; a systematic no-amplification identity
   would be useful positive evidence.

8. **Finite partition census of global minimizers.**  Exhaust or sample
   minimizers at the largest feasible orders and enumerate balanced
   partitions.  Record the Pareto frontier between cross overshoot, total
   induced-block norm, internal excess, and Hamming distance to local
   minimizers.  Use it to state the weakest plausible bicriteria theorem,
   while keeping all finite evidence explicitly numerical.

9. **Finite-temperature internal-block control.**  Replace each local maximum
   `Q(A[V_i])` by a log-sum-exp free energy.  Average the latter over random
   partitions using entropy/subadditivity, then take an order-dependent
   inverse temperature so the total max/free-energy error is subleading.
   Check whether this avoids the union bound that blocks direct local-norm
   concentration.

10. **Temporal fallback via cross-Gram/field dichotomy.**  Combine (10.555)
    with the field identity `p^T A n=<f^p,p*n>`: a Hamming-separated endpoint
    pair has either residual or row-field variance.  Try to convert the latter
    into negative centered demand over a contracted order window, explicitly
    testing A5, A7, and A8 before proposing any Hall lower bound.

## Ranking and selected routes

- First: idea 1, because it turns the corrected hybrid frontier into one
  concrete partition inequality and may admit a direct probabilistic proof.
- Second: idea 3, because cross control alone can be tautological (singletons)
  and the exact bicriteria recurrence must be understood before declaring the
  partition route successful.
- Third: idea 7, because an independent amplification search can falsify the
  whole route and the exact A9 seed is already available.

Ideas 2 and 5 are immediate successors if idea 1 isolates an operator-norm
barrier.  Ideas 4, 6, and 9 offer alternative ways around an exponential
state union bound.  Idea 8 supports conjecture selection but is numerical;
idea 10 preserves an independent temporal route.

Primary-source starting point for idea 2: Rudelson--Vershynin,
"Hanson--Wright inequality and sub-gaussian concentration",
https://arxiv.org/abs/1306.2872.  Any imported tail bound must be checked in
the sparse/balanced-color model rather than assumed from the iid scalar form.
