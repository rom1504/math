# Wave 46 idea ranking

Read before selection: `README.md`, Wave 45 `Updated frontier`
(ledger §10.98.4), and the early Wave 46 `STEERING.md` refresh.  Exclude
average incidence partitions, central half-cube Harnack, nearest-core hard
regeneration, strict pressure/high replicas, local harmonic Poincare, and
affine phase geometry without crossing-overlap control.

## Ranked ideas

1. **Noncentral block switching with explicit slack.**  Keep the exact
   transport
   `L_S(d^U)-L_S(d)=4C_d(U)-8C_d^S(U)` and condition on an active
   selector's slack `Gamma_S=L_S(d)-q_n`.  Derive a fixed-layer
   survival lower bound in terms of `Gamma_S`, `R_(B_S)(d)`, and
   `k/n`, then optimize the Harnack gain divided by
   `1-vartheta_k`.  The goal is a sharp slack-versus-row dichotomy in
   a genuinely noncentral regime, or a scalable obstruction showing that the
   remaining block route is also too strong.

2. **Hard common-core boundary identity and improved extraction.**  Rewrite
   `rho_ell` as one minus the incidence-size-biased escape probability
   of the hard favorable families under `K_ell`.  Use the exact spectral
   decomposition to sharpen (10.1176), identify the precise boundary deficit
   that must beat `1-lambda_ell`, and test whether near-identity
   `ell=m-Theta(n/L0)` can turn a clipped slack/coarea estimate into the
   needed degree.  Do not treat kernel iteration as new evidence.

3. **Tropical tight-principal-decomposition theorem.**  At low temperature,
   derive the endpoint exponent exactly.  The selector normalizer cancels
   inside the lifted likelihood, leaving `E_d-min_S K_S(d_S)` up to a global
   constant, where `K_S` is the best external completion energy.
   Affineness is
   equivalent to a full parent ground state that is simultaneously a lift of
   a ground child for a maximizing selector and an optimal external
   completion.  Prove this from exact minimality, or find an exact-minimizer
   counterexample; then quantify how its active crossing degree relates to
   selector exclusion.

4. **Mixtures of common-core kernels.**  Optimize a nonnegative mixture of
   `K_ell` over core scales.  Its retention and Johnson eigenvalues are
   linear in the mixture.  Formulate the exact finite LP and determine whether
   it yields a useful spectral excess when every nearest-core gap is negative,
   or whether the optimizer collapses to a single scale and adds no mechanism.

5. **Johnson coarea for clipped deficit.**  Express the hard boundary of
   `{S:d_S(z)<=H}` as a level-set integral of the clipped deficit.
   Seek a minimizer-specific inequality charging every escaping replacement
   to available cap slack or a row influence.  The required estimate must be
   averaged under hard incidence and reach the `1/L0` escape scale.

6. **Robust-incidence cloud from complement slack.**  For each covered
   selector, count block switches or selector replacements that retain the
   same complement witness because `Gamma_S` absorbs the exact transport
   error.  Double-count robust pairs to obtain one large row-good column.
   Ordinary polynomial Hamming balls and unconditioned first moments are
   insufficient; the cloud must have stretched-exponential reach.

7. **Scalar lower envelope across the row price.**  Study the Pareto curve
   `min_d{h_d+lambda R_2(d)}` as `lambda` varies.  Integrate its
   exact slopes and transition jumps, using complement cover constraints, to
   seek a price at which both mass and row are controlled.  First determine
   whether the state-count dilution (10.1170) makes every cover-only argument
   impossible.

8. **Weighted harmonic overlap from exclusion mass.**  The exact identity
   `sum_i tilde B_i=(n-m)f` constrains posterior selector exclusion.
   Try to charge every positive low-temperature base-surprise crossing to an
   exclusion event and bound the weighted overlap constant `D` without
   bounding the raw degree `Delta`.  Keep the zero-leading-slope and
   orientation costs explicit.

9. **Direct bare-tail count conditioned on principal norm.**  Couple the
   row-sensitive signed slice tail to the lower tail of `Q(A[S])` at the
   hard event, rather than through a soft center pressure.  Seek a
   dependence-sensitive lower count compatible with (10.1022); ordinary
   moments and fixed-center margins remain excluded.

10. **Adaptive favorable-Gram witness.**  Revisit the exact positive external
    Gram identity (10.1078) only with a joint selector rule that chooses a
    child ground and a low internal Gram witness simultaneously.  A useful
    theorem must prove the `n^(9/4-c)` internal scale and cap persistence;
    independent noise and arbitrary grounds are already obstructed.

## Wave 46 assignments

- Assignment A: ideas 1 and 6, noncentral complement slack transport.
- Assignment B: ideas 2, 4, and 5, hard common-core escape and spectral excess.
- Assignment C: ideas 3 and 8, tropical tight decomposition and harmonic
  crossing geometry.
