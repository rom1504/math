# Observable freeze: finite radial response-roof audit

Date frozen: 2026-08-17, before computing the response data.

Status: experimental protocol.  No numerical outcome was inspected when these
observables and comparison classes were fixed.

## Population

For each feasible order `7 <= n <= 11`, use three separately labelled classes:

1. repository matrices certified by their stored cap to be exact minimizers;
2. distinct one-edge perturbations of an exact minimizer, retaining only cap
   `M_n+2` matrices (one-step near-minimizers);
3. fresh uniformly random hollow symmetric sign matrices generated from seed
   `20260817`, without conditioning on their cap.

Deduplicate byte-identical matrices.  The target sample is at most four
matrices per class and order.  Selection within a class is deterministic and
does not use any response observable defined below.

## Exact contexts and response roof

Let `E=binom([n],2)`, let `a^F` reverse the edges in `F`, and let

```math
\mathcal B_r=\{F\subseteq E:|F|\le r\},\qquad r\in\{1,2,3\}.
```

The projective spin convention is `x_0=1`.  An augmented cut witness is
`z=(x,sigma)`, with edge coordinates `z_ij=sigma x_i x_j`.  Define

```math
g_z(F)=\langle a^F,z\rangle-Q(a),
\qquad R_a(F)=\max_z g_z(F)=Q(a^F)-Q(a).
```

All spins, all edit contexts, all ties, and all quantities below are
enumerated exactly.  Energies and errors are stored in unnormalized integer
units.

## Frozen observables

For each radius, record:

1. **RS shell size.**  With `eta=Q(a)-M_n`, count the augmented witnesses in
   `S_a(eta+2r)={z:Q(a)-<a,z><=eta+2r}`.  This is the common shell guaranteed
   by the radial theorem, not an asserted minimal state.
2. **Exposed witness count.**  Count witnesses attaining `R_a(F)` for at least
   one `F in B_r`.
3. **Exact and approximate response-cover number.**  For integer tolerance
   `delta in {0,2,4}`, let `C_delta` be the smallest subset `T` of exposed
   witnesses such that

   ```math
   \max_{z\in T}g_z(F)\ge R_a(F)-\delta
   \quad\hbox{for every }F\in\mathcal B_r.
   ```

   This is solved as a finite set-cover integer program.  Report the exact
   optimum when proved; otherwise report certified lower and feasible upper
   bounds and the solver status.
4. **Tie-invariant optimizer information.**  Draw `F` uniformly from
   `B_r`, then draw `Z` uniformly from the exact optimizer set at `F`.  Record
   `H(Z)`, `H(Z|F)=E log_2 |Argmax(F)|`, and
   `I(F;Z)=H(Z)-H(Z|F)`, together with effective support `2^{H(Z)}`.
5. **Response distribution.**  Record the histogram and Shannon entropy of
   the integer response `R_a(F)`, the mean optimizer multiplicity, and the
   maximum multiplicity.
6. **Contextual affine-function metric.**  On exposed witnesses use

   ```math
   d_r(z,z')=\max_{F\in\mathcal B_r}|g_z(F)-g_{z'}(F)|.
   ```

   Record its diameter and deterministic greedy packing lower bounds and
   covering upper bounds at tolerances `2` and `4`.  These greedy quantities
   are explicitly not claimed optimal.

Primary normalized comparisons are `log_2(size)/n`, `I(F;Z)/n`, and response
entropy divided by `log_2 |B_r|`.  No finite separation is interpreted as an
asymptotic theorem.

## Pre-registered discriminators and falsifiers

The compressibility hypothesis would receive finite support if exact and
one-step-near matrices had materially smaller normalized exposed/cover counts
and optimizer information than unconditioned random matrices, consistently
across orders and radii.

It is falsified at these orders if near-minimizers have comparable or larger
normalized response-cover complexity, or if the exact cover is already a
fixed fraction of the guaranteed shell.  Either outcome is only finite
evidence.  A solver timeout is not evidence in either direction.
