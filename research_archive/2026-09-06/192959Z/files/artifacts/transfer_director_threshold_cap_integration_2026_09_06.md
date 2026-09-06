# Direct integration of discrete-threshold entropy into the all-order upper proof

2026-09-06. **Proved, conditional only on the explicitly cited already-audited
annealed construction certificate.** This changes neither that certificate
nor its fixed-depth/all-order realization proof.

Let `p=31/32`, `t_weave=4`, `delta=2^-24`, `q=1-2delta`,
`v=4delta(1-delta)`, and

```math
 a_* = \frac{91470529542342299}{20460000000000000000},
 \qquad g_* = \frac{18725}{32768}v^2.
```

The threshold/resampling theorem proves an actual-signing partition lower
bound with credit at least `g_*` per vertex whenever
`||A||op/sqrt N<=33/32` and `8<=beta<=33/4`. Its errors are
`O(sqrt N)+log(N+1)+O(1)` with fixed parameters. For the retained weave,
`N=p m^2+O(m)`, `beta=8/sqrt p+o(1)`, and
`||A||op/sqrt N<=1/sqrt p+o(1)`. The displayed ranges thus hold for
all sufficiently large construction orders, with strict margins.

For every fixed `a'<a_*`, the existing Gaussian-boundary/conditional-variance
certificate permits a fixed finite depth with annealed pressure at most
`(4sqrt p-a'+o(1))m^2`. Applying Markov's
inequality exactly as in the standalone upper proof and then using the
deterministic threshold-law lower bound, and finally sending `a'` up to
`a_*`, gives

```math
 \limsup_{n\to\infty}\frac{M_n}{n^{3/2}}
 \le C_{\rm thr}:=\frac1{q^2}
 \left\{\frac12-\frac{a_*+p[h(\delta)+g_*]}{8\sqrt p}\right\}.
```

All-order realization uses the same relatively dense H2/H12 construction
orders and principal deletion. The extra threshold errors vanish after
division by `m^2`; they introduce no dependence on target residues or on
an unknown target minimizer. Keep the existing limit order: fixed strict
certificate margin, fixed sufficient recursion depth, construction order
tending to infinity, then margin tending to zero. The threshold parameters
are fixed throughout. No new interchange of these limits is used.

The clipped-law predecessor replaces `g_*` by
`g_clip=beta^2 v^2/128=v^2/(2p)`. Consequently its upper coefficient
strictly exceeds `C_thr`, by

```math
 \frac{\sqrt p}{8q^2}(g_*-g_{\rm clip})>0.
```

This small numerical gain is not the main result: the entropy bound is
uniform for actual bounded-operator sign matrices and does not substitute
continuous Gaussian KL for the discarded discrete information.
The safely displayed interval remains
`[0.4333221116640807,0.499432211)`.

## Proof dependencies and replay

- `transfer_reconstruction_standalone_2026_09_06.md`: construction,
  Gaussian replacement, conditional-variance supersolution, full-spin
  count, exact annealed certificate and all-order realization.
- `transfer_seed_threshold_entropy_resampling_2026_09_06.md`: exact
  discrete entropy bound, rational coefficient 22, pressure remainder,
  and credit `18725 v^2/32768`.
- `transfer_reconstruction_correlated_cluster_credit_2026_09_06.md`:
  predecessor only; used to compare the old coefficient, not to prove
  the new threshold lower bound.
- `computations/transfer_director_threshold_cap_integration_2026_09_06.py`:
  exact rational intervals, importing the tracked elementary log/square-root
  interval routines and the new rational threshold checker. Its generated
  JSON distinguishes exact intervals from illustrative decimals.

The original convergence question remains open. This theorem supplies an
all-order construction bound, not propagation of the unknown optimum.
