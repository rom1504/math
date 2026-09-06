# Wave 11: post-`K<=4` frontier

The current theorem is an unconditional endpoint-tree path cover with
congestion at most four.  The missing step is not another local capacity
inequality: it is a scale-sensitive use of that cover for globally minimizing
signings.

## Ten ranked routes

1. **Scale-weighted path minimax.** Replace the uniform path norm in (10.412)
   by weights depending on child order, block excess, or the target decrement.
   Prove a weighted version of the reservoir induction whose integrated
   all-successor cost is summable after `n^{-3/2}` normalization.
2. **Global internal-block replacement.** Simultaneously replace the disjoint
   sibling blocks discarded along an endpoint chain by optimal signings of the
   same orders.  Use global minimality of the parent to force a witness that
   charges the replacement excess to one of the selected successor layers.
3. **Exact laminar replenishment reduction.** Combine (10.340), (10.341), and
   the `K<=4` chain measure line by line and isolate the weakest additional
   inequality that would imply a Cauchy estimate for `q_n/n^{3/2}`.
4. **Sparse multiscale field peeling.** Charge only record-size replenishment
   gaps on geometrically separated batches, using endpoint-chain harvests
   between records so the same `3Q` budget is not reused at every time.
5. **Gibbs information on the endpoint tree.** Average (10.342) over the
   bounded-mass chain measure and seek a data-processing or KL telescope under
   one common reference law, removing the changing-parent-law term that
   stopped the earlier scalar entropy route.
6. **SDP-kernel block exchange.** Show that internal-block optimal replacement
   in a global minimizer forces cross mass through approximate kernels of the
   child SDP slacks, then convert that mass to Boolean layers using the endpoint
   chain rather than an unproved frame hypothesis.
7. **Prime puncture mean stability.** At prime orders `1 mod 4`, use the exact
   endpoint cover of favorable edges to prove the mean child defect
   `bar e_n=O(1)` or the summable variant (10.305), allowing growing affine
   dimension and bounded endpoint reuse.
8. **Guerra--Toninelli-style signing interpolation.** Revisit the canonical
   signing partition function with a smooth interpolation between one system
   and two independent subsystems.  The target must be fixed-temperature
   subadditivity and must explicitly evade the wrong characteristic of the
   already-failed Shearer recursion.
9. **Cut-code boundary-state recursion.** Use the exact short sequence obtained
   by restricting the augmented cut code to two vertex blocks, retaining the
   boundary switching state instead of a scalar covering radius.  Seek a finite
   or controlled-dimensional state whose min-plus product is scale-correct.
10. **Critical-lift minimizer diagnostics.** Locally improve or exchange edges
    in the regular-block critical obstruction and test which global-minimality
    inequalities it violates.  Use it as a falsifier and lemma generator, not
    as evidence that the obstruction itself is globally optimal.

## Selected independent attempts

- Route 1/3: derive the strongest formal scale-weighted or grouped consequence
  of the new path cover.
- Route 2: derive a genuinely non-vacuous internal-block exchange inequality
  for global minimizers.
- Route 8/9: search primary literature for an interpolation or code-recursion
  theorem and attempt an exact adaptation, explicitly checking against the
  ledger's stopped entropy and tensor routes.
