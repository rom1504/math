# Independent verification of the truncated latent-coreset theorem

**Verdict:** pass.

The conditioning, constants, full/deleted comparison, common-sample
selection, and asymptotics in
[`actual_child_truncated_latent_coreset.md`](actual_child_truncated_latent_coreset.md)
are correct.

On `G_H`, the reverse direction of GC.3 gives
`K_e<=e^(4t)K_0<=e^(4t)H` simultaneously for every edge.  Combining the
pointwise `32K_e/R` empirical-cavity bound with the deterministic squared
error bound four off `G_H` gives (TC.2).  For (TC.3), the relative sample
weight satisfies `W>=e^(-2td)` and has second moment at most `H` on `G_H`;
off `G_H`, both likelihood mixtures lie in
`[(1-rho)^d,(1+rho)^d]`, so their log ratio has magnitude at most `2td`.

At comparable splits, `t^2d=Theta(N)` and `td=Theta(N^(3/2))`.  The chosen
sample size gives good-set cavity error
`O(N^(1/2-zeta)/(log N)^2)`, escaping-set cavity error
`O(N^(1/2-zeta))`, and escaping-set scalar error `O(N^(1-zeta))`; all
remaining scalar terms are `o(N)`.  Markov's inequality converts the old
annealed condition to the new tail condition after a polynomial threshold
enlargement.

The strictness example is valid for ambient scalar random variables.  It
does not establish that the two conditions differ within actual optimizing-
child channels; the audit states that scope explicitly.
