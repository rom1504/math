# Independent verification of the diffuse collision falsifier

**Verdict:** pass.

The determinant MGF, projective lift, entropy-variational constants,
collision comparison, endpoint extraction from (37.184), and asymptotic
coefficient in
[`actual_child_latent_collision_generic_falsifier.md`](actual_child_latent_collision_generic_falsifier.md)
are correct.

In particular, the determinant in LC.4 has the right `n by n` dimension;
the factor-gauge lift preserves KL exactly; and at
`theta=1/(2sqrt(d))` the MGF cost is at most `||A||_F^2/(6d)`.  Thus the
coefficient in LC.5 is at least `1/(3sqrt(d))`.  The insertion estimate
gives `||M||^2>=||r||^2/2-4rho^2d`, Jensen upgrades pointwise `D_2` to the
annealed collision mean, and `K_e>=e^(-4t)K_0` gives the deleted result.

The proof of (37.184), rather than only its actual-child headline, confirms
that it applies to every rank-one prior.  Taking `delta=lambda` is valid for
the pointwise endpoint estimate when `lambda<1` and
`beta>beta_BG/(1-lambda)`; the resulting coefficient
`sqrt(gamma_0)eta_*/6` is exact.  The audit correctly labels the uniform
prior nonactual and does not claim an optimizer-specific obstruction.
