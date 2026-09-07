# Exact exclusion of sparse positive energy at a globally optimal signing

Let the objective be an expected log partition, with arbitrary quenched
randomness z and spin observables chi_e in {-1,1}. The deterministic
interaction on edge e is lambda_e A_e chi_e, lambda_e>=0. Suppose A is
a global minimizing signing (single-edge optimality is sufficient here).
Write R_e(z)=<chi_e>_(A,z). Flipping e gives the exact partition ratio

    Z_(A^e,z)/Z_(A,z)
       =cosh(2lambda_e)-A_e R_e(z)sinh(2lambda_e).

Optimality and Jensen imply

    0 <= E_z log[cosh(2lambda_e)
                       -A_e R_e(z)sinh(2lambda_e)]
       <= log[cosh(2lambda_e)
                       -A_e E_z R_e(z)sinh(2lambda_e)].

Consequently, for lambda_e>0,

    A_e E_z R_e(z) <= tanh(lambda_e).                 (1)

The zero-weight assertion below is automatic and needs no division.
This applies to the original absolute objective, to the width average
by regarding its two equal branches as a discrete quenched variable,
and to a quenched global orientation field. No pointwise choice of the
sign after observing z has been made.

For EVERY edge subset F, (1) yields

    sum_(e in F) [lambda_e A_e E_z R_e(z)]_+
       <= sum_(e in F) lambda_e tanh(lambda_e)
       <= sum_(e in F) lambda_e^2.                   (2)

Thus when lambda_e<=C/sqrt(N), the positive Gibbs energy carried by
o(N^2) edges is o(N), uniformly over the chosen subset. This explains
one precise optimizer property violated by the mesoscopic-clique
integrated-payment counterexample: that example stores a positive
order-N thermal energy in only order-N^(3/2) clique edges.

If the total averaged deterministic energy is nonnegative (as it is
under the homogeneous radial convexity argument used in this campaign),
the negative part is bounded by the positive part. Hence also

    sum_e lambda_e |E_z R_e(z)|
       <= 2 sum_e lambda_e tanh(lambda_e).            (3)

The absolute value in (3) is OUTSIDE the quenched average. Moving it
inside would be an additional, generally unjustified assertion.

## What this does not prove

The estimates refer to the parent's actual optimized Gibbs law. Child
reheating changes that law; its positive energy is not the restriction
of the parent energy. Randomly averaging a cut in (2) therefore does
not prove that the bridge pressure pays reheating. Moreover, negative
energy can still concentrate on a sparse edge set: (2) only bounds its
positive part. Deleting an arbitrary sparse frustrated set can increase
the pressure, so no two-sided sparse-deletion estimate is asserted.
