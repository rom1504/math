# A bounded adjacent test: move the Gaussian value inside the latent expectation

**Subsequent resolution in the same campaign:** the proposed supersolution
and a separate negative upper certificate were proved and independently
audited. See `continued_director_strict_all_order_upper_2026_09_06.md`.
The text below preserves the motivating test and its original evidentiary
distinctions; its formerly open factorization is no longer an obligation.

Date: 2026-09-06. The director's last bounded test of the missing Bellman
supersolution. The comparison inequalities below are proved; factorization
is OPEN. No improved signing bound follows.

With the existing g_t, Phi_t and latent envelope E_t, consider only as a
candidate majorant

`Ebar_t(nu)=sup_L [E_L g_t(Var(X|L))-I(X;L)]`.

No new dynamic state is introduced. Its purpose is to test whether a
classical conditional-envelope factorization is available when the original
g_t(E Var) expression has heterogeneous-temperature difficulties.

## Proved comparisons and exact sufficient implication

Convexity of g_t gives E_t<=Ebar_t. Also Ebar_t<=Phi_t. For the latter,
choose optimal conditional self-couplings given the finite label L. Their
mixture has the original source marginals; the mutual-information chain
rule gives

`F_t(nu)<=I(X;L)+E_L F_t(nu_L)`.

Thus `Phi_t(nu)>=E_L Phi_t(nu_L)-I(X;L)/2`, whereas
`g_t(Var(nu_L))<=Phi_t(nu_L)`. Subtracting I(X;L), rather than half of it,
proves the asserted upper comparison. Approximate optimizers suffice;
finite-second-moment conditional laws meet all self-transport hypotheses.

If this candidate satisfies `B Ebar_t<=Ebar_t` on finite reachable source
laws, the already proved stopped-tree argument yields

`lim_r B^r Phi_t(nu)<=Ebar_t(nu)`.

Indeed use Ebar as the stopped supersolution but stop on the compact gap
Phi-E, exactly as in Section6 of the terminal-gap theorem. Since Ebar>=E,
all stopping errors are still controlled. There is no need to prove that
Ebar equals the Bellman fixed point, or that it is a subsolution.
Consequently a certified negative offset value of Ebar at the ternary
source would suffice for a strict all-order upper improvement. Both the
supersolution and that numerical upper certificate remain OPEN here.

## Exact finite-source interpretation and ternary reduction

For a fixed finite source alphabet,

`H(nu)+Ebar_t(nu)=cav_nu [H(nu)+g_t(Var(nu))]`.

This is the ordinary concave envelope on its probability simplex: latent
posteriors have barycentre nu, and mutual information is source entropy
minus expected posterior entropy. It differs from the original E_t because
the Gaussian function now acts separately on each posterior variance.

For `nu_p=(1-p)delta_0+(p/2)delta_(+1/sqrt(p))`
` +(p/2)delta_(-1/sqrt(p))`, reflection symmetrization gives posterior
parameters z,s in [0,1], variance `(z-z^2 s^2)/p`, and posterior entropy
`h(z)+z h((1+s)/2)`. The sole mixing constraint is E z=p. Therefore the
offset criterion is EXACTLY

`-h(p)+t(1-sqrt(p))+cav_z max_s`
` [h(z)+z h((1+s)/2)+g_t((z-z^2 s^2)/p)] (p)`.

A two-point z mixture suffices by the one-moment extreme-point argument.
A 2001-by-2001 grid gives a LOWER witness about -.01505539 at p31/32,t4,
with z support approximately .9085,.997. This is not a rigorous upper
bound and does not rule out a positive unsampled optimum. The reproducible
script prints its mesh and scope; floating transcendental values are not
certificate evidence.

## Primary literature: mechanism, not an imported solution

The director and convergence researcher read Nair's
[Upper concave envelopes and auxiliary random variables](https://chandra.ie.cuhk.edu.hk/pub/papers/manuscripts/concenve.pdf),
especially the full proofs in Section4.1. Those results concern fixed
product-channel information functionals. Their chain-rule factorization
uses common coefficients and a favorable remaining mutual-information
term. No exact identification of the present variance functional with
one of those functionals has been established. Section4.2 of that paper
itself distinguishes a more complicated factorization conjecture from
its proved cases. The name 'concave envelope' therefore supplies no proof
of the open Bellman inequality.

The candidate should be rejected immediately if a finite Bellman pair
violates its supersolution property at a point relevant to the construction,
or if its full root offset is nonnegative. It is an adjacent falsifiable
test, not a replacement of the proved Gaussian-boundary recursion.
