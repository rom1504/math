# Conditional-exposure ideas retained after their limitations were identified

Date: 2026-09-06. Status: exploratory calculations, NOT additional theorems.
These notes preserve abandoned or incomplete lines from the director's
working analysis; they are not a new active research state.

For an orthogonal weave partition W=[[K,L],[L^T,J]], the exact identities
KL=-LJ and K^2=m^2I-LL^T put the nonendpoint spectrum of K in range L.
This suggested exposing the complement to choose a small candidate net
before an energy tail calculation. The following attempted extensions did
not establish the required tail.

- A small norm for the deleted block J would separate middle eigenvalues
  from a candidate mean field near the retained extremal value. But an
  L2 field-variance gap alone cannot yield an adaptive product-entropy
  improvement: the explicit hub example in the adaptive-cluster proof has
  nonvanishing field variance and vanishing L1 field dispersion.
- The middle spectral space contains actual Boolean columns of L. Its
  small dimension alone therefore does not exclude relevant spins.
- Revealing (L,J) reconstructs K exactly, and for odd deleted-fibre size
  L alone fixes it up to global sign. Those facts are now proved in
  `transfer_adversary_complement_exposure_2026_09_06.md`.
- One could try unioning candidate nets over all remaining outer seeds.
  The seed count costs roughly exp(m^2 log2/2), potentially less than
  a full-spin union. However the required uniform conditional m^2-speed
  energy tail is false for arbitrary exposed bases: the exact seed-only
  example costs only exp(-O(m)). It occurs at a spectral endpoint, so
  a separately proved interior-mean restriction would still be needed.
- Writing (B+mu)e=(s^2-mu^2)x for a nearly eigenvector-like Boolean spin
  was considered as a route from L2 to L1 field control. A residual
  concentrated on a small coordinate set can interact with a large-norm
  principal submatrix. No uniform inequality eliminating this possibility
  was proved, and no independence may be assumed for such a residual.

These are limitations of particular estimates or exposures, not proof that
all local-field or conditional-net methods fail. The corresponding exact
counterexamples and their parameters remain in their canonical artifacts.
