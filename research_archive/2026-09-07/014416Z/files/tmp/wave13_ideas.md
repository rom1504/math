# Wave 13: grouped compatibility after decrement tolling

Wave 12 removed the leading Q(A) term from the endpoint harvest, but only
after tolling each bucket by its edge decrement.  The missing theorem must
couple centered temporal demand to the residual capacities, localize negative
credits, or replace all diagonal blocks jointly with a near-minimizing hybrid.

## Ten ranked routes

1. **Endpoint-tie minimax for residual service.**  Jointly choose the
   positive/negative endpoint pair, the Section 10.62 allocation, and a
   temporal-to-bucket transport.  Dualize the max over endpoint ties and test
   whether averaging endpoint pairs defeats the prescribed-tree A5 wall.
2. **Causal signed Hall flow.**  Atomize every negative centered increment and
   terminal excess at its temporal location.  Permit a credit to pay only
   descendants or events in the same order window, derive the exact
   ancestor/prefix cut dual, and test whether telescoping forces those cuts.
3. **Residual-capacity exchange identity.**  Express
   `[c_b-(Q(U)-Q(X))]_+` as a boundary functional or global replacement
   excess, then seek a conservation law pairing zero residual with a negative
   temporal credit elsewhere in the same endpoint mosaic.
4. **Small-sibling endpoint stopping.**  Choose endpoint ties so that all but
   logarithmically many relevant splits peel siblings of size
   `o(sqrt(n))`; handle the exceptional macroscopic splits by an exact
   separate potential instead of charging them to (10.515).
5. **Make temporal peeling follow the endpoint tree.**  Rather than transport
   between unrelated trees, randomize which endpoint shore survives and
   choose its orientation so that the temporal decrement/replenishment
   identity and the residual allocation live on the same edge.
6. **Joint global hybrid minimax.**  Partition into mesoscopic blocks and
   choose all block minimizers simultaneously to minimize the hybrid excess
   `Q(tilde A)-q_n` in (10.524).  Use switching/permutation randomization only
   after restricting to the original minimizer's low-slack mosaic states.
7. **Hybrid block-size optimization.**  Quantify the random-switching hybrid
   excess for block size `n^alpha`, balance diagonal improvement against
   profile entropy, and determine whether any alpha gives `o(n^(3/2))`.
8. **Selected-child excess potential.**  Define the Bellman value of the best
   macroscopic induced child and seek a vanishing potential V satisfying
   (10.531), using near-top bound (10.527) only to land between macro steps.
9. **Cut-code repair for common mosaics.**  Reuse the compatible-child
   codeword repair identities around (10.307)-(10.308) to bound the global
   hybrid excess by disagreement among block closest-codeword profiles.
10. **Finite causal-transport dual from the literature.**  Reconstruct the
    filtration-respecting dual in Lassalle's causal transport framework and
    specialize it to the finite temporal/endpoint trees; keep it only if it
    yields a stronger potential than the elementary prefix-flow dual.

## Selected independent attempts

- Routes 1/3: optimize endpoint ties and residual service jointly, with exact
  finite minimax certificates and a candidate universal inequality.
- Routes 2/10: derive the honest causal signed-flow dual and prove or falsify
  its ancestor/window cut inequalities from the centered peeling identities.
- Routes 6/7: attack the joint global hybrid excess, beginning with an exact
  variational formula and finite/randomized falsification before any
  asymptotic claim.
