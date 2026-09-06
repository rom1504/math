# External inverse-theorem scale audit

Status: verified translation of the stated primary-source theorems; this is a
route audit, not a new theorem about `q_n`.

Primary sources:

- Balla--Hambardzumyan--Tomon, *Factorization norms and an inverse theorem for
  MaxCut*, arXiv:2506.23989, Theorems 1.1 and 1.7.
- Alon--Makarychev--Makarychev--Naor, *Quadratic Forms on Graphs*, Invent.
  Math. 163 (2006), author PDF.

The MaxCut inverse theorem says that a graph with `m` edges and MaxCut at most

```math
m/2+\alpha\sqrt m
```

contains a clique of order `2^{-O(\alpha^9)}\sqrt m`.  A competitive signing
only bounds its signed cut deviations at order `n^{3/2}`.  More explicitly,
let `G` consist of the negative edges, put `N=\binom n2`, `m=|E(G)|`, and
`T=\sum_{i<j}a_{ij}=N-2m`.  For a cut `S\sqcup S^c`, write
`C_A(S)=\sum_{S\times S^c}a_{ij}`.  Then

```math
|T|\le Q(A)/2,\qquad |C_A(S)|\le Q(A)/2,
```

and

```math
e_G(S,S^c)-m/2
=\frac{2|S||S^c|-N}{4}+\frac T4-\frac{C_A(S)}2
\le\frac n8+\frac{3Q(A)}8.
```

For a competitive signing, `m=\Theta(n^2)` and this corresponds to
`\alpha=O(\sqrt n)`, so the asserted
clique size becomes

```math
n\,2^{-O(n^{9/2})},
```

which is vacuous.  The useful regime of that theorem is cut surplus `O(n)`,
one square-root below the present scale.

The factorization-norm rectangle theorem has the same mismatch.  For an
`n`-by-`n` sign matrix, the elementary factorization through the identity
gives `\gamma_2(A)\le\sqrt n`; converting to a Boolean matrix changes this by
at most an absolute affine term.  The guaranteed homogeneous rectangle has
relative dimensions only

```math
2^{-O(\gamma_2(A)^3)}=2^{-O(n^{3/2})}.
```

No implication from `Q(A)=O(n^{3/2})` to bounded `\gamma_2(A)` is supplied by
the paper.  Therefore its current quantitative statements do not yield a
macroscopic block, response profile, or recursive partition here.

Finally, the graph Grothendieck framework does apply formally to the complete
support, but its complete-graph integrality gap grows logarithmically.  A
generic SDP/vector relaxation would therefore lose a non-negligible factor,
not the `1+o(1)` control required for convergence.  Either external route
would need additional structure specific to globally minimizing signings.
