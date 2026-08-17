# Independent audit: Gram-relative coherence blindness

**Verdict.** Pass at the symmetric weighted-contraction level.  This is not
an exact-sign child or a separated-optimum theorem.

For odd `p`, the normalized majority tail vector

```math
u_S={hat(Maj_p)(S)\over\sqrt{rho_p}}1_(|S|>=3),
\qquad rho_p\longrightarrow1-2/pi
```

is unit, and every coordinate square is `O(1/p)`.  On the orthogonal odd-
character pole basis, both `D_coh=uu^T` and
`D_diag=diag(u_S^2)` are PSD contractions.  Hence `R=I-D` is realized by a
symmetric PSD contraction on the row space.

The systems agree on every generator row/column, every individual product
deficit, trace, and average defect.  Yet their relative defects are `1` and
`O(1/p)`, and the prescribed majority selector losses tend respectively to
`1-2/pi` and zero.  With `N=2^(2p)` and `r=sqrt N`, this is a fixed
`Theta(N^(3/2))` separation of the selector certificates at
`q=sqrt(N)/2` active products.

The twirling theorem is also correct.  If every diagonal deficit is at most
`d`, then

```math
D\preceq(d+||A-Abar||_op)I.
```

For commuting involutory generators,

```math
||A-Abar||_op\le {1\over2}\sum_i||[A,rho(g_i)]||_op.
```

The factor `1/2` follows by averaging the Hamming length of a uniform group
word.  In the coherent example the equivariance defect is `1-O(1/p)`, so
the criterion detects rather than hides the obstruction.

The claim is deliberately limited to the full generator block and separate
product deficits; arbitrary low-degree cross-product data were not matched.
Another Boolean spin may repair the prescribed selector, and the realized
operator is weighted rather than hollow and entrywise signed.
