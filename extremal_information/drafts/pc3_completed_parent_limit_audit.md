# Independent audit: PC.3 completed-parent limit

**Verdict.** Pass.  The row law, support normalization, exact-sign
completion, and final `3/4` constant are correct.

The seed relative rows `(a odot b,a odot c)` occur with counts
`(4,8,4,0)` on `(++,+-,-+,--)`.  Tensor coordinates are independent, and
gauging by the anchor gives exactly

```math
(1,X_1,Y_1,\ldots,X_j,Y_j).
```

For

```math
Z=epsilon_0+\sum_t(alpha_tX_t+beta_tY_t),
```

one has `E X=1/2`, `E Y=0`, `E(XY)=-1/2` and

```math
Var(alpha X+beta Y)=7/4-alpha beta\le11/4.
```

Jensen and Cauchy--Schwarz therefore give

```math
1+j/2\le\max E|Z|\le1+j/2+\sqrt{11j/4},
```

so division by `2j+1` tends to `1/4`.  Exact product closure makes the
pre-completion cap

```math
r_jn_j/2+m_jn_j\max E|Z|=(3/4+o(1))r_jn_j.
```

Here `r_j=sqrt(n_j)`, the auxiliary count is at most `r_j`, and any hollow
sign completion costs at most `O(r_j^2)=O(n_j)=o(r_jn_j)`.  Also
`tr(H_j)=0`, so deleting the diagonal preserves every Boolean quadratic
energy.  Finally, total order `n_j+O(sqrt(n_j))` has the same `n_j^(3/2)`
normalization.  No hidden factor of two or diagonal term remains.

This audit verifies a structured signing sequence, not a statement about
the minima at those orders.
