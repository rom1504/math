# Independent audit: bounded-op cubic projection and direct input comparison

Date: 2026-09-05. Auditor: literature track. No Gaussian replacement of the
matrix entries is used: the only replacement is of the independent auxiliary
spins in a fixed polynomial response to the original deterministic signing.

## 1. Audit conclusion

The claims and proofs in
`fresh_limit_unmarked_cubic_projection_2026_09_05.md`, Sections 1--5, pass this
independent audit. In particular, the hypotheses are only a fixed bound on
`||B||op`, the original hollow symmetric sign entries, and a fixed finite
family of injective fully marked odd-degree tree fields. No off-diagonal
coherence assumption on `B^2` or higher powers is silently needed.

The conclusion is a weighted bilinear projection, not yet a lower-bound
improvement. A local weight depending on the old fields cannot be inserted
in the projection identity without further proof.

Here are the delicate checks.

* In the arbitrary-cut flattening proof, removing the external root leaves
  a connected tree. Every internal-edge component meets a crossing edge;
  therefore its crossing-isolated vertex count is at most its internal
  edge count. The copy maps used on crossing-incident vertices are
  isometries, while every genuinely isolated vertex contributes exactly
  the stated `sqrt(n)` factor. There is no missing exponential-in-n factor.
* The injectivity correction has squared Hilbert norm `O(1/n)` because
  forbidden tuples number `O(n^(d-1))` and each coefficient is
  `O(n^(-d/2))`. Symmetrization involves only a fixed number of slots.
* The repeated-slot formula for the transported cubic is exactly
  `K_i(a,a,b)=(Q_ib-B_ia B_ab)/(n-1)`. Its summed square is `O_L(1/n)`.
  This correction is essential: deleting repeated coordinates independently
  in every row and then summing absolute errors would not suffice.
* For the degree-three marked tree, normalization gives the covariance
  `sqrt(3) (BC)_ij` with `h3(G_j)`, where `C=B circ Q^(circ 2)`.
  The weighted sum is a diagonal entry of `BCQ`, not `B^3 C`.
  The row-sum bound `||C||op <= L^2/sqrt(n-1)` controls both the old/new
  field covariance and the exceptional third-derivative projection term.
* All first- and second-chain-rule classes in Section 4 have either at
  least two small non-edge derivatives, or an additional factor of
  `Q_ij`. The latter permits the squared-row-sum bound instead of the
  potentially expensive absolute sum. The third-chain-rule class must
  sum over `j` before taking norms, exactly as done there.
* The variance bounds are `0 <= v_i <= L^2`, not pointwise `v_i >= 1`.
  Only their average obeys the positive identity
  `mean(v_i)=n^-1 sum_ij Q_ij^4 >= 1`.

## 2. A stronger, whole-functional Rademacher comparison

The per-pair replacement in the audited artifact is valid. The following
alternative keeps the matrix weights intact and improves its error rate.

Let `X_i` be the fixed vector of injective marked-tree polynomials, let
`h3(t)=(t^3-3t)/sqrt(6)`, and define

`T(s) = n^-1 F(X(s))^T Q h3(Bs)`,

where `F(X(s))` is the vector with coordinates `F(X_i(s))`. Suppose first
that `F` is a fixed polynomial. All estimates below also hold for a fixed
smooth function whose needed derivatives have polynomial growth.

Consider independent hybrid inputs, each standard Gaussian or Rademacher,
and replace coordinate `k`. Since the actual tree kernels are injective,
every `X_i` is exactly affine in this coordinate. Coefficient counting
and fixed-degree hypercontractivity give, for every fixed finite `p`,

`sup_i ||partial_k X_i||_p <= C_p n^-1/2`.

The same bound holds in a hybrid product space. For example, condition on
which coordinates are Gaussian, apply the Rademacher and Gaussian
hypercontractive inequalities successively to the finite-degree polynomial,
and use that both families have orthonormal multilinear monomials.
Evaluation at a deterministic value `s_k=t` only adds a fixed polynomial
factor in `1+|t|`, with finite moments under both replacement laws.

Let `U(s)=F(X(s))` and `V(s)=h3(Bs)`. For `0 <= r <= 4`, the chain rule,
affineness, Holder, and the preceding derivative bound imply

`|| ||partial_k^r U||_2 ||_p <= C_p n^((1-r)/2)`.

The same bound holds for `V`; for `r=4` its derivative is zero. To check
the vector norm without assuming coordinate independence, use

`|| ||W||_2 ||_p <= (sum_i ||W_i||_p^2)^(1/2)`, for `p >= 2`.

Consequently, using only `||Q||op <= L^2`,

`E|partial_k^3 T|`

`<= n^-1 L^2 sum_(r=0)^3 binom(3,r)`

`   * || ||partial_k^r U||_2 ||_2`

`   * || ||partial_k^(3-r) V||_2 ||_2`

`<= C_L n^-3/2`.

Every product on the right has order `n^-1/2` before the factor `1/n`.
Taylor expansion through order two therefore gives total comparison error
`O_L(n^-1/2)` after replacing all `n` inputs. Uniform polynomial-in-`t`
versions of the displayed estimate justify the integral Taylor remainder.

If four derivatives are available, the same calculation gives

`E|partial_k^4 T| <= C_L n^-2`.

Gaussian and Rademacher variables also match their third moments, so
Taylor expansion through order three improves the total comparison to
`O_L(n^-1)`. Thus, for fixed polynomial `F`,

`E T(S) = E T(N) + O_L(n^-1)`.

This proves the needed Rademacher transfer directly at the weighted
functional level. It does not require any entrywise bound on `Q`, any
coherence condition, or joint convergence over all pairs of roots.

## 3. Scope of the resulting identity

Combining the independently checked Gaussian calculation with Section 2
gives, uniformly for `||B||op <= L`,

`n^-1 E F(X)^T Q h3(G)`

`= E[F(Z) h3(Z_edge)] * n^-1 sum_ij Q_ij^4 + o(1)`.

Here `Z` denotes the independent Gaussian limit of the old rooted fields,
not the transported cubic. Polynomial approximation and the operator-norm
Cauchy--Schwarz bound extend this identity to fixed bounded functions that
are continuous outside a Gaussian-null set, as in the audited artifact.

The unresolved next step is a local response/endpoint identity that
controls expressions such as

`n^-1 E [a(X) circ B h3(G)]^T B F(X)`.

The factor `a(X_i)` cannot be moved through `B`, and independence of the
single-root limit of `B h3(G)` from `X_i` alone does not prove this weighted
identity. No such inference is included in this audit.

## 4. Audit of the all-Hermite extension

The subsequently added Sections 7--8 of the algebra artifact also pass.
For every fixed `r >= 2`, the same argument gives

`n^-1 E F(X)^T Q h_r(G)`

`= E[F(Z)h_r(Z_edge)] * n^-1 sum_ij Q_ij^(r+1) + o(1)`.

The only new exceptional class is a tree of degree exactly `r`. Such a
degree must be odd. The collision-allowed contraction recursion

`phi_T^j = B diag(b_j) product_children(phi_child^j)`

is correct: the top marked vertex contributes the `b_j` factor, the root
edge the `B` multiplication, and disjoint child free-vertex sets give the
coordinatewise product. Its constants depend only on the fixed tree.
Every non-edge column has norm `O(n^-1/2)`. A non-star has a non-edge
child and consequently improves to `O(n^-1)` column norm; its whole
matrix therefore has Frobenius norm `O(n^-1/2)`. A star is handled by the
explicit matrix `B[B circ Q^(circ(r-1))]`, whose middle factor has
operator norm `O(n^-1/2)`. Thus every exceptional weighted diagonal
vanishes. Collision deletion is controlled against the bounded-norm
tensor `sum_j Q_ij b_j^(tensor r)`, not against individual terms.

The whole-functional replacement of Section 2 works for every fixed
polynomial `h_r`, including `r > 3`: its fourth derivative need not vanish,
but has the same `n^-3/2` vector-norm bound, so the result is unchanged.

For an odd response with conditional edge projection `f(g)`, removing its
linear part and normalizing yields a positive residual direction. The
series in Section 8 is genuinely termwise nonnegative because only odd
Hermite degrees occur and hence the exponents `r+1` are even. Its matrix
factors lie in `[1,L^2]`, giving uniform tail control. The stated scope is
important and correct: a positive direction at each fixed smooth slack
baseline does not itself bound the gain uniformly relative to a vanishing
smoothing loss, a vanishing slack, or a spectral-deletion loss.
