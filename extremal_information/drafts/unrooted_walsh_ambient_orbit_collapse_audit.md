# Independent audit: unrooted Walsh ambient-orbit collapse

**Verdict: PASS WITH A DECLARATION REPAIR.**  The ambient `2m`-dimensional
Witt-extension argument is correct and proves considerably more than equality
of scalar maxima: it gives one Boolean-coordinate permutation conjugating the
whole family of child and bridge matrices.  Arbitrary real scalar vertex and
edge weights are therefore covered.  The only repair needed before promotion
is to state the continuation syntax with an explicit port map (and, if desired,
explicit coefficient vectors for synchronously derived linear-combination
labels).  A new ambient-coordinate label or a fixed external Boolean field is
not covered.

## 1. The matrix embedding is exact

Let `V=F_2^m`, `E=V direct-sum V`, and let `W_E` be the order-`2^(2m)`
Walsh matrix.  With the coordinate order `(u,v)`,

```math
W_E=R\otimes R.
```

For `c_a=(0,a)` one has

```math
D_{c_a}=I\otimes D_a,
\qquad
D_{c_a}W_ED_{c_a}=R\otimes(D_aRD_a)=C_a.              \tag{UA.1}
```

Thus the draft has not changed the child family or its normalization.  It
has merely recognized the label modulation as a character of the full
Walsh coordinate space.  The trace of every displayed child is zero, so
retaining or deleting the diagonal makes no difference on the Boolean cube.

## 2. The even-dimensional extension step is valid

Suppose two ordered label tuples have the same Gram matrix and relation
kernel.  The induced map

```math
(0,a_i)\longmapsto(0,b_i)                              \tag{UA.2}
```

is well defined, bijective on the two presented spans, and preserves the
standard bilinear form on `E`.

The characteristic vector of `E` is

```math
\Omega_E=(\omega,\omega).
```

It is absent from both spans, including when the tuples contain zero labels
or are linearly dependent, because every vector in either span has first
coordinate zero whereas `omega` is nonzero.  Adjoining `Omega_E` and fixing
it preserves all pairings:

```math
\Omega_E\mathbin\cdot(0,a)=a\mathbin\cdot a.
```

The right side is preserved by the partial isometry.  Since `dim(E)=2m` is
even, this is exactly the delicate characteristic-two case of the rooted
Witt-extension lemma.  The explicit affine-symplectic proof in the audited
all-parity orbit theorem applies: a partial isometry which fixes the
characteristic vector extends to `O in O(E)`.  There is no illicit use of
the earlier Walsh corollary at a changed matrix size.

I reran

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_unrooted_walsh_ambient_collapse.py
```

which reports `103` exact checks and exhausts all tuples through length three
at `m=2`.  The independent all-parity verifier also reports `4763` checks.
These computations are diagnostics; the extension proof is uniform in `m`
and tuple length.

## 3. Coordinate-conjugacy orientation and weights

Choose `O` with `O(0,a_i)=(0,b_i)`, and let the coordinate permutation be

```math
(Pf)(z)=f(O^{-1}z).
```

Then

```math
PW_EP^T=W_E,
\qquad
PD_{(0,a_i)}P^T=D_{(0,b_i)},
\qquad
PC_{a_i}P^T=C_{b_i}.                                  \tag{UA.3}
```

The apparent transpose ambiguity is harmless but worth checking: orthogonality
gives `O^{-T}=O`, so the modulation character transforms by precisely the
same `O` that sends the embedded labels.  Applying the same `P` in every
block therefore conjugates every matrix

```math
\operatorname{diag}(h_vC_{a_{\sigma(v)}})
+(J_{uv}W_E)_{u\ne v}                                  \tag{UA.4}
```

for arbitrary real `h_v,J_uv`, with no sign, positivity, bipartiteness, or
integrality assumption.  It follows pointwise that upper and absolute
maxima, minima, histograms, and optimizer multiplicities all agree.

## 4. Exact continuation boundary

For promotion, the strongest clean unrooted language is the following.

- A syntax chooses any finite number `t` of variable blocks, a port map
  `sigma:[t]->[k]`, arbitrary real scalar vertex weights, and an arbitrary
  real symmetric bridge-weight matrix.  This includes synchronized repeated
  copies of an exposed label.
- More generally, a syntax may attach to a vertex a fixed coefficient vector
  `c_v in F_2^k` and use the synchronously derived label
  `sum_i(c_v)_i a_i`.  The same proof covers this because the partial
  isometry carries every such derived label to the identically presented
  combination on the other side.
- A syntax may **not** append a vector specified in the old ambient
  coordinates, such as a new fixed label `d in V`, unless its complete
  relative Gram/relation data are supplied as part of the context.  Such a
  label need not be carried correctly by `O`.
- A syntax may **not** add a fixed Boolean pole or linear field.  In
  particular the canonical `omega`-rooted field remembers the original
  `(u,v)` splitting, while the ambient orthogonal map may mix its two halves.

The draft currently writes (UA.4) only with one occurrence of each original
label and then describes repetitions and linear combinations in prose.
Writing the port map, or the coefficient vectors `c_v`, directly in the
definition removes that small scope ambiguity.  Also, `h_v` should be called
a **scalar onsite weight**, not a field; an arbitrary coordinate-dependent
linear field is outside the theorem.

## 5. Consequence for semantic minimality

The root collision is sound.  In odd label dimension, `omega` and `e_1`
have the same singleton Gram/relation state; in even dimension at least four,
the same is true of `omega` and `e_1+e_2`.  Repetition gives the asserted
even-parity relation kernel, while only the constant `omega` tuple has a
nonempty label-space root fibre.  Nevertheless (UA.3) conjugates every
unrooted weighted graph landscape.  The canonical rooted response separates
the singletons by `n^(3/2)/6`.

Therefore the correct query-relative conclusion is exact:

```math
\boxed{
\text{unrooted synchronized Walsh graphs factor through }(G,R),
\quad
\text{whereas the canonical rooted language can require }R_\omega.}
```

This proves sufficiency, not scalar minimality of `(G,R)`.  The finite wind
tunnel's collapse of ten ambient orbit states to three scalar signatures at
`m=2,k=2` is a genuine warning against upgrading the theorem.  The next
lower-bound question is which cross-Gram and relation/intersection bits are
actually exposed by scalar maxima; the present theorem neither assumes nor
answers that question.

