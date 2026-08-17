# Adversarial audit: odd-Walsh rooted amalgamation

**Verdict: PASS.**  The reconstruction formulas, rooted orbit conclusion,
associativity statement, three lower-bound families, and the one semantic
leading-scale separation are correct as stated.  I found no counterexample.
The draft also draws the crucial scope boundary correctly: its general
minimality and bit lower bounds concern an orbit-complete universal carrier,
not the scalar optimum of every fixed graph.  Two verifier-coverage comments
at the end are worthwhile follow-ups, but neither is a mathematical defect in
the theorem.

## 1. Reconstruction and converse minimality

Let the two quotient coefficient spaces be represented by the injective maps

```math
\bar\alpha:U_a\hookrightarrow V,
\qquad
\bar\beta:U_b\hookrightarrow V.
```

For the concatenated coefficient map

```math
\gamma(c,d)=\alpha(c)+\beta(d),
```

direct expansion gives

```math
B(\gamma(c,d),\gamma(c',d'))
=G_a(c,c')+G_b(d,d')
 +\kappa([c],[d'])+\kappa([c'],[d]).
```

There is no hidden choice of quotient representatives here: if an argument
is changed by an internal relation, its image in `V` is unchanged.  Likewise,

```math
\ker\gamma
=\{(c,d):([c],[d])\in J_{ab}\},
\qquad
\gamma^{-1}(\omega)
=\{(c,d):([c],[d])\in Z_{ab}^{\times}\}.
```

This proves (RA.12)--(RA.14), including tuples with internal dependencies.
The converse is also exact.  With the two coordinate blocks marked, the
off-diagonal block of the combined Gram form descends to `kappa`; the images
of the combined relation kernel and root fibre in `U_a\oplus U_b` are exactly
`J_ab` and `Z_ab^times`.  Hence, **given the isolated states**, the proposed
datum and the combined rooted orbit state are mutually recoverable.

The word "minimal" in RA.1 must therefore be read in its stated information
sense: any lossless summary required to recover the combined rooted orbit
must encode a value from which all three objects can be recovered.  It does
not assert that the displayed tuple is the shortest bit-level serialization,
nor that every scalar graph optimum separates all its values.  The draft
already says both things explicitly.

The appeal to the preceding odd-dimensional orbit theorem is valid.  In odd
dimension, `omega` is anisotropic and

```math
V=\langle\omega\rangle\perp\omega^\perp,
```

with `omega^perp` nondegenerate symplectic.  Equality of Gram form, relation
kernel, and root fibre induces an isometry of label spans that respects
`omega`; symplectic Witt extension then gives a global orthogonal map.
Simultaneously applying that map to the two Walsh coordinates is a coordinate
permutation, fixes the common Walsh kernel, and carries each labelled child
to its counterpart.  Thus the conclusion is about the full Boolean
landscape, not merely its coefficient matrix.

## 2. Associativity

The displayed quotient

```math
(U_a\oplus U_b)/J_{ab}
\simeq \operatorname{im}(\bar\alpha+\bar\beta)
```

is the first isomorphism theorem.  For three pieces, either parenthesization
is the presentation induced by the same map

```math
U_a\oplus U_b\oplus U_c\longrightarrow V,
\qquad (u,v,w)\longmapsto
\bar\alpha u+\bar\beta v+\bar\chi w.
```

Its pullback form, kernel, and inverse image of `omega` are independent of
the order of quotienting.  This proves RA.2.

There is an important and correctly stated limitation: this is associative
composition of **presented actual spans with newly supplied relative data**.
It is not a binary operation on the isolated orbit states alone.  RA.6 shows
why such a stronger assertion would be false.

## 3. Cross-form family RA.3

The dimension hypothesis `m >= 2(r+s)+1` supplies `r+s` symplectic pairs in
`H=omega^perp`.  In

```math
b_j^K=p_{r+j}+\sum_i K_{ij}q_i,
```

the private `p_(r+j)` coordinates imply all of the following:

- the `b_j^K` are independent for every `K`;
- their span is disjoint from `span(p_1,...,p_r)`;
- their internal Gram matrix is zero, independently of `K`.

Both spans lie in `H`, while `omega` does not, so all isolated and combined
root fibres are empty.  Finally

```math
B(p_i,b_j^K)=K_{ij}.
```

Thus the isolated rooted states, `J=0`, and root datum really are fixed while
all `2^(rs)` cross forms occur.  The `rs`-bit lower bound is exact for any
lossless orbit-complete encoding on this family.  There is no accidentally
varying isolated state.

I additionally checked this exhaustively for

```text
(r,s)=(1,1),(1,2),(2,1),(2,2),(2,3),(3,2).
```

## 4. Intersection family RA.4 and its bit count

For an invertible matrix `P`, both ordered tuples are bases of the same
totally isotropic `r`-space.  Consequently their isolated Gram matrices and
relation kernels are identical, the cross form vanishes, and no root fibre
occurs.  The mixed equality equation is exactly

```math
u=Pv,
```

so distinct `P` give distinct relation kernels.  Exhaustive independent
checks give respectively `1`, `6`, and `168` distinct graphs for
`r=1,2,3`, agreeing with `|GL(r,2)|`.

The asymptotic count is also correct:

```math
\log_2|GL(r,2)|
=r^2+\sum_{j=1}^r\log_2(1-2^{-j}).
```

The final sum converges to a finite negative constant, so the count is
`r^2+O(1)`, not merely `r^2-O(r)`.

## 5. Semantic visibility of a coincidence bit

For odd `m >= 5`, the singleton labels

```math
a=(1,1,1,0,\ldots,0),
\qquad b=(0,0,1,0,\ldots,0)
```

have the same singleton rooted state, and both pair to `1` with the
one-dimensional span represented by the repeated endpoint tuple `(a,a)`.
Neither combined span contains `omega`.  The only differing amalgamation
resource is therefore the coincidence correspondence: the middle `a` lies
in the endpoint span, while `b` does not.

The imported calculation from LG.1 is normalized consistently.  The word
`(a,b,a)` has a Boolean witness saturating three child terms of size
`n^(3/2)/2` and two bridge terms of size `n^(3/2)`, hence maximum
`7 n^(3/2)/2`.  For `(a,a,a)`, child and bridge involutions anticommute, and
the normalized block operator has norm `sqrt(3)`.  Its Boolean vector has
squared norm `3n`, giving the upper bound

```math
\frac{3\sqrt3}{2}n^{3/2}.
```

The stated gap `(7-3 sqrt(3))n^(3/2)/2` follows.  This proves semantic
necessity for this one intersection bit at leading scale.  It does **not**
prove that all `Theta(r^2)` intersection bits or all cross/root bits are
visible to one scalar maximum, and the draft does not claim that.

## 6. Root-fibre independence RA.5

Let `a=e_1+e_2`.  Since `m` is odd, both

```math
b^+=\omega+a,
\qquad b^-=e_3
```

are nonzero odd vectors.  Their singleton rooted states agree.  Both are
orthogonal to `a`, and both pairs are independent, so `kappa` and `J` agree.
But `a+b^+=omega`, whereas the four vectors in `span(a,b^-)` do not include
`omega` for `m>=5`.  Thus (RA.25) is correct and supplies two otherwise
identical unrooted amalgams.  I repeated this check for every ordered choice
of three distinct coordinates at `m=5,7,9`; no accidental isolated-state
variation occurs.

This is an orbit-complete root-bit lower bound only.  No scalar-response
separation is claimed for RA.5.

## 7. Pairwise-invisible ternary relation RA.6

The vectors `p_1,p_2,p_3` lie in a common totally isotropic subspace of
`H`.  In both triples, every singleton is a nonzero isotropic vector.  Every
two-position tuple consists of two distinct independent isotropic vectors,
has zero cross form and coincidence correspondence, and has empty root
fibre.  Nevertheless

```math
p_1+p_2+(p_1+p_2)=0,
```

while `(p_1,p_2,p_3)` is independent.  The global relation kernels differ.
The obstruction is therefore genuinely ternary and is not hidden in a
varying singleton or pair state.  Independent checks over all ordered choices
of three symplectic `p`-vectors in dimensions `m=7,9` reproduce it.

This example also validates the conceptual conclusion: pairwise edge labels
cannot be reused indefinitely.  The accumulated presented span is the
additional dynamic memory that makes later coincidence tests meaningful.

## 8. Complexity and scope

The `O(t^2)` storage claim is correct: a symmetric binary Gram form, a
row-reduced basis of a subspace of `F_2^t`, and at most one root-coset
representative all fit in that size.  Raw labels take `mt` bits.  Hence the
orbit carrier is `o(mt)` when `t=o(m)`.  This is a strict quotient of the
linear-label presentation and is vastly smaller than the full Boolean
energy landscape.

It is not, and is not advertised as, an efficient evaluator of the scalar
maximum.  Reconstructing an orbit representative can leave a difficult
Boolean optimization.  Likewise, the lower bounds count distinct rooted
orbits, so they apply to a universal orbit-complete carrier.  They cannot be
quoted as response-metric packing bounds without additional exposed-query
arguments.  Only the explicit path example currently supplies such an
argument.

## 9. Verifier audit

Running

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_odd_walsh_amalgamation.py
```

returns

```text
odd-Walsh rooted amalgamation checks passed: 6027
```

and the dependency verifier returns

```text
linear-label Walsh Gram/rooted-orbit checks passed: 181
```

The first program correctly checks exhaustive reconstruction at `m=3`,
random reconstruction at `m=5,7`, all `2^4` cross forms in RA.3, all six
members of `GL(2,2)` in RA.4, the coincidence-bit invariants, the root-bit
example, the ternary obstruction, and direct equality under both
parenthesizations.

Two nonblocking coverage caveats should be recorded:

1. Its associativity loop retains the raw concatenated labels and calls the
   two-piece reconstruction routine again.  It therefore checks the formulas
   under both parenthesizations, but does not implement a serializer for an
   abstract quotient presentation and then compose that serialization.  The
   first-isomorphism proof above supplies the missing abstract argument.
2. It imports rather than recomputes the large Boolean path gap; the separate
   LG verifier checks the bent witness and operator identity.  Running both
   verifiers is therefore the appropriate full audit.

Neither caveat weakens the rigorous theorem.  A future engineering revision
could strengthen the verifier by representing quotient bases abstractly and
by calling the LG structural checks, but no repair to the mathematical draft
is required.

