# Independent audit: query-local Walsh amalgamation

**Verdict: PASS AFTER SCOPE REPAIR AND ROOT-FIBRE SHARPENING.**  The local
amalgamation identities, the absence of a required global conjugacy, the
`O(wL)` serialization bound, and the `h`-bit coincidence packing are
mathematically correct.  I found no counterexample to either theorem after
the query class is made precise.  Two presentation repairs are important:

1. an active block with zero child weight but a nonzero incident bridge must
   not be discarded, and a rooted continuation must be genuinely
   componentwise (not share a hidden external variable between components);
2. the newer ambient-orbit result shows that `Z^times` is unnecessary for
   every **purely unrooted** weighted Walsh-graph component.  It is needed
   only on supports carrying a coordinate-rooted/external-field query.

The second point strengthens rather than weakens QL.2: its path queries are
unrooted, and they prove that the local coincidence correspondence `J`
survives the larger unrooted symmetry quotient.

## 1. Declared queries and connected supports

The intended support definition is sound once an active vertex is defined
as a marked block which participates in at least one nonzero term:

```math
h_i\ne0,
\qquad J_{ij}\ne0\text{ for some }j,
\qquad\text{or a declared rooted term uses block }i.       \tag{AQL.1}
```

The connected components are then those of the graph of nonzero bridge
weights on the active vertices; an active vertex with no bridge is a
singleton component.  The phrase "discard zero-weight vertices" is unsafe
as written.  If it means `h_i=0`, it incorrectly deletes an endpoint of a
nonzero bridge.  It should mean "discard blocks absent from every term" or
be replaced by (AQL.1).

For an upper maximum, a disconnected query with energy

```math
E(X)=\sum_D E_D(X_D)
```

and disjoint Boolean variables satisfies

```math
\max_X E(X)=\sum_D\max_{X_D}E_D(X_D).                    \tag{AQL.2}
```

Thus independent coordinate conjugacies may indeed be selected on different
components, and may differ from query to query.  A single global conjugacy is
not part of the response decoder.  If the declared scalar is the absolute
maximum, the literal sum in (AQL.2) must be replaced by

```math
\max\left\{\sum_D\max E_D,-\sum_D\min E_D\right\};       \tag{AQL.3}
```

the conclusion still follows because a local orbit isomorphism preserves
the entire component landscape.  Likewise a rooted field or pole is allowed
provided it is attached independently to one component and transforms
equivariantly there.  A shared external spin/query variable coupling two
nominal components would instead join their supports and is not covered.
The draft should say this explicitly; invariance of an otherwise undefined
"scalar continuation" is not by itself a compositional specification.

## 2. QL.1: exact local reconstruction

Fix one support `C`.  Passing from coefficient spaces to the isolated
quotients gives the local map

```math
\gamma_C(c,d)=\alpha(c)+\beta(d).
```

The three claims in (QL.7) follow directly:

* its pulled-back bilinear form has the two isolated diagonal blocks and
  off-diagonal block `kappa|_(P_C times Q_C)`;
* `gamma_C(c,d)=0` exactly when the quotient images form an element of
  `J cap (P_C direct-sum Q_C)`;
* `gamma_C(c,d)=omega` exactly when they form an element of
  `Z^times cap (P_C direct-sum Q_C)`.

Internal relations cause no ambiguity because they have already been
quotiented out.  Conversely, the off-diagonal Gram block descends to the
restricted cross form, and the images of the combined zero and root fibres
recover the displayed restrictions of `J` and `Z^times`.  Hence the claimed
coarseness is valid in the information order **among summaries required to
reconstruct every local rooted orbit triple**.  It is not semantic
minimality for one scalar maximum, and the draft correctly says so in the
theorem statement.

If `D subset C` is another component support, the data for `D` are obtained
by intersecting/restricting the data for `C` to `P_D direct-sum Q_D`.
Therefore retaining only inclusion-maximal supports loses nothing.  Supports
from different queries need not have compatible local orthogonal maps:
each answer is a scalar query evaluated separately.  Supports in two
components of the same query also need no compatibility because their
Boolean variables are disjoint, subject to the componentwise caveat in
Section 1.

The all-parity rooted relation-form theorem applies here, so QL.1 does not
need an odd-`m` restriction.  This was a possible edge case, but Theorem
21.14 and its exhaustive `m=1,2,3,4` audit cover it.

## 3. Serialization and the `O(wL)` count

For `p=dim P_C`, `q=dim Q_C`, and `d=p+q`, the stated direct encoding is
valid:

* `pq` bits store the cross bilinear form;
* at most `d^2` bits store a row-reduced basis of the subspace
  `J cap(P_C direct-sum Q_C)`;
* the intersection of the affine root fibre with this subspace is either
  empty or a coset of that restricted `J`, so an empty flag and one
  `d`-bit representative suffice.

Thus (QL.8) holds.  Since `p,q<=|C|`,

```math
\sum_C O(|C|^2)\le O\left(w\sum_C|C|\right)=O(wL),       \tag{AQL.4}
```

which proves (QL.9).  This counts the **relative** carrier; the isolated
local states are inputs and are not included in `B_Theta`.

One wording should be narrowed.  Bounded support size alone does not imply a
strict `O(t)` carrier: the family of all two-element supports has bounded
`w` but total incidence `L=Theta(t^2)`.  The linear conclusion is correct for
a partition, bounded-occurrence family, or any declared family with
`L=O(t)`.  In general the theorem gives `O(wL)`, and the union of local charts
can saturate the global quadratic information.

## 4. QL.2: construction and information lower bound

The construction checks in every coordinate.

For

```math
m=3h+r,
\qquad r=1\ (h\text{ even}),\quad r=2\ (h\text{ odd}),
```

`m` is odd and at least one coordinate is unused.  On disjoint
three-coordinate chunks,

```math
a_i=(1,1,1),\qquad c_i=(0,0,1)
```

have norm one and `a_i dot c_i=1`; different chunks are orthogonal.  Hence
every tuple `(b_i^sigma)` has Gram matrix `I_h`, zero relation kernel, and
empty root fibre.  The duplicated endpoint piece has a fixed state, and

```math
a_i\cdot b_j^\sigma=\mathbf1_{i=j}                      \tag{AQL.5}
```

fixes the complete cross form.  The unused coordinate forces `omega` out of
the first, second, and combined spans, so all combined root fibres are
empty.

Chunkwise independence also proves the exact intersection statement:

```math
\operatorname{span}\{b_i^\sigma\}
 \cap\operatorname{span}\{a_i\}
=\operatorname{span}\{a_i:\sigma_i=0\}.                \tag{AQL.6}
```

Thus distinct bit strings give distinct coincidence correspondences, and
the `i`th bit is tested by the three-block path with endpoint word
`(a_i,a_i)` and middle label `b_i^sigma`.

The normalization is also correct.  With `q=2^m` and child block order
`n=q^2`, Theorem 21.12 gives

```math
\max E_(a_i,c_i,a_i)={7\over2}n^{3/2},
\qquad
\max E_(a_i,a_i,a_i)\le {3\sqrt3\over2}n^{3/2}.         \tag{AQL.7}
```

The same coordinate permutation transports the audited bent witness to
every chunk.  Their gap is at least

```math
\Delta n^{3/2},\qquad \Delta={7-3\sqrt3\over2}>0.       \tag{AQL.8}
```

All unused blocks may simply be given zero weights.  Alternatively, if they
are retained as isolated child terms, each has the same isolated maximum and
contributes a bit-independent baseline.  The former convention makes clear
that each active query has only three blocks and that (AQL.8) is a fixed
leading-scale gap for the active query.  If one instead normalizes by the
order of all `3h` nonzero isolated blocks, the normalized fraction acquires
an `h^(-3/2)` factor; the draft's `n` is the order of one Walsh block, not
the order of that enlarged disconnected union.

For two distinct `sigma,tau`, some coordinate query differs by at least
`Delta n^(3/2)`.  If they shared a summary state, a decoder receiving the
query index would return the same answer for both, and the triangle
inequality would force one error to be at least
`Delta n^(3/2)/2`.  Therefore uniform error below this threshold requires
`2^h` states, or `h` bits.  This is a semantic lower bound for the declared
vector of ordinary unrooted scalar maxima, not merely an orbit count.

The result is sharp in the intended interface parameter: there are `t=3h`
marked labels, the maximal supports are `h` disjoint triples, and the exact
relative carrier costs `O(h)`.  Since `m=Theta(h)` and `n=2^(2m)`, this is not
a claim of a linear information rate in the Walsh block order `n`; it is a
linear rate in the number of declared bounded-support ports.

## 5. Verifier coverage

I ran

```text
./.venv/bin/python \
  extremal_information/experiments/verify_query_local_walsh_amalgamation.py
```

and obtained

```text
query-local Walsh amalgamation checks passed: 170
```

The program correctly checks:

* direct local Gram/zero-fibre/root-fibre reconstruction on 36 selected
  tuple/support cases;
* identical isolated second-piece rooted states, fixed cross forms, empty
  combined root fibres, the visible membership bits, and all `2^h` distinct
  coefficient-level coincidence graphs for `1<=h<=6`;
* positivity and arithmetic of `Delta`.

It does **not** itself optimize the three-block Boolean landscapes, construct
the bent witness, verify the anticommuting spectral ceiling, enumerate
orthogonal orbits, or test arbitrary query families.  Those ingredients are
legitimately imported from independently verified Theorems 21.12, 21.14,
and 21.15; the verifier should be described as an exact structural regression
test rather than as a standalone exhaustive proof of QL.1--QL.2.  I also ran
the existing linear-label and rooted-amalgamation verifiers, obtaining `181`
and `6027` passed checks respectively.

## 6. Sharpening from unrooted ambient-orbit collapse

The emerging ambient-orbit theorem changes one conclusion of the draft.
For an unrooted Walsh graph, embed a label as

```math
\iota(a)=(0,a)\in E=\mathbb F_2^m\oplus\mathbb F_2^m.
```

If two finite label tuples have equal Gram matrices and relation kernels,
the induced isometry between their `iota`-spans can be extended in the
ambient even-dimensional dot-product space while fixing its characteristic
vector.  The resulting element of `O(E)` is a coordinate permutation that
fixes the common Walsh matrix and simultaneously conjugates every exposed
child.  Consequently the **entire unrooted weighted graph landscape**, not
only its spectrum, depends on `(G,R)` and not on the label-space root fibre.
I checked the supplied explicit odd/even conjugacies and weighted
three-block identities (`54` passed checks).

This conclusion is query-relative.  It covers graphs on the already marked
tuple (and extensions presented intrinsically as repeats or linear
combinations of that tuple).  An arbitrary new label named in the original
ambient coordinate system need not be respected by the larger conjugacy;
such a label is additional rooted/relative data and must be included in the
support before applying the rule below.

Applied support by support, this yields the sharper mixed query-local rule:

```math
\begin{array}{c|c|c}
\text{component query class}&\text{isolated local state}&
   \text{relative local datum}\\ \hline
\text{unrooted weighted Walsh graph}&(G,R)&(\kappa,J)\\
\text{coordinate-rooted/external-field}&(G,R,Z)&(\kappa,J,Z^\times).
\end{array}                                             \tag{AQL.9}
```

For a mixed declared family, `Z^times` need be stored only on maximal
supports on which some root-sensitive continuation is allowed.  This does
not conflict with the orbit-complete coarseness assertion in QL.1: QL.1
asks to reconstruct a **rooted** orbit triple, whereas the unrooted semantic
query class has a larger symmetry group and a coarser orbit quotient.

QL.2 is unchanged and becomes the sharp part of the unrooted theory.  Its
varying resource is `J`; every root fibre is already fixed and empty.  Thus
it proves that, although `Z^times` vanishes from (AQL.9), the coincidence
correspondence cannot vanish with it.

The last bullet of Section 4 in the unaudited draft is now obsolete: the
"full metric entropy" of the characteristic-root fibre under the declared
unrooted weighted graph family is not open.  It is exactly zero conditional
on `(G,R)`.  What remains open is the semantic minimality of `(G,R)`, or
equivalently which Gram and relation bits unrestricted scalar graph maxima
actually expose.

## 7. Recommended disposition

Promote QL.1 and QL.2 after the following narrow edits:

1. define active vertices by participation in any nonzero term and define
   rooted continuations as independent component terms;
2. retain the rooted theorem as stated, but add the unrooted corollary
   (AQL.9) and delete the obsolete root-entropy claim;
3. qualify the linear-complexity prose by `L=O(t)` rather than bounded `w`
   alone;
4. state explicitly that `n` in QL.2 is one Walsh block order and preferably
   set all unused query blocks to zero weight;
5. describe the verifier's actual structural coverage and cite the imported
   path/orbit verifiers for the semantic step.

With these repairs, the draft gives a genuine strict, query-generated
amalgamation theorem and a matching semantic lower bound for one of its two
surviving unrooted resources.
