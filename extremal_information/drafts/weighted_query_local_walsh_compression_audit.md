# Independent audit: weighted query-local Walsh compression

**Verdict: PASS AFTER SCOPE AND DECODER REPAIRS.**  The interaction-mass
estimate, simultaneous local-state upper bound, `O(sum_C |C|^2)`
serialization, and path/dense regime separation are correct.  Signed bridge
weights and overlapping query partitions cause no failure.  The result is a
genuine combination of the semantic ambient-Witt quotient with a uniform
interaction truncation, although the truncation inequality itself is an
elementary Lipschitz argument.

Four repairs should be made before canonical promotion:

1. state explicitly that WQ.1 concerns the **upper maximum** and define the
   onsite normalization `H_a(x)=x^TC_ax/2`; if the absolute maximum is also
   claimed, use the two-sided component decoder in Section 3 below rather
   than the sum of component optima;
2. allow arbitrary real onsite coefficients (as the proof and verifier
   already do), and define active vertices so a zero-child endpoint of a
   nonzero bridge is not discarded;
3. require the partitions in the simultaneous theorem to be fixed public
   query data, and explain why a maximal local state restricts to every
   smaller part;
4. qualify rooted queries as componentwise after truncation, unless any
   shared auxiliary variable is added to the interaction support.

The dense estimate can also be sharpened by deleting the unnecessary `-t`
term.

## 1. Normalization and the one-query inequality

Let `E=F_2^(2m)`, `n=|E|=2^(2m)=q^2`, and let `W` be its unnormalized Walsh
matrix.  Then

```math
W^2=nI,
\qquad \lVert W\rVert_{2\to2}=\sqrt n.                 \tag{AWQ.1}
```

For Boolean vectors, `||x||_2=||y||_2=sqrt(n)`, so

```math
|x^TWy|\le n^{3/2}=q^3.                                \tag{AWQ.2}
```

This is exactly the normalization used in Theorems 21.9 and 21.18.  Equality
is possible because even-dimensional Walsh space has Boolean bent vectors.
The verifier uses `n=4`, where the numerical ceiling is `8`.

The natural fully weighted statement is

```math
E_{G,a}(x)=\sum_v h_v H_{a_v}(x_v)
            +\sum_{uv}w_{uv}x_u^TWx_v,
\qquad H_a(x)={1\over2}x^TC_ax,                        \tag{AWQ.3}
```

with arbitrary real `h_v,w_uv`.  If `D_P(x)` is the sum of the deleted
cross-part edge terms, then, pointwise,

```math
|D_P(x)|
 \le \sum_{uv:\,[u]_P\ne[v]_P}|w_{uv}|n^{3/2}
 =d_G(P)n^{3/2}.                                       \tag{AWQ.4}
```

The elementary inequality

```math
|\max f-\max g|\le\sup_x|f(x)-g(x)|                   \tag{AWQ.5}
```

proves WQ.6.  No positivity is used.  Thus negative bridge weights, negative
onsite weights, and zero onsite weights are harmless.  In particular, a
vertex with `h_v=0` but a nonzero incident bridge remains active; its label
may be irrelevant to that one term, but retaining it only makes the proposed
upper carrier redundant, not incorrect.

The same pointwise estimate also gives

```math
\left|\max_x|E_{G,a}(x)|-\max_x|E_{G[P],a}(x)|\right|
 \le d_G(P)n^{3/2}.                                    \tag{AWQ.6}
```

So the perturbation bound itself extends to the absolute objective.  The
decoder formula needs a separate repair, described below.

## 2. Why the local semantic carrier is sufficient

For a part `C`, Theorem 21.18 says that

```math
I_a(C)=(G(a|_C),R(a|_C))                               \tag{AWQ.7}
```

determines the complete unrooted Boolean landscape for arbitrary real
onsite and edge weights supported on `C`, up to a common coordinate
permutation.  This is stronger than equality of one optimum: it preserves
the maximum, minimum, absolute maximum, histogram, and multiplicities.

After cross-part edges are deleted, variables belonging to distinct parts
are disjoint.  For the upper maximum,

```math
\max_xE_{G[P],a}(x)
 =\sum_{C\in P}\max_{x_C}E_C(x_C).                    \tag{AWQ.8}
```

Consequently the list `(I_a(C))_(C in P)` determines the decoded value.
This is precisely where the ambient-Witt semantic collapse is used; without
it, one would have to store the full component landscapes or the raw
`m|C|` label bits.  WQ.4--WQ.6 then remove weak compatibility between those
exact local quotients.  In that limited but mathematically meaningful sense,
the theorem genuinely combines semantic orbit compression with
interaction-mass truncation.

The claim is not an information lower bound and should not be presented as
one.  Nor does the deletion lemma alone depend on Walsh algebra.  Its new
content in this setting is the sub-ambient exact carrier for every retained
part and the resulting query-local scaling laws.

## 3. Upper, absolute, and pairwise response errors

The phrase "returns the sum of the exact component optima" is correct only
for the upper maximum.  For an absolute maximum define, for each part,

```math
M_C^+=\max E_C,
\qquad M_C^-=\min E_C.                                \tag{AWQ.9}
```

Then the exact truncated decoder is

```math
D_P^{\rm abs}
=\max\left\{\sum_CM_C^+,-\sum_CM_C^-\right\}.        \tag{AWQ.10}
```

It is generally unsafe to sum the component absolute maxima: different
components can attain their larger absolute values with incompatible signs.
The state (AWQ.7) preserves the whole component landscape, so it supplies
both quantities in (AWQ.9), and (AWQ.6) proves the same one-copy error
`d_G(P)n^(3/2)`.

The factor accounting in the draft is otherwise correct.  If two systems
have the same stored carrier and are decoded to the same `D`, then

```math
|M(a)-M(b)|
 \le |M(a)-D|+|M(b)-D|
 \le2d_G(P)n^{3/2}.                                    \tag{AWQ.11}
```

There is no factor two in the error of either system relative to its own
decoded value.  This applies to upper or absolute maxima after choosing the
corresponding decoder.

## 4. Simultaneous and overlapping partitions

The simultaneous claim is valid provided every `P_theta` is selected from
the declared graph query and error budget, not adaptively from hidden label
data.  Equivalently, the query family and its partitions are public decoder
data.  If adaptive partitions were permitted, their description would have
to be included in the code length.

Every part lies inside an inclusion-maximal part occurring in the finite
family.  If `C subseteq C'`, then

```math
G(a|_C)=G(a|_(C'))|_(C\times C),
\qquad
R(a|_C)=R(a|_(C'))\cap F_2^C.                         \tag{AWQ.12}
```

Thus `I_a(C')` computes `I_a(C)`.  Overlap between two incomparable maximal
parts is also harmless.  Their stored charts may duplicate data, but each
scalar query is decoded independently and no common coordinate conjugacy
between two disconnected parts is needed.  This is the same query-local
principle as Theorem 21.19.

For a maximal part of size `k`, its binary Gram matrix and a row-reduced
basis of its relation kernel use `O(k^2)` bits.  Therefore

```math
O\left(\sum_{C\ \mathrm{maximal}}|C|^2\right)         \tag{AWQ.13}
```

is a valid direct upper bound even for overlapping parts.  It is not a
minimality statement; overlap consistency or semantic coincidences may make
a smaller serialization possible.

For rooted states, an empty flag and one affine root-fibre representative
add only `O(k)` bits per part.  Independent rooted fields fixed by each
component's equivariant coordinate relabelling are covered.  A shared
external Boolean variable or any other auxiliary degree of freedom joining
two nominal parts makes them one support and is not covered by separate
component decoding.

## 5. Path and dense regimes

The path estimate is correct.  Put

```math
s=\lceil1/\eta\rceil+1.
```

Cut after each consecutive run of `s` vertices.  The number of deleted unit
edges is

```math
\left\lceil{t\over s}\right\rceil-1
 <{t\over s}<\eta t,                                  \tag{AWQ.14}
```

while

```math
\sum_C|C|^2\le s\sum_C|C|=st=O(t/\eta).               \tag{AWQ.15}
```

Hence the additive upper- or absolute-response error is at most
`eta t n^(3/2)`.  This is a linear-in-block-count normalized error statement,
not an `O(1)` total-error theorem.

For `K_t`, the cross-part edge count is exactly

```math
{1\over2}\left(t^2-\sum_C|C|^2\right).                \tag{AWQ.16}
```

There is no diagonal correction.  If every edge has `|w_e|>=c` and
`d_G(P)<=eta t^2`, then

```math
\sum_C|C|^2
 \ge \left(1-{2\eta\over c}\right)t^2.               \tag{AWQ.17}
```

This sharpens WQ.10 by removing `-t`.  For fixed `eta<c/2`, the particular
component-deletion carrier is therefore quadratic.  As the draft correctly
warns, (AWQ.17) is not a lower bound for arbitrary semantic response codes;
signed edge cancellations may permit a different architecture even though
the total-variation deletion budget is large.

## 6. Adversarial cases

I found no scope-breaking example in the requested cases.

- **Signed weights.**  They are absorbed by absolute values in (AWQ.4), and
  exact local conjugacy preserves their signs.
- **Zero child weights.**  They do not affect the perturbation argument.  A
  bridge-active zero-child vertex must remain in the support; a completely
  absent vertex can be discarded.
- **Absolute maxima.**  The same error bound holds, but the corrected
  decoder is (AWQ.10), not a sum of local absolute optima.
- **Overlapping query partitions.**  The restriction identity (AWQ.12)
  validates the union-of-maximal-parts carrier.  Public fixed partitions are
  essential for the stated description count.

The most plausible hidden failure would have been a need to synchronize the
ambient-Witt conjugacies across parts.  No such synchronization is needed
for the declared scalar queries because deletion leaves disjoint Boolean
variables.  It would reappear immediately for a shared rooted variable or a
later continuation that reconnects the parts; that continuation must be
included as a new query and its connected support stored accordingly.

## 7. Verifier coverage

I ran

```text
./.venv/bin/python \
  extremal_information/experiments/verify_weighted_query_local_walsh_compression.py
```

and obtained

```text
weighted query-local Walsh checks passed: 1874
```

The exhaustive `m=1`, four-block loop checks WQ.6 for three signed integer
graph families, all sixteen linear-label words, and all fifteen canonical
partitions.  Its onsite coefficients include negative values, but not zero.
The final loop checks the elementary complete-graph cross-edge identity
through `t=7`.

The script does **not** verify:

- ambient-Witt equivalence for two different tuples with equal `(G,R)`;
- recovery of smaller states from overlapping maximal parts;
- the absolute decoder (AWQ.10) or the factor-two pairwise response bound;
- a zero-child, bridge-active endpoint;
- the path construction and its `O(t/eta)` serialization;
- the weighted dense lower bound, rooted extension, or arbitrary real
  weights.

Those omissions do not undermine the proof, because the semantic step is
already Theorem 21.18 and the remaining estimates are explicit.  The script
should nevertheless be described as an exact regression test for the
pointwise deletion bound, not as an exhaustive verifier of WQ.1--WQ.2.  A
stronger regression test should add zero onsite weights, compute both minima
and maxima for (AWQ.10), use two overlapping partitions, and check
(AWQ.14)--(AWQ.17) arithmetically.

## 8. Recommended canonical disposition

Promote the result as a short theorem/corollary immediately after the exact
query-local Walsh theorem, after making the four repairs listed at the top.
The canonical statement should distinguish:

1. a general interaction-deletion Lipschitz lemma;
2. its Walsh consequence, where every retained part has the exact semantic
   state `(G,R)` (or `(G,R,Z)` for a genuinely rooted component);
3. the path upper bound and the dense **architecture-specific** obstruction.

This is theorem-level progress: it turns the exact component quotient into
an approximate response carrier whose complexity is controlled by
interaction mass and yields a nontrivial `O(t/eta)` path law.  Its correct
interpretation is an upper theorem for weighted query-local response
complexity, not semantic minimality of the carrier and not an information
lower bound for dense Walsh queries.
