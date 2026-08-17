# Scale-sensitive spectral bridge compression

Status: rigorous task-local theorem, independently audited after the scope
repairs below.  This note
connects the exact finite-rank roof algebra to sparse and synchronized
full-rank bridges at the `n^(3/2)` scale.  It does not assert compression for
a generic dense sign bridge.

## 1. Uniform spectral-tail replacement

Let `X subseteq {+-1}^p`, `Y subseteq {+-1}^q`, let `H:X->R` and `K:Y->R`
be arbitrary internal landscapes, and let

```math
\mathcal V_R(H,K)
=\max_{x\in X,y\in Y}\{H(x)+K(y)+x^TRy\}.               \tag{ST.1}
```

### Theorem ST.1 (numerical rank is the scale-sensitive interface)

For every two real `p by q` bridges `R,S`,

```math
|\mathcal V_R(H,K)-\mathcal V_S(H,K)|
\le\sqrt{pq}\,||R-S||_{2\to2}.                           \tag{ST.2}
```

The **replacement bound** is uniform over the two internal landscapes and
over every later common max-type future on `(x,y)`, because it is pointwise.
This arbitrary-future statement does not apply to the compressed roof state
introduced in Section 2. If `S=R_r` is a best rank-`r` spectral truncation,
then

```math
|\mathcal V_R-\mathcal V_{R_r}|
\le\sqrt{pq}\,\sigma_{r+1}(R).                           \tag{ST.3}
```

The constant in (ST.2) is sharp: for Boolean `x_0,y_0`, the rank-one matrix

```math
R-S=\sigma{x_0\over\sqrt p}{y_0^T\over\sqrt q}          \tag{ST.4}
```

attains equality when the two landscapes pin `(x_0,y_0)` (or their state
sets are singletons).

#### Proof

For every Boolean pair,

```math
|x^T(R-S)y|
\le||x||_2||R-S||_{2\to2}||y||_2
=\sqrt{pq}\,||R-S||_{2\to2}.                            \tag{ST.5}
```

A pointwise `delta`-bound changes a maximum, even after adding an arbitrary
shared future, by at most `delta`.  Eckart--Young gives (ST.3), and direct
substitution proves sharpness. `square`

At balanced size `p=q=n` and target error `epsilon n^(3/2)`, only singular
directions above `epsilon sqrt(n)` can matter uniformly.  Define the target-
scale numerical rank

```math
r_\epsilon(R)=#{j:\sigma_j(R)>\epsilon\sqrt n}.         \tag{ST.6}
```

Then every bridge is response-equivalent, up to `epsilon n^(3/2)`, to a
rank-`r_epsilon(R)` bridge.

## 2. The induced roof state

Factor a rank-`r` truncation as

```math
R_r=U\Sigma V^T,
\qquad
\phi(x)=\Sigma^{1/2}U^Tx,
\quad
\psi(y)=\Sigma^{1/2}V^Ty.                               \tag{ST.7}
```

The exact interaction is `<phi(x),psi(y)>`.  Hence Theorem 18.2 applies to
the two upper roofs over these `r` features; no optimizer or full
**state-specific** response landscape is stored. The upper roof is precisely
the full response to the declared retained-feature fields. Treat the bridge
and its factorization as shared query data.

Now restrict continuations to depend on a child through its retained feature
and to have total dual field radius at most `P`. This includes one opposite
child of feature radius `P`. A fixed number `d` of such ports has radius
`O_d(P)`, with all constants below depending on the declared port bound. It does
not include an arbitrary future that can distinguish states in one feature
fibre. If

```math
||R||_{2\to2}\le C\sqrt n,                              \tag{ST.8}
```

then both feature radii are at most

```math
P=\sqrt C\,n^{3/4}.                                     \tag{ST.9}
```

For `0<epsilon<=1`, quantize each feature ball at radius

```math
\eta={\epsilon n^{3/4}\over4(C+1)}                      \tag{ST.9a}
```

and retain roof heights at mesh `epsilon n^(3/2)/4`.  Since

```math
2P\eta+\eta^2\le {5\epsilon\over16}n^{3/2},             \tag{ST.9b}
```

the roof-quantization error is below
`epsilon n^(3/2)`.  This is additional to the spectral-tail error in (ST.3),
so a prescribed total budget must be divided between them. Buckets more than
the maximum possible field oscillation
below the top height may be clamped, so only `O_C(1/epsilon)` height levels
are response-relevant.  The number of feature buckets is

```math
\left(1+O_C(1/\epsilon)\right)^r.                        \tag{ST.10}
```

Thus the complete table description has

```math
\exp(O_C(r\log(1/\epsilon)))                            \tag{ST.11}
```

quantized entries, plus one unrestricted additive baseline such as `max H`.
This is constant for fixed numerical rank, polynomial
when `r=O(log n)`, and subexponential in `n` whenever `r=o(n)`.  The
worst-case `2^(Omega(r))`-bit response family of Theorem 18.3 shows that an
exponential dependence on `r` cannot be removed for general finite-feature
roofs. That lower bound has not been embedded into every Boolean SVD port.

The statement is deliberately about an explicit table realization.  A
particular structured child may have a much smaller symbolic roof.

## 3. Graph composition

### Corollary ST.2 (local factor presentation for spectrally compressible ports)

Put arbitrary landscapes on `k` Boolean blocks of equal size `n`, and put a
bridge `R_e` on every edge of a graph `G`.  Choose rank-`r_e` matrices `S_e`
with

```math
||R_e-S_e||_{2\to2}\le\delta_e\sqrt n.                  \tag{ST.12}
```

Then replacing every bridge simultaneously changes every global maximum and,
pointwise, every later shared max-type response by at most

```math
n^{3/2}\sum_{e\in E(G)}\delta_e.                         \tag{ST.13}
```

At vertex `v`, concatenate the singular features of its incident truncated
bridges.  This gives an exact roof algebra for the truncated graph of local
dimension at most

```math
d_v\le\sum_{e\ni v}r_e.                                 \tag{ST.14}
```

For bounded degree, uniformly bounded operator scale, and a declared fixed
port/field budget, its quantized **local factor** table has
`exp(O(d_v log(1/epsilon)))` entries. Gluing while retaining every exposed
port is associative because all cross terms are bilinear in those features;
errors are charged once per omitted physical edge, not once per
parenthesization.

This is not a bounded global quotient for arbitrary bounded-degree graphs.
Eliminating a region creates a joint factor on its boundary ports, whose
dimension is governed by the cut or treewidth and can be linear even when
all vertex degrees are bounded. For `O(k)` physical edges, a total error
`epsilon k n^(3/2)` is obtained by assigning an `O(epsilon)` replacement and
quantization budget per edge/block.

#### Proof

Apply (ST.2) pointwise on every edge and sum.  Factoring every `S_e` and
concatenating its endpoint features makes the retained Hamiltonian a sum of
internal roofs and explicit bilinear pairings.  The exact roof algebra and
the quantization bound then apply vertexwise. `square`

## 4. What the theorem predicts

The result unifies several previously separate bridge observations at the
declared scale.

1. A bridge `alpha I+beta J` with `alpha=O(1)` has only one singular
   direction above `epsilon sqrt(n)` for large `n`; the identity residual has
   total response at most `O(n)=o(n^(3/2))`.  Thus its leading-scale state is
   rank one even without invoking the stronger exact rearrangement theorem.
   When `beta` is constant the leading singular value is order `n`, so the
   bounded-operator table estimate (ST.8) does not apply even though the
   numerical-rank conclusion does.
2. Every bipartite bridge with uniformly bounded row and column degrees has
   operator norm `O(1)` and is uniformly
   negligible at scale `n^(3/2)`, although it may be information-heavy at
   its natural scale `n`.  This resolves the apparent tension with the
   matching-bridge lower bound.
3. For every threshold strictly inside the iid sign-matrix spectral bulk, an
   iid dense sign bridge typically has order `n` singular directions at
   scale `Theta(sqrt(n))`. Spectral truncation therefore predicts an
   extensive port rather than a compression.  The response-packing theorems
   show that this negative prediction is semantically real for broad child
   classes.

Consequently the transition is not literal algebraic rank.  It is the count
of singular directions visible at the error scale, supplemented by any
nonlinear synchronization that can collapse the retained roof.  The theorem
does not compress an arbitrary dense sign bridge and does not use information
specific to near-minimizing sign quadratics.
