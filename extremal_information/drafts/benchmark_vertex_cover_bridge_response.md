# Solution-hidden benchmark: vertex-cover states for Boolean bilinear bridges

Status: derived from contextual future equivalence before consulting the
separator, vertex-cover, or bridge-response literature; exact factorization,
worst-case minimality, and rate--distortion proofs are accompanied by a small
exact verifier.

## 1. Operational problem and convention

Let `L` and `Rgt` be finite index sets and use sign-valued Boolean variables

```math
x\in\{-1,1\}^{L},\qquad y\in\{-1,1\}^{Rgt}.                 \tag{VC.1}
```

Fix a real matrix `R` and an entirely arbitrary internal child landscape
`F:{-1,1}^L -> R`.  The bridge message is

```math
(P_RF)(y)=\max_x\{F(x)+x^TRy\}.                              \tag{VC.2}
```

An arbitrary future landscape `G` on the right returns

```math
\operatorname{Opt}_R(F,G)
=\max_y\{(P_RF)(y)+G(y)\}.                                  \tag{VC.3}
```

Nothing is assumed about the factorization, symmetry, or regularity of `F`.
This is a theorem about the exact outgoing interface once the child landscape
has been supplied, not a promise that an arbitrary oracle table `F` can be
read or optimized in fixed-parameter time.

The `{0,1}` convention is equivalent here.  Substituting sign variables
`2z-1` changes the bilinear matrix only by a nonzero scalar and introduces
unary terms and a constant; the unary terms can be absorbed into the arbitrary
child and future landscapes without changing the support graph.

## 2. Contextual equivalence first

For two child landscapes define their absolute contextual distance by

```math
d_R(F,F')=\sup_G
|\operatorname{Opt}_R(F,G)-\operatorname{Opt}_R(F',G)|.      \tag{VC.4}
```

### Proposition VC.1 (bridge-response isometry)

```math
d_R(F,F')=\|P_RF-P_RF'\|_\infty.                             \tag{VC.5}
```

Modulo a global additive score, the distance is

```math
\bar d_R(F,F')
=\frac12\operatorname{osc}(P_RF-P_RF').                     \tag{VC.6}
```

#### Proof

The maximum in (VC.3) is one-Lipschitz in the sup norm, giving `<=` in
(VC.5).  Conversely, fix a right assignment `y_0` and choose `G(y_0)=0`, with
all other future values sufficiently negative that `y_0` uniquely maximizes
both messages.  The response difference is then the coordinate difference at
`y_0`.  Pinning a coordinate of maximum absolute difference proves equality.
Best approximation of a finite vector by a constant is its midrange, proving
(VC.6). `square`

Thus the canonical exact state for a fixed bridge is its realizable response
`P_RF`; equality of a guessed statistic is relevant only if it forces equality
of this response.  The next section factors that canonical state through a
small support cover.

## 3. The `2^k` cover table

Let `H_R` be the bipartite support graph of `R`: an edge `(i,j)` is present
exactly when `R_ij` is nonzero.  Let a vertex cover split as

```math
C=A\mathbin{\dot\cup}B,qquad
A\subseteq L,\quad B\subseteq Rgt,quad |A|+|B|=k.           \tag{VC.7}
```

Put `U=L\A` and `V=Rgt\B`.  Since `C` covers every support edge,

```math
R_{UV}=0,
\qquad
R=\begin{pmatrix}R_{AB}&R_{AV}\\R_{UB}&0\end{pmatrix}.     \tag{VC.8}
```

For assignments `a` on `A` and `beta` on `B`, define

```math
Q_F^C(a,\beta)
=a^TR_{AB}\beta
 +\max_{u\in\{-1,1\}^{U}}
   \{F(a,u)+u^TR_{UB}\beta\}.                              \tag{VC.9}
```

### Theorem VC.2 (exact cover factorization)

For `y=(beta,v)`,

```math
(P_RF)(\beta,v)
=\max_{a\in\{-1,1\}^{A}}
 \{Q_F^C(a,\beta)+a^TR_{AV}v\}.                             \tag{VC.10}
```

Consequently the fragment has an exact composable representative with

```math
2^{|A|+|B|}=2^k                                             \tag{VC.11}
```

real coefficient slots.  A global baseline can be carried separately, leaving
a normalized projective shape with at most `2^k-1` real degrees of freedom.

#### Proof

For `x=(a,u)` and `y=(beta,v)`, (VC.8) gives

```math
x^TRy
=a^TR_{AB}\beta+a^TR_{AV}v+u^TR_{UB}\beta.                 \tag{VC.12}
```

In (VC.2), first maximize `u` at fixed `(a,beta)`, producing (VC.9), and then
maximize `a`, producing (VC.10). `square`

The representative is composable in the strong replacement sense.  For every
future `G`, gluing gives

```math
\max_{a,\beta,v}
\{Q_F^C(a,\beta)+a^TR_{AV}v+G(\beta,v)\}.                   \tag{VC.13}
```

Thus the original private variables may be discarded.  The table `Q`, the
known residual bridge `R_AV`, and the private auxiliary assignment `a` form a
replacement fragment.  Associativity of maximum and addition makes the same
replacement valid inside every later max-sum contraction.

For a fixed numerical `R`, all `2^k` coefficients need not be observable.
Define the envelope operator

```math
(\mathcal E_CQ)(\beta,v)
=\max_a\{Q(a,\beta)+a^TR_{AV}v\}.                            \tag{VC.14}
```

The genuinely coarsest exact coefficient state is

```math
Q\sim_RQ'\quad\Longleftrightarrow\quad
\mathcal E_CQ=\mathcal E_CQ',                               \tag{VC.15}
```

or equality up to a constant projectively.  Coefficients that never appear on
the finite Boolean upper envelope can be deleted.  Therefore `2^k` is a
universal support-based upper bound, not a claim of pointwise minimality for
every degenerate matrix or landscape.

Taking a minimum vertex cover gives the support width

```math
\tau(R)=\min\{|C|:C\text{ covers }H_R\},                    \tag{VC.16}
```

and an exact representative with `2^tau(R)` slots.

## 4. Worst-case exact minimality

The cover bound cannot be improved as a support-only guarantee.

### Theorem VC.3 (matching selector)

Let `L=Rgt=[k]`, let

```math
R=W I_k,\qquad 2W>D,                                        \tag{VC.17}
```

and let `Q:{-1,1}^k -> [-D,0]` be arbitrary with `max Q=0`.  Regard `Q` as
the child landscape.  Then

```math
(P_RQ)(y)=kW+Q(y)\qquad(y\in\{-1,1\}^k).                   \tag{VC.18}
```

Hence arbitrary futures recover every coordinate of `Q`, the normalized
response is exactly `Q`, and the projective state has `2^k-1` independent
real coordinates.

#### Proof

For any `x,y`,

```math
x^TRy=W(k-2d_H(x,y)).                                       \tag{VC.19}
```

The intended state `x=y` scores at least `kW-D`.  Every other state has
Hamming distance at least one and scores at most `kW-2W`.  The strict
inequality `2W>D` makes `x=y` the unique optimizer and proves (VC.18).
Proposition VC.1 and coordinate-pinning futures prove minimality. `square`

The support is a matching of `k` disjoint edges.  Every vertex cover needs at
least one distinct endpoint of every edge, and choosing all left endpoints
attains this bound, so its minimum cover size is exactly `k`.

If normalized coefficients are restricted to

```math
\{0,-\eta,\ldots,-M\eta\},\qquad D=M\eta,                   \tag{VC.20}
```

then the matching bridge has exactly

```math
(M+1)^{2^k}-M^{2^k}                                        \tag{VC.21}
```

distinct normalized response classes: all `2^k`-entry grid tables except
those with no zero entry.  Thus even the finite exact-state count is sharp.

## 5. Sharp approximate rate

Assume a normalized cover table has spread at most `D`:

```math
\max Q=0,\qquad -D\le Q(a,\beta)\le0.                       \tag{VC.22}
```

The envelope is sup-norm nonexpansive,

```math
\|\mathcal E_CQ-\mathcal E_CQ'\|_\infty
\le\|Q-Q'\|_\infty,                                       \tag{VC.23}
```

because every maximized term changes by at most the coefficient error.
Proposition VC.1 then shows that an exact future does not amplify a one-shot
table error.

Let `N=2^k`.  Gridding the `N` faces `max Q=0` gives a projective
`epsilon`-response code with at most

```math
N(\lceil D/\epsilon\rceil+1)^{N-1}                          \tag{VC.24}
```

codewords.  Conversely, apply the matching selector and fix one coefficient
to zero while placing every other coefficient on a `5epsilon`-grid in
`[-D,-D/2]`.  Two tables have a zero coordinate difference and another
difference of magnitude at least `5epsilon`; by (VC.6) their projective
distance is at least `5epsilon/2>2epsilon`.  For `0<epsilon<D/20`, this gives

```math
\left\lfloor{D\over10\epsilon}\right\rfloor^{N-1}
\le \mathcal N_{\rm proj}(\epsilon)
\le N(\lceil D/\epsilon\rceil+1)^{N-1}.                     \tag{VC.25}
```

Therefore the worst-case semantic rate is

```math
\log_2\mathcal N_{\rm proj}(\epsilon)
=(2^k-1)\log_2(D/\epsilon)+O(2^k)                           \tag{VC.26}
```

up to endpoint conventions and constants.  This is a genuine lossy theorem:
the distortion is tested by every future and the matching supplies a converse,
not merely an algorithm that runs exact dynamic programming on a smaller
instance.

The spread assumption is necessary.  With arbitrary unbounded real
coefficients, no finite uniform-error code exists.  Likewise no strict
contraction can be inferred from a small cover: the matching selector is an
isometry.  Requantizing at several compositions may therefore accumulate new
errors even though a single stored error survives every exact future unchanged.

For a fixed matrix the rate may be much smaller when many coefficients are
permanently dominated or several rows induce the same Boolean affine form.
The more precise instance parameter is the metric entropy of the reachable
envelope image in (VC.14); `tau(R)` is the sharp graph-support-only bound.

## 6. Why the degree-one obstruction is consistent

The degree-one bridge `R=I_k` has `k` disconnected edges, maximum degree one,
and treewidth one when considered in isolation.  Nevertheless its live cut
contains one independent channel per edge, and its minimum vertex cover has
size `k`, not one.  The bound above therefore predicts `2^k` coefficient
slots and exponential-in-`k` response bits.

The stronger macroscopic construction uses a Hamming code
`C subset {-1,1}^k`.  At query `y=c`, the matching score loses twice the
Hamming distance, so a bounded bonus attached independently to each
`c in C` remains exposed.  Since `|C|=2^{Omega(k)}`, arbitrary binary bonuses
give `2^{2^{Omega(k)}}` separated response functions and require
`2^{Omega(k)}` bits.  This agrees with (VC.26): the cover parameter itself is
extensive.

There is no contradiction with small-treewidth dynamic programming.  Such an
algorithm contracts matching edges one at a time only when the factors on the
two sides also decompose compatibly.  Here the child landscape is explicitly
arbitrary and may couple every left spin, while the future may inspect every
right assignment.  The bridge graph alone is not the full factor graph, and
its isolated treewidth is not the size of the declared live interface.

For bipartite graphs, the classical Koenig theorem identifies minimum vertex
cover size with maximum matching size.  Maximum matching is therefore an
equivalent support parameter and makes the obstruction particularly visible;
it is not a smaller width.

## 7. Exact verifier

Run

```bash
python3 extremal_information/experiments/verify_vertex_cover_bridge_response.py
```

The verifier uses integer and rational arithmetic.  It checks the block-cover
factorization against direct maximization for arbitrary small landscapes,
contextual pinning, gauge composition, matching selection and cover size,
the finite lattice count, envelope nonexpansiveness, and matching isometry.

## 8. Post-freeze classical comparison

Only after freezing (VC.9)--(VC.26) was the construction compared with the
literature and the existing bridge notes.

Classically, (VC.9) is cutset/separator conditioning: retain a table on a live
set and eliminate variables separated from the future.  Bucket and junction-
tree algorithms similarly create a function on the current separator, with
time and space exponential in its width.  Recursive conditioning makes the
same decomposition principle explicit as a time--space tradeoff.  Tensor-
network contraction likewise relates contraction cost to a width of the
network rather than to maximum degree.

Primary sources:

- [Dechter, *Bucket elimination: A unifying framework for
  reasoning*](https://ics.uci.edu/~dechter/publications/r76A.pdf)
- [Darwiche, *Recursive
  conditioning*](https://doi.org/10.1016/S0004-3702(00)00069-2)
- [Markov--Shi, *Simulating quantum computation by contracting tensor
  networks*](https://arxiv.org/abs/quant-ph/0511069)
- [Bodlaender et al., *Parameterized Complexity of Binary CSP: Vertex Cover,
  Treedepth, and Related Parameters*](https://arxiv.org/abs/2208.12543)

The classical comparison validates the `2^k` separator-table form.  The added
benchmark conclusions are semantic: (i) arbitrary futures identify the exact
envelope quotient, (ii) matching selectors prove worst-case minimality of the
cover table, and (iii) the contextual metric entropy has the sharp rate
(VC.26).  No priority claim is made for the decomposition itself.

## 9. Benchmark verdict

**Pass, independently predicted.**  A support vertex cover of size `k` gives
an exact `2^k`-coefficient tropical replacement for every internal child
landscape.  Rich futures show that the coefficient count is worst-case sharp,
and the same selector produces a matching approximate rate theorem.  The
degree-one obstruction is not an exception: its cover, matching number, and
live semantic interface are all extensive.
