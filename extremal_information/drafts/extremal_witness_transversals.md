# Extremal witness transversals and coordinate-query complexity

**Status.** Rigorous task-local draft.  This note quantifies a resource left
open by the exact coordinate compiler: the number of coordinate contexts
needed to expose every member of a switched family.

The result is a group-covering theorem, not a lower bound for arbitrary
quadratic futures.  It nevertheless shows that exact linear-order metric
compilation and a small reusable query language are logically independent.

## 1. Witness transversals

Let `G={+-1}^k` under coordinatewise multiplication and let
`f:G->R` have `Q=max_x|f(x)|>0`.  For `0<alpha<=1`, put

```math
W_alpha(f)=\{x:|f(x)|\ge\alpha Q\},
\qquad p_alpha={|W_alpha(f)|\over2^k}.               \tag{WT.1}
```

For every switch `s in G`, write `f_s(x)=f(sx)`.  A coordinate library
`X subset G` is `alpha`-extremal for the orbit if

```math
\max_(x\in X)|f_s(x)|\ge\alpha Q
\quad\hbox{for every }s\in G.                       \tag{WT.2}
```

### Theorem WT.1 (query size is extremal transversal number)

Every `alpha`-extremal coordinate library satisfies

```math
|X|\ge {1\over p_alpha}.                            \tag{WT.3}
```

Conversely, there is such a library with

```math
|X|\le
\left\lceil{k\log2+1\over p_alpha}\right\rceil.    \tag{WT.4}
```

Thus, up to a factor linear in `k`, the coordinate-query complexity of the
full switching orbit is the reciprocal mass of one extremal witness set.

#### Proof

Condition (WT.2) says that `X` meets every translate `sW_alpha`.  Equivalently,
for every `s` there are `x in X,w in W_alpha` with `w=sx`, so

```math
G=W_alpha X.
```

The product-set bound `|W_alpha X|<=|W_alpha||X|` proves (WT.3).

For the converse, sample `q` independent uniform elements of `G`.  A fixed
translate is missed with probability at most `exp(-p_alpha q)`.  A union
bound over the `2^k` translates is below one when
`q>(k log2)/p_alpha`, proving (WT.4). `square`

This is an operational identity rather than a new name for entropy: the
same tail mass that controls a one-landscape rare event determines, within a
factor `k`, the exact number of coordinate futures needed after all switches
are declared.

## 2. An exponential flat quadratic example

Let `k=s^2`, with `s` even, and partition the coordinates into `s` blocks of
size `s`.  Put

```math
D_(ij)=\begin{cases}
2,&i,j\hbox{ are distinct vertices in one block},\\
0,&\hbox{otherwise}.
\end{cases}                                        \tag{WT.5}
```

Then, with block magnetizations `M_b`,

```math
f(x)=H_D(x)=\sum_(b=1)^s(M_b^2-s),
\qquad Q=s^3-s^2.                                   \tag{WT.6}
```

### Theorem WT.2 (flat switched quadratics need exponentially many pins)

For every fixed `0<alpha<1`,

```math
p_alpha\le\exp\{-\alpha k/4+O(\sqrt k)\}.          \tag{WT.7}
```

Consequently every coordinate-pin library exposing an `alpha` fraction of
the absolute maximum of every switch `f_s` has

```math
|X|\ge\exp\{\alpha k/4-O(\sqrt k)\}.                \tag{WT.8}
```

The example may be realized as differences of spectrally flat exact
signings.  There are exact hollow sign matrices `A,A'` such that

```math
A-A'=D,
\qquad ||A||_(2->2)+||A'||_(2->2)=O(\sqrt k).       \tag{WT.9}
```

All switched pairs `(D_sAD_s,D_sA'D_s)` retain these bounds and have the
same target-scale contrast `f_s`.

#### Proof

Equation (WT.6) is immediate, and the negative tail is absent for large `s`
because `f>=-s^2` whereas `alpha Q` has order `s^3`.  For
`0<lambda<1/2`, Gaussian integration and
`cosh t<=exp(t^2/2)` give

```math
E\exp\{\lambda M_b^2/s\}
\le(1-2\lambda)^(-1/2).                             \tag{WT.10}
```

Independence of the blocks and Chernoff with `lambda=1/4` yield

```math
Pr\{f\ge\alpha(s^3-s^2)\}
\le\exp\{-\alpha s^2/4+O(s)\},                    \tag{WT.11}
```

which is (WT.7).  Theorem WT.1 gives (WT.8).

For (WT.9), prescribe `A=+1` and `A'=-1` on every within-block edge, and
give them the same symmetric random signs on all between-block edges.  The
within-block matrices have operator norm `s-1`.  A standard symmetric
Rademacher matrix bound gives one common between-block completion of norm
`O(s)`.  The triangle inequality proves (WT.9). `square`

## 3. Coding benchmark: ordinary covering codes reappear

The same theorem independently recovers a classical composable state in a
different model.  Identify `G` with the binary Hamming cube and take

```math
f_s(x)=d_H(x,s),
\qquad Q=k.                                         \tag{WT.12}
```

For `alpha>1/2`, a library `X` satisfies

```math
\max_(x\in X)d_H(x,s)\ge\alpha k
\quad\hbox{for every }s                             \tag{WT.13}
```

if and only if the antipodal library `-X` is a binary covering code of
radius `floor((1-alpha)k)`.  Therefore the sphere-covering bound and the
standard random-cover upper bound give

```math
{2^k\over\sum_(j\le(1-alpha)k){k\choose j}}
\le |X|
\le
\left\lceil{(k\log2+1)2^k\over
 \sum_(j\le(1-alpha)k){k\choose j}}\right\rceil.   \tag{WT.14}
```

In particular,

```math
\log|X|=k\{\log2-h(1-alpha)\}+O(\log k).            \tag{WT.15}
```

This is exactly WT.1 because the extremal witness set is one Hamming ball
about the antipode.  Thus reciprocal extremal mass is not merely a quadratic
artifact: on the simplest code-distance landscape it is the classical
covering-code rate.

## 4. Meaning for response compression

The rank-one coordinate compiler maps the response at query `x` to a
disjoint complete-sign parent with only twice the child order.  WT.2 shows
that, even for bounded-operator exact sign children, retaining a fixed
fraction of every switched contrast through **coordinate pins** may require
exponentially many predeclared contexts.

This does not imply that an arbitrary nonlinear future language needs that
many contexts.  A single context may aggregate many coordinate witnesses.
It does prove that the exponential language in the exact pin construction is
not removable by a generic spectral-flatness or large-maximum argument.  A
smaller compiler must exploit a different query algebra, not merely sample
fewer old configurations.
