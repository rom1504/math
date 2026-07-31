# Edge-twisted Hadamard order multiplier

Status: explicit non-Kronecker candidate, exact finite-state audit, and a
proved scalable obstruction to **recursive use from two fixed seeds**.  The
one-step all-spin inequality for fresh growing-order minimizers remains open,
but its product channel already exposes a principal-restriction obligation.

## 1. Two incompatible order-four channels

Let

```math
H=
\begin{pmatrix}
1&1&1&1\\
1&-1&1&-1\\
1&1&-1&-1\\
1&-1&-1&1
\end{pmatrix},
\qquad
D=\operatorname{Diag}(-1,-1,-1,1),
\qquad K=DHD.                                         \tag{ET1}
```

Both `H` and `K` are symmetric Hadamard matrices.  Exact enumeration of the
16 Boolean vectors gives

```text
Ext(H) = {---+, -+++, +---, +++-},
Ext(K) = {----, -++-, +--+, ++++},
```

where `Ext(G)={v in {+/-1}^4:Gv=+/-2v}`.  In particular their Boolean
extremal sets are disjoint.

The disjointness remains strong after appending an arbitrary common
Hadamard factor.  Direct multiplication gives

```math
HK=
\begin{pmatrix}
2&2&2&2\\
2&2&-2&-2\\
2&-2&2&-2\\
-2&2&2&-2
\end{pmatrix},                                       \tag{ET2}
```

with characteristic polynomial

```math
(\lambda-4)^2(\lambda^2+4\lambda+16).                \tag{ET3}
```

There is no `-4` eigenspace, since `det(HK+4I)=1024`, and exact enumeration
shows that the `+4` eigenspace contains no Boolean vector.

Let `L` be any symmetric Hadamard matrix.  If a Boolean `x` were extremal for
both `H tensor L` and `K tensor L`, with eigenvalue signs `alpha,beta`, then

```math
(HK\otimes I)x=4\alpha\beta x.                       \tag{ET4}
```

View `x` as four rows, or equivalently as columns of length four indexed by
the `L` coordinate.  Every such column is Boolean and would have to lie in
the `4 alpha beta` eigenspace of `HK`.  Equations (ET2)--(ET3) and the exact
Boolean check exclude both signs.  Thus

```math
\boxed{\operatorname{Ext}(H\otimes L)
\cap\operatorname{Ext}(K\otimes L)=\varnothing.}     \tag{ET5}
```

This carefully handles all four choices of the two eigenvalue signs.

There is also an exact Boolean uncertainty statement, but it is weaker than
extremal-set disjointness might suggest.  Put

```math
\operatorname{Dom}(G)={x\in\{\pm1\}^4:Gx/2\in\{\pm1\}^4}.
```

Exact enumeration proves that both domains have size eight and partition the
Boolean cube.  Moreover `|x^TGy|<=8`, with the sharper bound `<=4` whenever
either endpoint is outside `Dom(G)`.  Thus every microstate is good for
exactly one channel, but a state bad for a channel still retains half of its
maximum bilinear response.  This finite-state fact explains why incompatibility
alone need not produce a small lift cap.

## 2. Explicit order-four multiplier

Let `A` be a signing of order `n>=4`.  Choose the deterministic balanced bit
vector

```math
b_i=0\quad(1\le i\le\lfloor n/2\rfloor),
\qquad b_i=1\quad(i>\lfloor n/2\rfloor),             \tag{ET6}
```

and write `G_0=H`, `G_1=K`.  Define the order-`4n` signing
`L_b(A)` in four-by-four blocks by

```math
(L_b(A))_{ij}=
\begin{cases}
G_{b_i}-\operatorname{diag}(G_{b_i}),&i=j,\\
a_{ij}G_{b_ib_j},&i\ne j.
\end{cases}                                          \tag{ET7}
```

The lower off-diagonal blocks are the transposes of the upper blocks.  Since
`H` and `K` are symmetric full sign matrices, (ET7) is symmetric, has zero
diagonal, and has a sign on every other entry.

The macro edge coloring uses `H` on every edge incident with the zero part
and `K` within the one part.  It therefore contains both channel types and is
not a common-block tensor `A tensor G`.  More importantly, (ET5) proves that
there is no Boolean microvector that can carry every old macro spin through
the multiplier by simultaneous Hadamard extremality.  The Boolean-channel
persistence theorem and its small-perturbation extension consequently do not
apply to (ET7) through the old seed state.

At successive orders one reapplies (ET7) with the new balanced partition of
all current vertices.  The edge coloring changes with the current scale; the
next matrix is not obtained by appending one common fiber generator to the
previous mosaic.

## 3. Exact all-spin inequality required

For a Boolean parent state `x=(x_1,...,x_n)`, `x_i in {+/-1}^4`, its energy is

```math
H_{L_b(A)}(x)
=\sum_{i<j}a_{ij}x_i^TG_{b_ib_j}x_j
 +{1\over2}\sum_i x_i^T
 (G_{b_i}-\operatorname{diag}G_{b_i})x_i.             \tag{ET8}
```

The precise target is the following uniform Boolean inequality, for an
absolute constant `C` and every exact minimizer `A` (or, more strongly, every
signing in an objective-independent structured family containing one):

```math
\boxed{
|H_{L_b(A)}(x)|\le8\operatorname{cap}(A)+Cn
\quad\hbox{for every }x\in\{+1,-1\}^{4n}.}           \tag{ET9}
```

Unlike a product-spin check, (ET9) quantifies over all states entangled
between the macro vertices and the two incompatible fiber channels.

The defect in (ET9) is strong enough for the convergence architecture.  Put
`p=cap(A)` and use the project's uniform lower scale `p>=c n^(3/2)`.
Concavity gives

```math
(8p+Cn)^{2/3}
\le4p^{2/3}+{2Cn\over3(8p)^{1/3}}
=4p^{2/3}+O(\sqrt n).                                \tag{ET10}
```

Thus (ET9) would give a `b`-scale multiplier

```math
b_{4n}\le4b_n+O(\sqrt n),                            \tag{ET11}
```

whose relative defect is `O(n^(-1/2))` and is geometrically summable under
iteration.  A separate bounded-residue order-filling operation would still
be needed to pass from the native multiplier orders to every integer order;
(ET9) supplies the previously missing entangled-spin part.

## 4. Exact finite audits: mixed evidence

Take the all-negative triangle, whose cap is three, and use

```math
b=(0,1,1).
```

Then its three cross blocks have channel colors `(H,H,K)` and its diagonal
blocks are `(H-diag H,K-diag K,K-diag K)`.  Exhaustive enumeration of all
`2^11` projective Boolean states proves

```math
\min H=-20,
\qquad\max H=20,
\qquad\operatorname{cap}(L_b(A))=20.                 \tag{ET12}
```

Consequently

```math
{20\over12^{3/2}}=0.481125224324688\ldots<\frac12,  \tag{ET13}
```

and also `20<8 cap(A)=24`.  Thus this particular mixed Boolean triangle does
**not** immediately recreate a normalized cap above `1/2`.  This is an exact
finite positive fact, not evidence that (ET9) holds uniformly.

The next two exact minimizer audits are negative.  For the canonical
root-gauged order-four minimizer (code 1, cap 4), all six balanced two/two
fiber partitions were enumerated.  Five lifts have cap 42 and one has cap 48:

```math
{42\over16^{3/2}}=0.65625,
\qquad {48\over16^{3/2}}=0.75.                       \tag{ET14}
```

For the canonical order-five minimizer (root-gauged code 13, cap 4), all ten
balanced two/three partitions have cap 56:

```math
{56\over20^{3/2}}=0.6260990337\ldots.                \tag{ET15}
```

These are exact exhaustive cap computations.  They show that merely mixing
the incompatible channels does not suppress finite entangled states.  They
do not yet falsify (ET9): the order-five excess is
`56-8*4=24=4.8n`, so it is consistent only with `C>=4.8`.

Larger exact witness checks are more negative.  From the certified order-ten
minimizer of cap 13, the saved order-40 state

```text
little-endian positive-bit hex e13ae1eeee
```

has exact energy `-170`, hence normalized absolute energy
`0.6719840027857806`.  Freshly twisting that order-40 signing again, the saved
order-160 state

```text
7889877787767889788987777888788878887888
```

has exact energy `1488`, hence normalized energy `0.7352295559891482`.
Its 40 fibers comprise 21 copies of `---+`, 15 of its antipode `+++-`, three
copies of `+--+`, and one of its antipode `-++-`.  Thus the large witness is
almost a two-type product channel, not a diffuse 16-state phenomenon.

For the saved order-14 conference-derived signing of cap 21, the order-56
state `07810007778967` has exact energy `238` and normalized energy
`0.5679301390639018`.  These states came from heuristic searches, but every
displayed energy is independently recomputed by exact integer arithmetic.

## 5. Exact product channel and scalable recursive obstruction

The near-product pattern has an exact explanation.  Set

```math
v=(-1,-1,-1,1).
```

Direct calculation gives

```math
v^THv=8,\qquad v^TKv=4,
\qquad {1\over2}v^T(H-\operatorname{diag}H)v=4,
\qquad {1\over2}v^T(K-\operatorname{diag}K)v=2.       \tag{ET16}
```

If `S={i:b_i=1}` and `z` is any parent spin, substitution of
`x_i=z_iv` in (ET8) gives the exact identity

```math
H_{L_b(A)}(z\otimes v)
=8H_A(z)-4H_{A[S]}(z_S)+4n-2|S|.                    \tag{ET17}
```

For a balanced split the last term is `3n`.  Hence even the simplest product
channel in (ET9) requires control of a signed combination of full and
principal-half energy.  Channel incompatibility did not eliminate the old
restriction-profile obligation; it merely assigned the restriction a
different coefficient.

There is a fully scalable obstruction to using this multiplier recursively
from a fixed seed.  Let `A_0=A` have even order `n`, let
`A_{r+1}=L(A_r)` with the aligned first/second-half partition, put
`N_r=n4^r`, and define `z_{r+1}=z_r tensor v`.  Write

```math
E_r=H_{A_r}(z_r),\qquad
R_r=H_{A_r[S_r]}((z_r)_{S_r}),                       \tag{ET18}
```

where `S_r` is the second half.  Equation (ET17) gives

```math
E_{r+1}=8E_r-4R_r+3N_r.                             \tag{ET19}
```

The important channel in the restriction recurrence is `K`, not `H`:
`A_{r+1}[S_{r+1}]` is the common-`K` lift of `A_r[S_r]`.  Therefore

```math
R_{r+1}=4R_r+N_r,qquad
R_r=4^r\left(R_0+{nr\over4}\right).                 \tag{ET20}
```

Dividing (ET19) by `8^(r+1)` and summing the geometric series proves

```math
\boxed{\lim_{r\to\infty}{E_r\over8^r}
=E_0-R_0+{n\over2}.}                                \tag{ET21}
```

The audit exhausts every projective base spin for the fixed balanced split.
For the order-ten exact minimizer it finds `E_0=11`, `R_0=-6`, giving

```math
\lim_{r\to\infty}{|E_r|\over N_r^{3/2}}
={22\over10^{3/2}}=0.6957010852\ldots .             \tag{ET22}
```

For the order-14 signing it finds `E_0=19`, `R_0=-5`, giving

```math
\lim_{r\to\infty}{|E_r|\over N_r^{3/2}}
={31\over14^{3/2}}=0.5917927499\ldots .             \tag{ET23}
```

This obstruction is not an artifact of the displayed vertex ordering.  The
audit also exhausts all `binom(10,5)=252` and `binom(14,7)=3432` choices of
the `K`-colored half, maximizing (ET21) over every projective spin for each
choice.  Even the best partition has worst limiting normalized energy

```math
{18\over10^{3/2}}=0.5692099788\ldots,
\qquad
{29\over14^{3/2}}=0.5536125725\ldots,                \tag{ET23a}
```

respectively.  Switching the seed merely bijects the maximized spin set, and
permutation merely changes the chosen subset, so (ET23a) covers all switching
and permutation preprocessing of these two seeds.

Explicit matrix construction independently verifies the first two steps:
the product energies are `142,1312` for the order-ten seed and `214,1904`
for the order-14 seed.  Thus fresh retwisting does **not** kill compatible
Boolean descendants; a fixed four-spin state survives through the aligned
`K` half, with a computable geometric correction.

This is a proved scalable obstruction to the proposed recursive construction
from these two seeds.  Its scope is precise: it does not disprove a one-step
inequality applied only to a newly chosen exact minimizer at each growing
order, because the recursively produced parents need not themselves be
minimizers.  Such a one-step theorem remains open, but (ET17) shows that it
must solve a principal-restriction correlation problem rather than follow
from the disjoint extremal channels alone.  Accordingly (ET7) is not a live
standalone composition mechanism absent a new restriction ingredient.

## 6. Immediate proof-route obstruction and current classification

Every order-`N` signing has Frobenius norm squared `N(N-1)`, hence operator
norm at least `sqrt(N-1)`.  For `N=4n`, the ordinary spectral cap certificate

```math
\operatorname{cap}(L_b(A))
\le {N\over2}\|L_b(A)\|_{op}                         \tag{ET24}
```

can therefore never have a right side asymptotically below
`(1/2)N^(3/2)`.  If the seed family has normalized cap below `1/2`, a proof of
(ET9) cannot come from operator norm, Frobenius norm, or a triangle inequality
that first bounds the two channel contributions separately.  It must use
the incompatibility in (ET5) directly at the Boolean level.

The Frobenius observation is only a proof-method obstruction.  Equations
(ET17)--(ET23), by contrast, are an exact scalable Boolean witness theorem.
Together with the finite caps (ET14)--(ET15) and saved heuristic witnesses,
they rule out further investment in the fixed recursive multiplier.  A
one-step growing-minimizer version of (ET9) should be reconsidered only if a
new theorem controls the full/principal energy combination in (ET17) without
reintroducing an equivalent restriction bridge.

## Reproduction

```bash
.venv/bin/python computations/audit_edge_twisted_hadamard_multiplier.py \
  --output computations/results/edge_twisted_hadamard_multiplier.json \
  > computations/logs/edge_twisted_hadamard_multiplier.log
```

The program checks both Hadamard identities, enumerates the extremal and
transport domains, checks the `HK` obstruction, exhausts the order-12,
order-16, and order-20 audits, verifies the saved order-40, order-56, and
order-160 witness energies, exhausts the base spins and every balanced
partition in (ET22)--(ET23a), and constructs the first two product descendants
to verify (ET19)--(ET21).
