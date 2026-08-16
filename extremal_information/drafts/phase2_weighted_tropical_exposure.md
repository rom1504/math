# Query-weighted tropical exposure

**Status:** independently audited theorem draft.  The live version includes
the injective-anchor hypothesis required for disjoint witness rectangles.  It
turns the qualitative “mass is necessary” conclusion into a finite inequality
but is not promoted to the surface because its canonical code application has
vanishing macroscopic exposure.

## 1. Why a weighted theorem is needed

A tropical fooling set proves uniform-error rank lower bounds, but the
zero-diagonal/one-off-diagonal matrices show that the witness cells may have
vanishing query mass.  The correct average-error statement must combine:

1. the four-cell crossing gap;
2. the probability assigned to those four cells; and
3. the extent to which every `k`-channel assignment leaves disjoint
   same-channel witness pairs.

The next definition records exactly those three quantities.

## 2. A weighted monochromatic-matching parameter

Let `M` be a finite real matrix and let `I` index distinguished cells
`(x_i,y_i)`, with the `x_i` pairwise distinct and the `y_i` pairwise distinct.
Let `E` be a graph on `I`.  For `ij in E`, assume the crossing gap

```math
G_{ij}:=
M(x_i,y_j)+M(x_j,y_i)
-M(x_i,y_i)-M(x_j,y_j)>0.                           \tag{WE.1}
```

Let `mu` be a probability distribution on matrix cells.  For an edge `ij`,
write `C_ij` for its four witness cells

```math
C_{ij}=\{(x_i,y_i),(x_j,y_j),(x_i,y_j),(x_j,y_i)\}.
```

Define its exposure weight by

```math
w_{ij}=
{G_{ij}^2\over
 \displaystyle\sum_{c\in C_{ij}}1/\mu(c)},         \tag{WE.2}
```

with `w_ij=0` if any of the four cell masses is zero.  For `k>=1`, define

```math
\mathfrak m_k(M,E,\mu)=
\min_{c:I\to[k]}
\max\left\{
 \sum_{ij\in F}w_{ij}:
 F\text{ is a matching in }E,
 c(i)=c(j)\text{ for every }ij\in F
\right\}.                                         \tag{WE.3}
```

This is a weighted anti-coloring parameter.  It is zero precisely when some
`k`-color assignment can avoid every positively exposed matching witness.
The matching restriction prevents double counting matrix cells.

## 3. Average-error rank theorem

### Theorem WE.1 (query-weighted tropical exposure bound)

Let `Mtilde` have min-plus factorization rank at most `k`.  Then

```math
\boxed{
\mathbb E_{c\sim\mu}
[\widetilde M(c)-M(c)]^2
\ge \mathfrak m_k(M,E,\mu).}                       \tag{WE.4}
```

For max-plus factorization, apply the theorem to `-M` and `-Mtilde`; the
crossing gaps and `mathfrak m_k` are then recomputed for the reversed
four-cell contrast.  One must not reuse the min-plus weights of `M`
unchanged.

#### Proof

Choose, for every distinguished cell `(x_i,y_i)`, one factor term tight at
that cell.  Its term label defines a coloring `c:I->[k]`.

If `ij` is monochromatic, the separability/crossing calculation gives

```math
\sum_{z\in C_{ij}}
|\widetilde M(z)-M(z)|\ge G_{ij}.                  \tag{WE.5}
```

Indeed, the common term is tight at the two distinguished cells, majorizes
`Mtilde` at the crossed cells, and (WE.1) leaves exactly the four absolute
errors as slack.

Weighted Cauchy--Schwarz gives

```math
\sum_{z\in C_{ij}}\mu(z)
 [\widetilde M(z)-M(z)]^2
\ge
{G_{ij}^2\over\sum_{z\in C_{ij}}1/\mu(z)}
=w_{ij}.                                           \tag{WE.6}
```

For a matching `F`, the cell sets `C_ij` are disjoint: different edges have
disjoint endpoint indices, so neither their distinguished nor crossed cells
coincide.  Summing (WE.6) over a monochromatic matching bounds the total
mean-square error from below by its weight.  Maximize over such matchings and
then use the minimum over all `k`-colorings in (WE.3). `square`

### Interpretation

The theorem is neither a mutual-information bound nor a generic rank
continuity result.  It states the exact additional datum absent from an
unweighted tropical rank certificate: how much declared query mass remains
on disjoint witnesses after the factor terms partition the exposed anchors.

It also has an immediate falsifier.  If `mathfrak m_k=0`, no positive
mean-square lower bound follows from this witness system; another witness
graph or another query law is required.

## 4. The zero-density diagonal example

For the matrix `D_r` with zero diagonal and one off diagonal, take all
diagonal cells as anchors and the complete witness graph.  Every gap is two.
Under the uniform law on all `r^2` cells,

```math
w_{ij}={1\over r^2}.                                \tag{WE.7}
```

Any `k`-coloring has a monochromatic matching of size at least

```math
\left\lceil{(r-k)_+\over2}\right\rceil.            \tag{WE.8}
```

To see this, take a maximum matching inside every color class and sum
`floor(|C|/2)`.  At most one vertex per nonempty color class is left
unmatched, giving (WE.8).  Therefore

```math
\mathbb E_{\rm unif}(\widetilde D-D_r)^2
\ge {\lceil(r-k)_+/2\rceil\over r^2}.              \tag{WE.9}
```

For `k=1` this has the correct order `1/r`, matching the rank-one all-one
approximant up to a factor of two.  Thus the weighted theorem detects, rather
than conceals, the vanishing-mass obstruction.

## 5. What would make this generative

Theorem WE.1 reduces an average-error tropical lower bound to a weighted
monochromatic-matching problem.  Its value for the wider theory depends on
whether `mathfrak m_k` can be bounded from model structure without listing an
exponential response table.

The next falsifiable target is:

> For a natural trellis, code, or finite-CSP conditional response family and
> a non-adversarial declared query measure, construct a crossing witness graph
> whose exposure weights and monochromatic matching profile are computable
> from polynomial-size structural data and remain nonvanishing at the target
> distortion scale.

A proof would connect tropical feature growth to mean-square response
distortion.  A family for which every structurally describable witness graph
has `mathfrak m_k=o(1)` would show that average response complexity again
requires a different object.
