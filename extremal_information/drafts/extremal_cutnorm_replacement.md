# Extremal cut-norm replacement

Status: literature-grounded orthogonal theorem, independently checked at the
level of constants used below.  Portfolio judgment: **promote narrowly** for
dense leading-scale replacement.

## Uniform all-future replacement

Let `A,B` be real matrices on the same labeled vertex interface and use the
unnormalized cut norm

```math
||A-B||_square=max_(S,T subseteq V)|(A-B)(S,T)|.              \tag{CR.1}
```

For `q` labels, a pair reward `J in R^(q times q)`, and an arbitrary future
conditional response `F:[q]^V->R`, put

```math
M(A;J,F)=max_(sigma in [q]^V)
\left\{F(sigma)+sum_(u,v)A_(uv)J_(sigma_u,sigma_v)\right\}.    \tag{CR.2}
```

### Theorem CR.1 (labeled cut-norm extremal replacement)

```math
|M(A;J,F)-M(B;J,F)|
<=||J||_1||A-B||_square.                                      \tag{CR.3}
```

#### Proof

For fixed `sigma`, let `V_i=sigma^(-1)(i)`.  Its two pair energies differ by

```math
sum_(i,j)J_(ij)(A-B)(V_i,V_j),
```

whose absolute value is at most the right side of (CR.3).  Apply the pointwise
bound to an optimizer for each matrix. `square`

The future `F` can be the optimum of an arbitrarily large private fragment,
or it can pin one exceptional labeling.  Thus (CR.3) preserves rare
extremizers, not only a typical energy distribution, and assumes no knowledge
of the optimizer.

For bounded dense matrices, Frieze--Kannan weak regularity supplies a
`k`-block average representative with

```math
k<=2^(O(epsilon^(-2))),
\qquad ||A-B||_square<=epsilon n^2.                            \tag{CR.4}
```

Its response depends only on the `k times q` occupancy table.  Independent
rounding of fixed block weights adds only `O(n^(3/2))=o(n^2)` cut error, so
the representative can be realized by finite simple weighted/unweighted
graphs at the dense leading scale.

Primary sources:

- [Frieze--Kannan, *Quick Approximation to Matrices and Applications*](https://doi.org/10.1007/s004930050052)
- [Borgs--Chayes--Lovasz--Sos--Vesztergombi, *Convergent Sequences of Dense
  Graphs II*](https://annals.math.princeton.edu/wp-content/uploads/annals-v176-n1-p02-p.pdf)
- [Braides--Cermelli--Dovetta, *Gamma-limit of the cut functional on dense
  graph sequences*](https://www.numdam.org/item/10.1051/cocv/2019029.pdf)

## Benchmark with an exponentially large exact response

Let `A_uv` be independent continuous `[0,1]` weights on a symmetrized complete
graph and let `h_A(S)` be its cut weight.  Almost surely, all
`2^(n-1)` projective cut coordinates are distinct, and ordinary private
Max-Cut attachments can expose them.  Nevertheless the constant
representative `B_uv=1/2` obeys

```math
||A-B||_square=O(n^(3/2))                                     \tag{CR.5}
```

with high probability, by a rectangle union bound.  Hence, uniformly over
every future conditional profile,

```math
max_S\{h_A(S)+F(S)\}
=max_S\left\{{1\over2}|S|(n-|S|)+F(S)\right\}+o(n^2).        \tag{CR.6}
```

The exact contextual state is exponentially large, while one density and
one magnetization count suffice at its dense normalized scale.  This is a
benchmark where replacement succeeds and finite exact quotienting does not.

## A checkable falsifier

For a proposed representative set `D=A-B`.  The Alon--Naor Grothendieck SDP
finds rectangles `S,T` with

```math
|D(S,T)|>=rho||D||_square,
\qquad rho>0.56.                                               \tag{CR.7}
```

Encode the four-state label
`sigma_v=(1_S(v),1_T(v))` and use, for symmetric `D`,

```math
J_((a,b),(c,d))={ad+bc\over2}.                                \tag{CR.8}
```

The energy difference at this labeling is exactly `D(S,T)`.  A unary future
`-K d_H(sigma,sigma^*)`, with `K` larger than every possible one-vertex pair
energy change, pins it for both matrices.  Their optimized responses then
differ by exactly `|D(S,T)|`.  Thus a cut residual of order the declared
scale is an explicit future-response falsifier, not merely failure of one
sufficient premise.

Primary source: [Alon--Naor, *Approximating the Cut-Norm via Grothendieck's
Inequality*](https://doi.org/10.1137/S0097539704441629).

## Scale boundary

The theorem is useful only after naming the leading scale `L_n`:

```math
||A_n-B_n||_square=o(L_n)                                     \tag{CR.9}
```

is required.  At `L_n=n^2`, weak regularity gives a finite block state for
fixed accuracy.  At `L_n=n^(3/2)`, it requires
`epsilon=o(n^(-1/2))`; the generic bound in (CR.4) then permits exponentially
many blocks.  Ordinary cut regularity therefore does not compress the
motivating signing scale.  Promotion beyond the dense regime requires a new
structured signed replacement lemma with subexponential state at
`o(n^(3/2))` cut error.
