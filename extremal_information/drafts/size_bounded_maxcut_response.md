# Size-bounded Max-Cut response entropy

**Status.** Audited benchmark draft. The proof uses classical hyperplane-
arrangement and finite-dimensional covering arguments. No novelty is claimed
for those mechanisms. The project-specific conclusion is that component size,
unlike boundary load alone, bounds coarse future-response information even
when edge weights have unlimited precision.

## 1. Resource model

Let

```math
X_w=\{+1,-1\}^w/\{s\sim-s\}
```

with projective Hamming distance. A pure weighted Max-Cut component `G` has
labelled boundary `[w]`, at most `p` private vertices, at most `m` nonzero
edges, and arbitrary nonnegative real edge weights. Its conditional response
is

```math
h_G([s])=\max_z\operatorname {Cut}_G(s,z).
```

The outer load at boundary vertex `i` is the total weight of incident edges.
We impose unit load at every boundary vertex. Responses are compared modulo
constants, in the operational shape metric

```math
d_{sh}([f],[g])={1\over2}\operatorname {osc}(f-g). \tag{1.1}
```

Quotienting is necessary: a private component disjoint from the boundary can
add an arbitrarily large constant without changing any boundary behavior.
Write `C_(w,p,m)` for this response-shape class.

The resources are semantic: `(m,p)` counts shared weight parameters and
private optimization variables. No weight-magnitude or precision bound is
assumed.

## 2. Shared-parameter arrangement theorem

### Proposition 2.1 (max-affine presentation entropy)

Let `X` have `q` elements. For each `x in X`, let `A_x subset R^d` have at
most `r` elements, and put

```math
f_\theta(x)=\max_{a\in A_x}\langle a,\theta\rangle,
\qquad \theta\in\Theta\subseteq\mathbb R^d.       \tag{2.1}
```

Suppose all `[f_theta]` lie in the radius-`R` ball of
`R^X/R1` under `d_sh`. Then, with internal cover centers,

```math
\log_2\operatorname {Cov}_\delta\{[f_\theta]:\theta\in\Theta\}
\le d\log_2\!\bigl(4(qr^2+1)\bigr)
 +d\log_2\!\left(1+{2R\over\delta}\right).       \tag{2.2}
```

If every witness vector is binary, `A_x subset {0,1}^d`, then the first
term can be replaced by

```math
d\log_2\!\left(4\left({3^d-1\over2}+1\right)\right), \tag{2.3}
```

so the entropy is

```math
O\left(d^2+d\log\left(1+{R\over\delta}\right)\right), \tag{2.4}
```

independently of `q`, `r`, and the precision of `theta`.

#### Proof

For every `x` and pair `a,b in A_x`, introduce the comparison hyperplane

```math
\langle a-b,\theta\rangle=0.                     \tag{2.5}
```

There are at most `q binom(r,2)` listed hyperplanes. On every relatively open
face of their arrangement all comparisons have fixed signs. A fixed
tie-breaking rule therefore selects one maximizer for every `x`, and (2.1)
is one linear map on that face.

An arrangement of `H` hyperplanes in `R^d` has at most

```math
\sum_{j=0}^d2^j\binom Hj\le[4(H+1)]^d            \tag{2.6}
```

faces, including lower-dimensional tie faces. For the image of one face,
take a maximal `delta`-separated subset in its quotient linear span. The
disjoint `delta/2` norm balls lie in the ball of radius `R+delta/2`, giving
at most `(1+2R/delta)^d` centers. This argument does not require the image to
be convex or the parameter set to be bounded.

If the witnesses are binary, every nonzero comparison normal belongs to
`{-1,0,1}^d`, and opposite normals give the same hyperplane. Hence there are
at most `(3^d-1)/2` distinct hyperplanes, proving (2.3). `square`

## 3. Pure-Max-Cut corollary

### Theorem 3.1 (precision-free size bound)

For `w>=2`, define

```math
H_*:=\min\left\{2^{w+2p-2},{3^m-1\over2}\right\}. \tag{3.1}
```

Then

```math
\begin{aligned}
\log_2\operatorname {Cov}_\delta(C_{w,p,m})
\le{}&\log_2(m+1)+2m\log_2(w+p)\\
&+m\log_2(4(H_*+1))
+m\log_2\left(1+{w\over2\delta}\right).         \tag{3.2}
\end{aligned}
```

Consequently, for fixed `epsilon>0`,

```math
\log_2\operatorname {Cov}_{\epsilon w}(C_{w,p,m})
=O_\epsilon\bigl(m^2+m\log(w+p)\bigr).           \tag{3.3}
```

#### Proof

Fix a labelled topology with `e<=m` edges and weight vector `c`. For boundary
word `s` and private assignment `z`, let `a_(s,z) in {0,1}^e` be its cut-
incidence vector. Then

```math
h_c(s)=\max_z\langle a_{s,z},c\rangle.           \tag{3.4}
```

Proposition 2.1 applies. The direct count uses
`|X_w| binom(2^p,2)<2^(w+2p-2)` comparisons; the ternary-normal count uses at
most `(3^e-1)/2` distinct hyperplanes.

Unit outer load makes `h_c` one-Lipschitz for projective Hamming distance:
hold the private cut fixed while boundary spins are changed, then compare
maxima in both directions. Since the projective diameter is at most `w/2`,

```math
\|[h_c]\|_{sh}\le w/4.                           \tag{3.5}
```

There are at most `(m+1)(w+p)^(2m)` labelled simple topologies. Parallel
edges have first been combined and zero edges deleted. Unioning the covers
from Proposition 2.1 over these topologies proves (3.2).

Modulo constants, delete every private connected component not meeting the
boundary. A spanning forest in each remaining component shows `p<=m`.
Thus one may replace (3.3) by

```math
\log_2\operatorname {Cov}_{\epsilon w}(C_{w,m})
=O_\epsilon\bigl(m^2+m\log(w+m)\bigr).           \tag{3.6}
```

Unbounded internal weights cause no gap: parameter space may be unbounded,
but (3.5) bounds its response image in the quotient, which is the only set
used in the volumetric argument. `square`

## 4. Exponential edge cost of a universal compiler

The normalized distance-shell theorem realizes the entire class

```math
L_w=\operatorname {Lip}_1(X_w,d_{proj})/\mathbb R \tag{4.1}
```

by unit-load Max-Cut components, without a size bound. Its Hamming-code
packing gives, for fixed `0<eta<1/4`,

```math
\log_2\operatorname {Cov}_{\eta w}(L_w)
\ge2^{(1-H_2(2\eta)+o(1))w}.                     \tag{4.2}
```

If the `m=m(w)` edge class were an `epsilon w`-net of `L_w`, then for every
fixed `tau>0` with `epsilon+tau<1/4`, a `tau w`-cover of that class would be
an `(epsilon+tau)w`-cover of `L_w`. Equations (3.6) and (4.2) imply

```math
\liminf_{w\to\infty}{\log_2m(w)\over w}
\ge {1-H_2(2\epsilon)\over2},
\qquad 0<\epsilon<1/4.                            \tag{4.3}
```

Thus polynomial-size graphs cannot approximate all unit-load responses at
macroscopic accuracy, even with arbitrary real weights.

## 5. Why a bare bit budget is not the right resource

A bare `B`-bit codebook has at most `2^B` outputs, and no stronger uniform
statement is possible without restricting its decoder. Select any `2^B`
members of a separated subset of `L_w`; the distance-shell construction
realizes all of them at unit load, possibly with exponential component size.
A hardwired decoder therefore attains the cardinality bound.

For a conventional graph encoding with `b` bits per weight, the bound
`O(m log(w+p)+mb)` is ordinary representation counting. The nontrivial
content of Theorem 3.1 is that `(m,p)` bounds coarse response information
even when `b` is infinite. The useful global resource is therefore a
restricted max-affine presentation or grammar, not an abstract bit label.

## 6. Second-model validation

The same theorem applies to a fixed weighted Boolean-CSP topology with `m`
weighted predicates. Conditional on an interface assignment, every private
assignment produces a satisfaction-incidence vector in `{0,1}^m`; its
conditional optimum is (2.1). Any promise bounding response oscillation by
`2R` therefore gives the precision-free entropy ceiling (2.4), plus the
logarithm of the number of allowed topologies.

Equivalently, consider a finite acyclic max-plus network with `m` weighted
arcs and any number of declared source-terminal queries. Each path has a
binary arc-incidence vector, so all path-value responses together obey

```math
\log_2\operatorname {Cov}_\delta
=O\left(m^2+m\log\left(1+{R\over\delta}\right)\right) \tag{6.1}
```

whenever their shapes have radius at most `R`. The number of paths and
queries does not enter. This is not a new shortest/longest-path algorithm;
it validates that shared finite-alphabet presentations, rather than raw
landscape cardinality, control approximate future-response information.

## 7. Remaining sharp question

The `m^2` term is the generic arrangement ceiling for ternary normals. It is
unknown here whether Max-Cut realizes `Omega(m^2)` *robustly separated*
response-cell entropy or whether its cut structure forces an
`O_epsilon(m log m)` law. Resolving that dichotomy is the next non-counting
question.
