# Grassmannian packing for macroscopic syndrome responses

**Status.** Proved and independently audited draft.  The construction gives
a quadratic-bit packing for the unrestricted appended-fragment
covering-radius experiment.  It improves
the linear direct-sum packing of Theorem 8.3, but it does not match the
exponential-in-`w` description size of the current general upper net.

## 1. Response metric

Let `G=F_2^w`.  For a spanning support
`S subseteq G\{0}`, write

```math
\lambda_S(g)=\min\{|A|:A\subseteq S,\ \sum_{a\in A}a=g\},
\qquad
r(S)=\max_{g\in G}\lambda_S(g).                 \tag{GP.1}
```

Repeated binary columns never help, so this is the covering radius of any
full-rank parity-check fragment with distinct column support `S`.  The future
response and complete future-response metric are

```math
R_S(U)=r(S\cup U),
\qquad
d_{\rm fut}(S,T)=\sup_{U\text{ spanning}}|R_S(U)-R_T(U)|. \tag{GP.2}
```

The query support is allowed to depend on the pair being separated, as it is
in every uniform response-packing lower bound.

## 2. A dense carrier and a sparse quotient

Fix `w=2d`.  For every `d`-dimensional subspace `W<=G`, choose an arbitrary
linear complement `V_W`, and let `C_W` be a basis of `V_W`.  Define

```math
S_W=(W\setminus\{0\})\cup C_W.                  \tag{GP.3}
```

Thus `S_W` contains one dense `d`-dimensional carrier and only a basis in the
transverse directions.  It spans `G` and has

```math
|S_W|=2^d-1+d.                                   \tag{GP.4}
```

### Lemma GP.1 (self and cross radii)

For any `d`-subspaces `W,W'<=G`, with
`t=dim(W cap W')`,

```math
r(S_W)=d+1,                                      \tag{GP.5}
```

and

```math
r(S_W\cup S_{W'})\le t+2.                       \tag{GP.6}
```

Consequently the valid self-query `U=S_W` gives

```math
d_{\rm fut}(S_W,S_{W'})\ge d-t-1.               \tag{GP.7}
```

#### Proof

Write `G=W direct-sum V_W`.  A vector `g=a+b`, with `a in W` and
`b in V_W`, has a unique coordinate representation of `b` in `C_W`.
The `W`-component costs zero generators when `a=0` and exactly one generator
when `a!=0`.  Hence

```math
\lambda_{S_W}(a+b)=\operatorname{wt}_{C_W}(b)+1_{a\ne0}.
```

Taking `b` to be the sum of every member of `C_W` and taking `a!=0` proves
(GP.5).

For the union, put `L=W+W'`.  Its codimension is

```math
\operatorname{codim}L
=2d-(2d-t)=t.                                    \tag{GP.8}
```

The image of `C_W` spans `G/L`, because `W union C_W` spans `G`.  Choose
`t` members of `C_W` whose images form a quotient basis.  Any `g in G` can
therefore be changed into an element of `L` using at most `t` generators
from `C_W`.  Every element of `L=W+W'` is `a+a'` with `a in W` and
`a' in W'`, hence uses at most two further generators from the two dense
carriers.  This proves (GP.6).  In (GP.2), choose `U=S_W`.  Its response from
`S_W` is (GP.5), while its response from `S_W'` is bounded by (GP.6), proving
(GP.7). `square`

The proof does not require compatible choices of the complements.  The
macroscopic separation is carried by the relative position of two global
subspaces, not by a declared direct sum of constant-size independent blocks.

## 3. A quadratic-bit constant-dimension packing

Write `[n choose k]_2` for a Gaussian binomial coefficient.

### Lemma GP.2 (elementary Grassmann packing)

For every `d>=2` and integer `1<=r<=d`, there is a family `A_(d,r)` of
`d`-subspaces of `F_2^{2d}` such that

```math
d-\dim(W\cap W')\ge r\quad(W\ne W'),             \tag{GP.9}
```

and

```math
|A_{d,r}|\ge {2^{d^2-r(2d-r)}\over16d}.          \tag{GP.10}
```

When the right-hand side is below one, the assertion has only its trivial
meaning; its useful regime is `r<d` by a fixed linear margin.

#### Proof

There are

```math
N={2d\brack d}_2\ge2^{d^2}                       \tag{GP.11}
```

`d`-subspaces.  The inequality follows factor by factor from

```math
{2^{2d-i}-1\over2^{d-i}-1}\ge2^d.
```

For a fixed `W`, the number of `d`-subspaces `W'` with intersection dimension
exactly `j`, or injection distance `s=d-j`, is

```math
{d\brack j}_2^2 2^{(d-j)^2}.                     \tag{GP.12}
```

The elementary uniform estimate

```math
{n\brack k}_2\le4\,2^{k(n-k)}                   \tag{GP.13}
```

follows from the product formula and
`prod_(i>=1)(1-2^{-i})>1/4`.  Therefore (GP.12) is at most

```math
16\,2^{2j(d-j)+(d-j)^2}=16\,2^{d^2-j^2}.        \tag{GP.14}
```

In terms of `s=d-j`, the exponent in (GP.14) is

```math
d^2-(d-s)^2=s(2d-s).
```

This increases for `0<=s<=d`.  The total number with `s<r`, including `W`
itself, is therefore at most

```math
16d\,2^{r(2d-r)}.                                \tag{GP.15}
```

Greedily choose one subspace and delete this bad neighborhood.  Equations
(GP.11) and (GP.15) leave at least (GP.10) choices. `square`

For completeness, the product lower bound used in (GP.13) is elementary:
the first three factors have product `21/64`, and the tail is at least
`1-sum_(i>=4)2^{-i}=7/8`, so the infinite product exceeds `147/512>1/4`.

### Theorem GP.3 (quadratic macroscopic response information)

Let `w=2d>=4` be even and let `1<=r<=d`.  There is a family of full-rank
syndrome fragments, all of common length

```math
n_w=2^{w/2}-1+w/2,                               \tag{GP.16}
```

whose complete future-response maps contain a packing of size

```math
{2^{w^2/4-r(w-r)}\over8w}                        \tag{GP.17}
```

at pairwise response distance at least

```math
r-1.                                             \tag{GP.18}
```

Consequently, fix `epsilon<1/4` and choose any constant

```math
2\epsilon<\gamma<1/2.
```

For all sufficiently large even `w`, any deterministic summary which answers
every unrestricted appended-fragment covering-radius query with uniform
additive error at most `epsilon*w` requires at least

```math
\left((1/2-\gamma)^2-o(1)\right)w^2              \tag{GP.19}
```

bits in the worst case.

#### Proof

Use Lemma GP.2 and the supports (GP.3).  They have the common size (GP.16)
and are full rank.  For distinct packed subspaces, Lemma GP.1 gives

```math
d_{\rm fut}(S_W,S_W')\ge r-1,
```

which proves (GP.17)--(GP.18).  Two sources assigned the same deterministic
message have the same decoded answer to every query and therefore have true
response distance at most `2*epsilon*w`.  Put `r=ceil(gamma*w)`.  Since

```math
2\epsilon w<r-1
```

for all sufficiently large `w`, the summary must be injective on the packing.
Finally,

```math
{w^2\over4}-r(w-r)
=\left({1\over2}-\gamma\right)^2w^2-O(w),
```

so taking the base-two logarithm of (GP.17) proves (GP.19). `square`

Because `gamma` may be chosen arbitrarily close to `2*epsilon`, the
quadratic-rate coefficient may be made arbitrarily close to
`(1/2-2*epsilon)^2`.  This optimization is not asserted at the boundary,
where the strict decoder separation and integer rounding matter.

The same packing gives a Fano/mutual-information statement for expected
global sup-response distortion, exactly as in Theorem 8.3: put the uniform
prior on the packing and decode the closest response map.  We omit it because
the deterministic metric-entropy conclusion is the new point.

## 4. Scope and comparison with the previous block theorem

1. **This is not the old direct-sum bit source.**  The latent state is a point
   of a Grassmannian with `Theta(w^2)` separated choices in its logarithm.
   The exposing query is the complete support of one moving carrier.  The
   gap is the subspace-injection distance `d-dim(W cap W')`, up to one.

2. **Each source still has internal quotient structure.**  Its dense carrier
   `W` and sparse complement make the self-radius calculable.  Thus the result
   does not say that generic supports are hard; it says that arbitrary
   supports already contain a global, non-product, quadratic-bit packing.

3. **The fragments are exponentially long.**  Their length is
   `Theta(2^{w/2})`, unlike the linear-length fragments in Theorem 8.3.  The
   stronger information lower bound is therefore not free.  It is still far
   below the exact support state's `Theta(2^w)` bits and below the current
   general upper net at fixed normalized distortion.

4. **The construction transfers.**  In native language it is a packing of
   diameters of elementary-abelian Cayley graphs under future generator-set
   union.  The proof therefore applies verbatim to a query interface for
   finite-group word metrics: append generators and ask the diameter.  The
   dense-carrier argument uses the vector-space quotient and is specifically
   strongest for elementary abelian groups.

## 5. General-theory lesson

Composition can amplify the information in a **moving global carrier**.  A
single source has a dense low-cost subspace and a sparse expensive quotient.
Appending another carrier changes the quotient dimension by

```math
\dim W-\dim(W\cap W'),
```

so one scalar extremal response measures a macroscopic Grassmann distance.
The response complexity is therefore governed not only by how many local
features are present, but by the metric entropy of their possible rooted
positions inside the interface.

This is more than ordinary dynamic programming: the exact syndrome support
algebra was already known, but it did not imply that a one-parameter family
of legal future contexts embeds a constant-dimension-code metric into the
scalar radius response.  It is also not a generic rate--distortion restatement:
Lemmas GP.1--GP.2 supply the nontrivial response embedding and its packing.

The result does not settle the unrestricted dichotomy.  The remaining gap is
between quadratic response bits and an upper description with exponentially
many bits in `w`.  Reaching the latter scale would require exponentially many
macroscopically independent carrier positions without their unions collapsing
the Cayley diameter; the present one-carrier architecture cannot provide that.
