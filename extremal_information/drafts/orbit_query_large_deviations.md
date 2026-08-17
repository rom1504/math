# Orbit-query complexity is a large-deviation rate

**Status.** Rigorous task-local draft.  This theorem joins two branches of
the program: contextual query complexity and rare-event compactness.  For a
translation orbit under repeated product composition, the exact exponential
rate of coordinate contexts is the Cramer rate of one extremal tail.

## 1. Product landscapes and their switched queries

Let `G_0` be a finite abelian group of order `q` and let
`f:G_0->R` be nonconstant.  On `G_n=G_0^n`, put

```math
F_n(x)=\sum_(i=1)^n f(x_i).                         \tag{LD.1}
```

For a threshold `a` with

```math
E f<a<\max f,                                      \tag{LD.2}
```

let

```math
W_n(a)=\{x:F_n(x)\ge an\},
\qquad p_n(a)={|W_n(a)|\over q^n}.                 \tag{LD.3}
```

A coordinate library `X_n subset G_n` answers every switched upper-tail
query at level `a` if

```math
\max_(x\in X_n)F_n(s+x)\ge an
\quad\hbox{for every }s\in G_n.                    \tag{LD.4}
```

Write `L_n(a)` for the smallest size of such a library.

### Theorem LD.1 (exact orbit-query exponent)

Define

```math
\Lambda_f(\theta)=
 \log\left({1\over q}\sum_(z\in G_0)e^(\theta f(z))\right),
\qquad
I_f(a)=\sup_(\theta\ge0)\{\theta a-\Lambda_f(\theta)\}.
                                                               \tag{LD.5}
```

Then

```math
{1\over p_n(a)}\le L_n(a)
\le\left\lceil{n\log q+1\over p_n(a)}\right\rceil,
                                                               \tag{LD.6}
```

and

```math
\boxed{\lim_(n\to\infty){1\over n}\log L_n(a)=I_f(a).}
                                                               \tag{LD.7}
```

Thus the macroscopic information rate of the future query language is
neither the number of configurations nor the maximum alone.  It is the
large-deviation cost of finding one threshold witness.

#### Proof

The translates of `W_n(a)` all have density `p_n(a)`.  A library satisfying
(LD.4) is exactly a transversal of every translate.  Product-set counting
and random sampling give (LD.6), as in Theorem WT.1.

It remains to compute the exponential tail.  Chernoff gives, for every
`theta>=0`,

```math
p_n(a)\le\exp\{-n(\theta a-\Lambda_f(\theta))\},    \tag{LD.8}
```

hence the lower bound `I_f(a)` on the decay rate.  For the reverse bound,
group equal values of
`f` into a finite alphabet.  A type with empirical law `nu` has probability

```math
\exp\{-nD(\nu||\mu)+O(\log n)\},                   \tag{LD.9}
```

where `mu` is the uniform pushforward law of `f`.  Minimize
`D(nu||mu)` subject to `E_nu f>=a` and approximate a minimizer by rational
types.  Finite-dimensional entropy duality identifies this minimum with
the Legendre transform in (LD.5).  Therefore

```math
-{1\over n}\log p_n(a)\longrightarrow I_f(a).      \tag{LD.10}
```

The polynomial factor in (LD.6) disappears after division by `n`, proving
(LD.7). `square`

## 2. Heterogeneous composition law

Let `f` and `g` be landscapes on finite groups and compose `m` copies of `f`
with `n` copies of `g`.  If `m/(m+n)->lambda` and the threshold `a` lies
strictly between the weighted mean and weighted maximum, the normalized
log-moment state is

```math
\Lambda_\lambda(\theta)
=\lambda\Lambda_f(\theta)+(1-\lambda)\Lambda_g(\theta),        \tag{LD.11}
```

and the switched-query exponent at level `a` is

```math
I_\lambda(a)=\sup_(\theta\ge0)
 \{\theta a-\Lambda_\lambda(\theta)\}.             \tag{LD.12}
```

If `L_(m,n)(a)` denotes the smallest coordinate library meeting every
translate of this heterogeneous upper tail, then

```math
\lim_(m+n\to\infty){1\over m+n}\log L_(m,n)(a)
=I_\lambda(a).                                      \tag{LD.12a}
```

Hence this response complexity has an exact composition algebra: log-moment
functions add, and the exposed query rate is their Legendre dual.  It is a
strict quotient of the full product landscape whenever the alphabet of
`f,g` is fixed.

The state is query-relative.  Retaining all thresholds `a` requires the
whole convex function `Lambda`; a single declared threshold needs only its
supporting tilt.  This distinction prevents the theorem from being mistaken
for a universal finite statistic.

## 3. Benchmarks

For a Rademacher base `f in {+-1}` with equal masses,

```math
I_f(a)=D\left({1+a\over2}\middle\|{1\over2}\right),
\qquad 0<a<1.                                      \tag{LD.13}
```

After complementing coordinates, `X_n` is precisely a binary covering code
of relative radius `(1-a)/2`; (LD.7) is the classical sphere-covering
exponent.

For a finite random-energy atom table, the same formula says that repeated
independent composition converts its one-block energy histogram into the
exact future-query information rate.  No overlap geometry is needed because
the allowed composition and switching action are coordinatewise products.
Adding cross-block couplings destroys (LD.11); that is exactly where a richer
extremal state becomes necessary.

## 4. Theory judgment

This theorem is generative but scoped.  It predicts a compositional response
complexity from a standard rare-event object, applies to code distance and
finite random-energy products, and explains why exponentially rare witness
sets force exponential query languages.  It does not preserve optimizer
geometry under interacting composition.  Its falsifier is immediate and
useful: introduce a coupling whose energy depends on pair overlaps while
keeping the one-block histograms fixed; `Lambda_f` then no longer determines
future response.
