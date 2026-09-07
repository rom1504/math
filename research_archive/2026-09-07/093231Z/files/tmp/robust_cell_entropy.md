# Robust cells in binary max-affine response grammars

## Outcome

The raw `O(m^2)` ternary-hyperplane term is **not yet proved necessary** for
arbitrary binary max-affine grammars.  I found, however, an exact robust
construction showing that it cannot in general be reduced below
`Theta(m log m)` in the exponent.  The construction uses facets of a single
`0/1` polytope and has a fixed, dimension-free response separation after
normalization.  For this canonical exposed-face mechanism, classical face
bounds also give the matching `O(m log m)` exponent, so the raw
`Theta(m^2)` chamber count is definitely an overcount there.

This applies to unrestricted Boolean max-affine presentations with `m`
shared real parameters.  It is not currently a unit-load pure-Max-Cut
construction.

## Set-up

For a nonempty `A subset {0,1}^m`, write

```math
h_A(\theta)=\max_{a\in A}\langle a,\theta\rangle.
```

For a finite query set `X`, a binary max-affine response is

```math
f_\theta(x)=h_{A_x}(\theta),\qquad x\in X.
```

We put the projective response metric

```math
d_{\rm sh}([f],[g])={1\over2}\operatorname{osc}_{x\in X}(f(x)-g(x)).
```

Thus constants in all query responses are invisible.

## Facet deletion theorem

**Theorem.**  Let `P=conv(V)` be a full-dimensional `0/1` polytope in
`R^m`, and let `Fac(P)` be its set of facets.  There is a binary max-affine
presentation with `m` shared parameters and `1+|Fac(P)|` queries containing
a family

```math
\{[f_{\theta_F}]:F\in Fac(P)\}
```

such that

```math
\|[f_{\theta_F}]\|_{\rm sh}=\frac12,
\qquad
d_{\rm sh}([f_{\theta_F}],[f_{\theta_G}])=1
\quad(F\ne G).                                      \tag{1}
```

Moreover, for every fixed `eta>0`, the parameters may be chosen off all
optimizer-tie hyperplanes, while replacing the right sides in (1) by at
most `1/2+eta` and at least `1-2eta`, respectively.

**Construction.**  Include a base query `0` with `A_0=V`.  For every facet
`F`, regard `F` also as its set of vertices and include the deletion query

```math
A_F=V\setminus F.
```

This set is nonempty because `P` is full-dimensional and `F` is proper.

Choose an outward exposing normal `u_F` and write

```math
\beta_F=h_V(u_F),
\qquad
\gamma_F=\beta_F-\max_{v\in V\setminus F}\langle v,u_F\rangle>0,
\qquad
\theta_F={u_F\over\gamma_F}.
```

At `theta_F`, the base response is `c_F=beta_F/gamma_F`.  The `F`-deletion
response is `c_F-1`.  If `G` is a distinct facet, then `F` is not contained
in `G`: otherwise the two maximal proper faces would be equal.  Hence some
vertex of `F\setminus G` remains in `A_G`, and the `G`-deletion response is
`c_F`.  Consequently

```math
f_{\theta_F}=c_F\boldsymbol 1-e_F.                 \tag{2}
```

Equation (1) follows immediately: `osc(e_F)=1`, whereas
`osc(e_F-e_G)=2`.  Support functions of finite sets are continuous.  A
generic parameter can be chosen arbitrarily close to each `theta_F`, so
(2) persists to arbitrary prescribed accuracy while all optimizer ties are
removed.  This proves the robust version.

The same proof works for any antichain of nonempty proper faces of `P`, not
only facets: use `A_F=V\setminus F` and a relative-interior exposing normal
for each face.  Equal-dimensional faces form an antichain.

## Quantitative consequence

Gatzouras, Giannopoulos and Markoulakis prove that, for an absolute `c>0`
and all sufficiently large `m`, some full-dimensional `0/1` polytope in
`R^m` has at least

```math
\left({c m\over(\log m)^2}\right)^{m/2}             \tag{3}
```

facets.  Applying the theorem above gives, inside the fixed projective ball
of radius `1/2`, a `1`-packing of that cardinality.  Therefore, for every
fixed `delta<1/2`,

```math
\log \operatorname{Cov}_{\delta}
\ge {m\over2}\bigl(\log m-2\log\log m-O(1)\bigr)
=\Omega(m\log m).                                  \tag{4}
```

This is semantic response separation, not finite-precision parameter
counting: the facet gaps are normalized away individually by
`theta_F=u_F/gamma_F`.

The primary source is:

* D. Gatzouras, A. Giannopoulos, N. Markoulakis, *Lower bound for the
  maximal number of facets of a 0/1 polytope*, Discrete & Computational
  Geometry **34** (2005), 331--349, DOI
  `10.1007/s00454-005-1159-1`, arXiv `math/0406125`.  Its main theorem is
  exactly (3).

For comparison, Fleiner--Kaibel--Rote prove that every `m`-dimensional
`0/1` polytope has at most `C(m-2)!` facets, and more generally only
`exp(O(m log m))` faces in each dimension.  Thus the facet-deletion grammar
has

```math
\log \operatorname{Pack}_{\delta}=\Theta_{\rm exponent}(m\log m)          \tag{5}
```

up to constants and lower-order logarithms; it cannot realize
`exp(Omega(m^2))` robust cells.  The primary upper-bound source is:

* T. Fleiner, V. Kaibel, G. Rote, *Upper bounds on the maximal number of
  facets of 0/1-polytopes*, European Journal of Combinatorics **21**
  (2000), 121--130, DOI `10.1006/eujc.1999.0326`.  Corollary 8 gives
  `O((m-2)!)` facets; Theorem 7 gives
  the corresponding bounded-coordinate face estimates.

## What this says about the `m^2` term

The construction separates three notions that the arrangement proof can
conflate.

1. The arrangement of all ternary comparison normals can have
   `exp(Theta(m^2))` relatively open faces.
2. A single `0/1` support polytope has only `exp(O(m log m))` faces, yet its
   facets already generate `exp(Omega(m log m))` macroscopically separated
   response shapes.
3. Hence optimizer-face counts can be genuinely superexponential while
   still falling exponentially short of the raw common-refinement chamber
   count.

The remaining gap is precise.  An arbitrary presentation uses many support
polytopes `P_x=conv(A_x)`.  Their common normal fan can have
`exp(Theta(m^2))` cells, but a cell change can occur in a query whose support
gap is arbitrarily small compared with the projective response radius.
The hyperplane proof charges that change as a full new cell.  To make its
`m^2` term sharp, one must exhibit `exp(Omega(m^2))` common-fan cells whose
**response vectors**, after projective normalization, remain pairwise a
fixed distance apart.  The facet-deletion construction shows how to turn
one exposed face into one robust coordinate, but the 0/1 face upper bound
prevents doing this `exp(Omega(m^2))` times with one support polytope.

Conversely, (5) is not a proof of an `O(m log m)` bound for arbitrary
multi-polytope max-affine presentations.  Such a proof would need a theorem
showing that fixed-scale projective response separation can be charged to
only `exp(O(m log m))` exposed faces, rather than to the full common normal
fan.  I do not have that theorem, and the classical 0/1-polytope facet bound
alone does not imply it.

## Director-level verdict

There is a rigorous new lower bound and a decisive correction to the naive
route:

```math
\boxed{
\Omega(m\log m)\ \le\ \log \operatorname{Cov}_{\delta}
\ \le\ O(m^2)
}
```

for general binary shared-parameter grammars at fixed response distortion,
with the lower bound realized by a fixed-radius robust family.  The
`O(m^2)` term is therefore not mere finite-precision bookkeeping, but its
necessity is unproved; raw chamber count alone is not evidence for it.  For
the natural exposed-face/deletion subclass the correct exponent is
`Theta(m log m)`.

The clean next theorem is a **robust common-fan charging theorem**: either
show that every fixed-scale separated response family injects, with bounded
multiplicity, into faces of `O(1)` bounded-coordinate polytopes (giving
`O(m log m)`), or construct a multi-polytope family violating exactly that
claim.  This is the point at which the two possible exponents genuinely
diverge.
