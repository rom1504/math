# Robust entropy of binary max-affine presentations

## Status and scope

This draft records a rigorous lower bound for the shared-parameter response
class in Theorem 16.8.  It constructs a **binary max-affine presentation
with tie-free generic parameter witnesses**

```math
f_\theta(x)=\max_{a\in A_x}\langle a,\theta\rangle,
\qquad A_x\subseteq\{0,1\}^m,
```

with arbitrary real parameters.  It is **not** presently a unit-boundary-load
pure-Max-Cut construction.

The result shows that fixed-scale response entropy can be
`Omega(m log m)`.  It does not show that the current `O(m^2)` upper bound is
sharp.  The unresolved gap is isolated at the end.

## Response metric

For a nonempty `A subset {0,1}^m`, put

```math
h_A(\theta)=\max_{a\in A}\langle a,\theta\rangle.
```

For functions on a finite query set `X`, use the projective shape metric

```math
d_{\rm sh}([f],[g])
={1\over2}\operatorname {osc}_{x\in X}(f(x)-g(x)).       \tag{1}
```

It quotients out a common additive response constant.

## Facet-deletion exposure theorem

### Theorem

Let `P=conv(V) subset R^m` be a full-dimensional `0/1` polytope, with
`V subset {0,1}^m`, and let `Fac(P)` be its facets.  There is a binary
max-affine presentation with `m` shared parameters and
`1+|Fac(P)|` queries such that it contains response shapes

```math
\{[f_{\theta_F}]:F\in\operatorname {Fac}(P)\}
```

satisfying

```math
\|[f_{\theta_F}]\|_{\rm sh}={1\over2},
\qquad
d_{\rm sh}([f_{\theta_F}],[f_{\theta_G}])=1
\quad(F\ne G).                                      \tag{2}
```

For every `eta>0`, the parameters can instead be chosen away from every
optimizer-tie hyperplane, positively rescaled to radius exactly `1/2`, and
chosen with pairwise distance greater than `1-2eta`.

### Proof

Use one base query with witness set

```math
A_0=V.
```

For every facet `F`, write `V_F=V cap F` and use the deletion query

```math
A_F=V\setminus V_F.                                \tag{3}
```

This set is nonempty because a facet of a full-dimensional polytope is a
proper face.

Choose an outward normal `u_F` exposing `F`, and define

```math
\begin{aligned}
\beta_F&=h_V(u_F),\\
\gamma_F
&=\beta_F-\max_{v\in V\setminus F}\langle v,u_F\rangle>0,\\
\theta_F&={u_F\over\gamma_F}.
\end{aligned}                                      \tag{4}
```

Write `c_F=beta_F/gamma_F`.  At `theta_F`, the base response is `c_F` and
the `F`-deletion response is `c_F-1`.  If `G\ne F` is another facet, then
`F` is not contained in `G`: two maximal proper faces related by inclusion
must coincide.  Hence a vertex in `F\setminus G` remains in `A_G`, so the
`G`-deletion response is still `c_F`.  Therefore, in the coordinates indexed
by the base query and all facets,

```math
f_{\theta_F}=c_F\boldsymbol 1-e_F.                 \tag{5}
```

Now `osc(e_F)=1`, whereas `osc(e_F-e_G)=2` for distinct facets.  Equations
(1) and (5) prove (2).

All the support functions involved are continuous and there are only
finitely many optimizer-tie hyperplanes.  A generic parameter can therefore
be chosen arbitrarily close to each `theta_F`.  The corresponding response
shape is arbitrarily close to `[-e_F]`.  Positive rescaling preserves
tie-freeness and restores shape norm exactly `1/2`, proving the robust
generic version.  This gives response-metric robustness, not a uniform
parameter-space margin from the tie hyperplanes. `square`

The same argument works for any antichain of nonempty proper faces: use
their deletion queries and relative-interior exposing normals.  In
particular, all faces of any fixed dimension form such an antichain.

## Macroscopic entropy consequence

Gatzouras, Giannopoulos and Markoulakis proved that an absolute constant
`c>0` exists such that, for all sufficiently large `m`, some full-dimensional
`0/1` polytope in `R^m` has at least

```math
\left({c m\over(\log m)^2}\right)^{m/2}            \tag{6}
```

facets.  Applying the theorem gives a non-strictly `1`-separated family
inside the fixed shape ball of radius `1/2` (equivalently, a `rho`-packing
for every `rho<1`).  Thus, for every fixed `delta<1/2`, even with external
cover centres,

```math
\log \operatorname {Cov}^{\rm ext}_\delta
\ge {m\over2}
   \bigl(\log m-2\log\log m-O(1)\bigr)
=\Omega(m\log m).                                  \tag{7}
```

This is semantic response separation, not a parameter-precision count: the
possibly small facet gap is removed by the normalization in (4).

The primary lower-bound source is:

* D. Gatzouras, A. Giannopoulos and N. Markoulakis,
  [*Lower bound for the maximal number of facets of a 0/1 polytope*](https://doi.org/10.1007/s00454-005-1159-1),
  *Discrete & Computational Geometry* **34** (2005), 331--349;
  [arXiv:math/0406125](https://arxiv.org/abs/math/0406125).

For comparison, Fleiner, Kaibel and Rote proved that every
`m`-dimensional `0/1` polytope has `O((m-2)!)` facets.  Their bounded-coordinate
face estimates also control the whole face lattice at `exp(O(m log m))`
scale.  Consequently this **single-polytope, one-shape-per-face deletion**
mechanism has response packing exponent
`Theta(m log m)`, up to constants and lower-order logarithms; it cannot
produce `exp(Omega(m^2))` robust shapes.

The primary upper-bound source is:

* T. Fleiner, V. Kaibel and G. Rote,
  [*Upper bounds on the maximal number of facets of 0/1-polytopes*](https://doi.org/10.1006/eujc.1999.0326),
  *European Journal of Combinatorics* **21** (2000), 121--130.

## Raw chambers versus robust response shapes

For clarity, let `C_m(R,delta)` denote the supremum, over finite query sets
and binary witness families, of

```math
\log\operatorname {Cov}_\delta
\{[f_\theta]:\|[f_\theta]\|_{\rm sh}\le R\}.
```

The radius restriction is essential: homogeneity makes the unrestricted
parameter image unbounded, so its fixed-radius covering number is infinite.

The elementary arrangement bound allows `exp(O(m^2))` relatively open cells
for all ternary comparison normals.  The theorem above proves only

```math
\Omega(m\log m)
\le C_m(1,\delta)
\le O(m^2)                                         \tag{8}
```

for fixed `delta<1/2`, among radius-bounded binary shared-parameter response
classes.  It gives two concrete cautions about interpreting the upper bound:

1. A change of optimizer may occur across a comparison hyperplane whose
   support gap is negligible compared with the projective response radius.
   Counting that chamber as a macroscopically new response can overcount.
2. A single `0/1` support polytope cannot witness the full ternary chamber
   count through exposed faces: its entire face lattice has only
   `exp(O(m log m))` members.

An arbitrary presentation uses many support polytopes
`P_x=conv(A_x)`.  The available arrangement estimate permits as many as
`exp(O(m^2))` common cells, and the single-polytope face bounds do not
control that common refinement.  Therefore (7) is not an `O(m log m)` upper
bound for general presentations.

One sharp pair of routes through the remaining **robust common-fan** problem
is:

> Either show that every fixed-scale separated response family can be
> charged, with bounded multiplicity, to only `exp(O(m log m))` exposed
> bounded-coordinate faces, or construct a multi-polytope family with
> `exp(Omega(m^2))` projectively normalized response shapes separated by a
> fixed constant.

Raw chamber cardinality alone resolves neither side.
Intermediate exponents are also possible; the two displayed routes are
useful decisive targets, not an asserted exhaustive dichotomy.
