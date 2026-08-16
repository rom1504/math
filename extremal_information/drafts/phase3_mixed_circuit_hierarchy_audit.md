# Independent audit: mixed-circuit holonomy hierarchy

**Status.**  Fresh adversarial reconstruction of
`phase3_mixed_circuit_hierarchy.md`.  The mathematical cores of MC.1--MC.6
survive.  In particular, the gluing torsor has exactly `2^(D kappa)` points,
the nullity bound is valid and sharp for the antipodal endpoint, and the
one-channel packing has the claimed exponential size.  No counterexample
was found.  Three scope qualifications at the end should accompany surface
promotion.

The originating verifier checks 24,551 exhaustive lifted families, 2,000
random fragment partitions, the arity construction through seven sources,
the sharp nullity family, and a 256-point response packing.  The arguments
below do not rely on that code.

## 1. Elimination and circuits

For a selected indexed set `R`, the quotient coordinate is fixed as `q_R`
and the unique systematic-basis correction to target `(u,q_0)` has weight
`|u+a_R|`.  This proves MC.9 directly.  At `u=t,q_0=0`, complementation in
the fixed basis gives

```math
|t+a_R|=D-|a_R|,
```

so minimization is exactly `D-max_R(|a_R|-|R|)`.  MC.1 has no hidden
spanning or simplicity assumption.

Every binary matroid cycle is a disjoint union of circuits: remove a
minimal nonempty subcycle and iterate on the complementary cycle.  For such
a partition, Hamming subadditivity bounds total excess by the sum of the
positive circuit excesses.  Disjoint nonempty circuit indicators are
linearly independent, so their number is at most the cycle-space nullity.
This proves all inequalities in MC.13.  Replacing each positive excess by
`eta(C)` proves MC.14.  The parallel-pair construction has `nu` disjoint
circuits with disjoint kernel holonomies, and therefore attains `nu eta`;
the nullity factor is not a proof artifact.

## 2. Complete gauge invariant

Write

```math
q:\mathbb F_2^E\longrightarrow U,
\qquad
a:\mathbb F_2^E\longrightarrow W.
```

If two assignments agree on `Z=ker q`, their difference vanishes on the
kernel of `q` and hence factors uniquely through `U`.  Conversely every
factor `Lq` vanishes on `Z`.  Thus restriction `a|Z` is exactly the complete
invariant under kernel-fixing linear shears.  MC.3 is correct for labeled,
indexed quotient columns, including loops and parallel columns.

## 3. The MC.4 torsor and exact cardinality

For fragment `j`, let `U_j=im q_j`.  After fixing one representative of
each local gauge class, every other compatible assignment is obtained by a
tuple

```math
(L_1,\ldots,L_m)\in\bigoplus_j\operatorname{Hom}(U_j,W).
```

Two such tuples are globally gauge-equivalent precisely when they are the
restrictions of one `L in Hom(U,W)`, where `U=sum_j U_j`.  Therefore the
space of gluing classes is the cokernel of the injective restriction map

```math
0\longrightarrow\operatorname{Hom}(U,W)
\longrightarrow\bigoplus_j\operatorname{Hom}(U_j,W).          \tag{MA.1}
```

On the primal side there is an exact sequence

```math
0\longrightarrow Z/Z_{\rm loc}
\longrightarrow\bigoplus_jU_j
\overset{\sum}{\longrightarrow}U
\longrightarrow0.                                               \tag{MA.2}
```

Every linear map from the left-hand kernel to `W` extends to the middle
vector space.  Restriction identifies the cokernel in (MA.1) with
`Hom(Z/Z_loc,W)`.  Its dimension is

```math
D\left(\sum_j\dim U_j-\dim U\right)=D\kappa,
```

so the number of points is exactly `2^(D kappa)`.  For two fragments,
dimension inclusion--exclusion gives `kappa=dim(U_1 cap U_2)`.

This count is an exact count of **labeled gauge gluing classes**, not a
lower bound on bits required for one scalar response.  The originating note
states this qualification and should retain it.

## 4. Arity hierarchy

The quotient columns in MC.23 form one circuit and have no proper
dependence.  Every proper offset assignment is therefore removable by a
linear shear defined on an independent family.  In the full family, the
sole nonempty cycle has holonomy `t`, cardinality `r`, and excess `D-r`.
MC.1 gives length `r`.  Hence MC.5 proves a genuine strict arity hierarchy;
it is not merely failure of a particular local certificate.

## 5. MC.6 response packing

For the two parallel lifts `(0,q),(v,q)`, quotient-zero selections are
empty or the two-column circuit, giving

```math
F_v(u)=\min\{|u|,2+|u+v|\}.
```

At least half of the binary cube has weight at least `D/2`.  Greedy packing
inside that set while deleting Hamming balls of radius `<D/4` leaves

```math
2^{(1-H_2(1/4)-o(1))D}
```

centres at mutual distance `D/4-O(1)`.  At query `u=v`, the own profile is
two, while every other packed profile is at least
`min(D/2,2+D/4-O(1))`.  Thus their uniform separation is
`D/4-O(1)`, and decoder error below `D/8` requires distinct states.  MC.6
and its `Omega(D)` information conclusion are correct.

## 6. Relation to the affine-circuit-rank quotient

There is no collision with TC.5 of
`phase3_transversal_composition_growth.md` once scopes are separated.

* MC.4 counts exact labeled mixed holonomies.  For complete linear graph
  fragments, a new labeled map direction carries `D dim(Q)` raw holonomy
  bits.
* TC.5 forgets the labeled list and retains only its affine subspace in map
  space.  That feature state updates exactly but decodes responses only to
  `O(r)` error, where `r` is affine rank.
* MC.6 proves that one independent mixed channel can carry `Omega(D)` rooted
  response bits.  It does not contradict TC.5, whose rank-one interval may
  have constant width while the **oriented affine line itself** still
  carries `Theta(D)` bits.  In the `dim Q=1` instance, the TC.5 state is
  precisely the line `span{v}` (and over `F_2` its nonzero point recovers
  `v`).  Low affine rank bounds decoder distortion, not the number of
  possible oriented rank-`r` states.

The phrase “structured encoding of mixed holonomies” is safest as
“structured **quotient** of mixed holonomies”: the affine state intentionally
forgets which labeled maps occurred.

## 7. Scope corrections

No theorem needs mathematical repair, but three phrases should remain
narrow.

1. MC.13--MC.14 control `Delta` and hence the single antipodal response
   `ell(t)`.  They are not a uniform word-profile or all-future replacement
   theorem.  “Approximate replacement” should be read as “approximate
   antipodal defect bound.”
2. MC.6 is a packing for the declared family of rooted kernel-endpoint
   queries.  It is not by itself a packing in the appended-fragment
   covering-radius response metric.
3. `2^(D kappa)` is the exact number of compatible labeled **gauge classes**.
   A response quotient may collapse many of them; MC.6 supplies a linear
   lower bound only on its one-channel subfamily, not `D kappa` bits in full
   generality.

With these qualifications, the hierarchy is rigorous and complementary to
the affine-rank response law.  Its main new conclusion is that composition
creates gauge-invariant mixed holonomy on the quotient cycle space, and
cycle-space nullity is the exact amplification parameter for the antipodal
defect.
