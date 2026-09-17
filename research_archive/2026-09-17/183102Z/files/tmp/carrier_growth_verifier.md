# Independent audit: compositional carrier--exposure law

**Scope.** I audited `/home/math/quadra/tmp/carrier_growth_unification.md`,
with particular attention to Theorems 3.1 and 4.1, projective baselines,
external covers, irrational generators, nonnegative histogram constraints,
and the finite-grid mean-field constants.

## Verdict

The two main mathematical arguments are sound. In particular:

* the exact arithmetic exponent really is
  `r_Z=rank_Z <phi_j-phi_d>` and can exceed ordinary real query dimension;
* the lower polynomial packing respects histogram nonnegativity;
* the upper metric cover is legitimately **external**, so its real and
  possibly negative last coordinate is not an error;
* irrational atom responses do not invalidate the exact exponent--they are
  precisely why arithmetic rank is needed--while `alpha=0` or `sigma=0`
  correctly prevents an unjustified robust lower bound;
* for the literal, anchored mean-field response norm,
  `alpha>=Delta/4` and `sigma=Delta` are correct; and
* the covering constants in Theorem 3.1 are conservative but valid.

I found two formal corrections and one scope clarification. None destroys
the theorem.

## Corrections

### 1. Theorem 4.1 must say explicitly that (4.3b)--(4.3c) are on a fixed-mass slice

The map `T` in (4.2) is defined only on

```math
V={z:sum_j z_j=0}.
```

Thus `T(c-c')` is defined as written only when `c,c'` have the same total
mass. The proof immediately uses fixed mass and (4.3c) counts
`H_(n,d)`, so the intended result is correct, but the theorem's first
sentence currently reads as an unrestricted contextual-equivalence claim.
It should read:

> On every fixed-mass slice `H_(n,d)`, contextual equivalence is exactly
> `c~c' iff T(c-c')=0`, and ...

Across different masses one must either retain mass, extend `T` to all of
`R^d`, or specify that the future experiment quotients the corresponding
baseline. The fixed-mass exponent is unaffected.

### 2. The blanket projective statement does not preserve the displayed constants

Section 2 says that everything applies projectively after replacing the sup
norm by `osc/2`, and that a fixed anchor converts the projective version to
the literal one. The structural results do transfer, but anchoring is not an
isometry: for `g(q_0)=0`,

```math
||[g]||_(osc/2) <= ||g||_infty <= 2||[g]||_(osc/2).
```

Consequently `R`, `alpha`, `sigma`, and covering radii can change by a
factor two unless the carrier/readouts are formulated intrinsically in the
quotient norm. Individual readouts `D_q` are also not invariant under
adding a constant; one should either select an anchored representative and
accept the factor, or use difference readouts `D_q-D_(q_0)`.

This matters concretely in Section 5. For the **literal anchored** hinge
responses, the report correctly proves

```math
alpha>=Delta/4,\qquad sigma=Delta.
```

If instead responses are quotiented by constants with norm `osc/2`, moving
one site between adjacent bins has response range `[-Delta,0]`, so

```math
sigma_proj=Delta/2,
```

not `Delta`. A conservative repetition of (5.5) gives
`alpha_proj>=Delta/8` (and it may be sharpened). Section 5 should simply
declare that (5.5)--(5.8) use literal responses; if a projective corollary is
wanted, restate its rescaled constants.

### 3. Clarify the regime behind the “mesoscopic law” sentence

The bounds in (4.5) are correct for every `epsilon>0`. The prose saying they
give `(d-1)log(n/epsilon)` “throughout the mesoscopic range” should specify
the regime (for fixed `Phi`, roughly `alpha lesssim epsilon lesssim n`, up
to constants). Below the lattice scale the entropy saturates at the exact
state count in (4.6), and if `alpha` varies with `n` the constants are not
uniform. This is a wording issue, not a defect in (4.5).

## Detailed checks

### Theorem 3.1

An interval of length `2R` has an external radius-`epsilon` net with at most
`1+2R/epsilon` centres (the displayed bound is looser than necessary).
Taking the Cartesian product and using one-Lipschitz readouts proves (3.1).
The approximate-carrier error adds exactly as stated.

The exposed-cube grid has spacing `4epsilon/alpha` and therefore

```math
(1+floor(alpha a/(2epsilon)))^k
```

points separated in response by at least `4epsilon`. An external
radius-`epsilon` ball cannot contain two such points, so (3.2)--(3.3) are
valid. External rather than internal centres do not weaken this argument.

For the homogeneous binary grammar, all comparison normals lie in
`{-1,0,1}^g` modulo sign, hence at most `(3^g-1)/2` central hyperplanes occur.
Their total face count is `exp(O(g^2))`; the report's coarse
`[4(N_g+1)]^g` bound is safe. On each face the response is one affine image
of dimension at most `g`, and the radius bound gives the stated volumetric
cover. The common witness-independent offset `b(q)` is essential, as the
report notes.

There is a typographical `left\lfloor` in (3.2), but no mathematical error.

### Arithmetic-rank exponent in Theorem 4.1

At fixed mass,

```math
F_c=n phi_d+sum_(j<d)c_j(phi_j-phi_d).
```

A finitely generated subgroup of the real function space (or its quotient
by constants) is torsion-free. Choosing a `Z`-basis bounds all generator
coordinates, so the number of attainable sums is `O_Phi(n^r)`. Conversely,
choose `r` displayed differences independent over `Q`, let each associated
count range from `0` to `floor(n/r)`, put all unused mass in type `d`, and
set all other counts to zero. Every histogram is nonnegative, and rational
independence makes all `(floor(n/r)+1)^r` sums distinct. Thus
`N_n=Theta_Phi(n^r)` including the rank-zero case.

This proof remains valid for irrational atom values. For example, several
rationally independent real scalars can have large arithmetic rank even in
one real query coordinate. In that situation real conditioning can vanish
and `sigma` can be zero by Diophantine approximation, which is exactly why
the report does not infer a robust packing from arithmetic rank alone.

### Metric bounds in Theorem 4.1

The lower packing uses

```math
s=1+floor(2epsilon/alpha),
```

so `alpha s>2epsilon`. Its box of first `d-1` counts has total at most `n`,
and the final count is nonnegative. The upper cover rounds only the first
`d-1` counts. Their total absolute error plus the compensating last-coordinate
error is at most `(d-1)h`, giving response error `epsilon`. The decoded
centre need not be a realizable or nonnegative histogram because the cover
is external; this is explicitly and correctly stated.

If `sigma>0`, distinct integer histograms are more than `2epsilon` apart at
`epsilon<sigma/2`, so even external balls need one centre per point. Under
the theorem's stronger `alpha>0` hypothesis, `sigma>=alpha`; retaining the
separate `sigma` clause is harmless and highlights the possible lattice-only
weakening.

### Mean-field constants

At `lambda=-gamma_k`,

```math
G_z/Delta=S_k=sum_(j>k)z_j(j-k).
```

Tail sums are `S_k-S_(k+1)`, and each coordinate of `z` is a first or second
difference of the `S_k`; hence `||z||_infty<=4 max|S_k|`, proving
`alpha>=Delta/4` in the literal norm. If `z` is a nonzero integer vector,
some integer `S_k` is nonzero, so the norm is at least `Delta`. The adjacent
move `z=e_j-e_(j+1)` has literal response norm exactly `Delta`, proving
`sigma=Delta`. Substitution into (4.5) gives the constants in (5.7), and
(5.8) follows.

## Bottom line

The carrier/exposure theorem and additive arithmetic-rank law pass audit.
Before repository integration, restrict Theorem 4.1 explicitly to fixed
mass and either scope Section 5 to literal anchored responses or rescale its
projective constants. The nonnegative-simplex, irrational-precision, and
external-centre concerns do not yield counterexamples.
