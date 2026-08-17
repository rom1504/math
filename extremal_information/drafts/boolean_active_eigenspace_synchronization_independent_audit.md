# Independent audit: Boolean synchronization in an active eigenspace

**Audited files.** `boolean_active_eigenspace_synchronization.md` and
`verify_boolean_active_eigenspace_synchronization.py`.

## Verdict

**PASS WITH TWO THEOREM-DOMAIN REPAIRS AND SCOPE QUALIFICATIONS.**  The cap
constant in BS.1, the necessity/sufficiency decomposition, extraction
constant four, spectral-gap estimate, `2^d` exact-eigen count, and `1/32`
floor are all correct in their intended domains.

Required repairs:

1. BS.3 must assume `d>=2`, or its cardinal bound must be replaced by
   `max{2,(C/sqrt(epsilon))^(d-1)}`.  As written it is false for `d=1`.
2. BS.10 should assume `0<eta/beta<=1/2`; at `eta=0` its displayed quotient
   is undefined.  Zero gap can instead be ruled out separately for `d>=2`.

For the angular certificate, state explicitly that `0<=alpha<=1` and that
only witnesses with nonzero projection are normalized.  These are natural
but used by the inequality in BS.26.

The main scope qualification is that the lower bound prices uniform response
to **every direction of the active sphere at a fixed field strength**.  A
single finite-port SA.3 instance generally exposes only finitely many field
directions, with channel-dependent strengths.  Thus BS.10 is not
automatically an information lower bound for the entire original Boolean
cap or for every compressed response algorithm.

## 1. Spherical normalization and gap decomposition

For `v in S(V)`, the spectral and Cauchy bounds give

```math
{z^THz\over2rn}+\beta{v^Tz\over\sqrt n}
\le {1\over2}+\beta.                                 \tag{ABES.1}
```

The vector `z=sqrt(n)v` lies in the `+r` eigenspace and attains both terms,
so the formula for `S_beta(v)` is exact.

For every Boolean `x`,

```math
q(x)\le1,
\qquad
\langle v,u(x)\rangle
={v^TPx\over\sqrt n}={v^Tx\over\sqrt n}le1.        \tag{ABES.2}
```

If `x` is a Boolean optimizer for direction `v`, direct subtraction gives

```math
S_\beta(v)-B_\beta(v)
={1-q(x)\over2}
 +\beta\{1-\langle v,u(x)\rangle\}.                 \tag{ABES.3}
```

Both summands are nonnegative.  Therefore a gap at most `eta` implies

```math
1-q(x)\le2\eta,
\qquad
\langle v,u(x)\rangle\ge1-\eta/\beta.              \tag{ABES.4}
```

This verifies the necessity constants.  Conversely, a witness with defect
at most `delta` and support deficit at most `epsilon` has gap at most
`delta/2+beta epsilon`, proving sufficiency.  There is no hidden maximizer or
factor-two issue.

The premise `Gamma_beta<=eta` supplies ABES.4 for every sphere direction.
To invoke BS.1 and divide by `eta`, BS.10 should read
`0<eta/beta<=1/2`.  If `eta=0`, ABES.4 would demand an exact unit support in
every direction; the exact-eigen argument in Section 6 shows this is
impossible for finite `n` and `d>=2`.

## 2. Spherical-cap constant in BS.1

Any `u` serving BS.4 has positive inner product.  Replacing it by
`u/||u||` only enlarges that support, so the normalized centres cover
`S^(d-1)` by angular caps with

```math
\cos\theta=1-\epsilon,
\qquad
\sin\theta=\sqrt{2\epsilon-\epsilon^2}
\le\sqrt{2\epsilon}.                                \tag{ABES.5}
```

For `epsilon<=1/2`, one has `theta<=pi/3`.  Orthogonal projection of such a
cap to its tangent `(d-1)`-ball is injective and has surface Jacobian at most
`1/cos(theta)<=2`.  Its fraction of total sphere area is therefore at most

```math
{2v_{d-1}\over d v_d}\sin^{d-1}\theta.              \tag{ABES.6}
```

The prefactor is at most one for every `d>=2`, so one cap has fraction at
most `(2epsilon)^((d-1)/2)`.  Covering measure one requires

```math
|U|\ge(2\epsilon)^{-(d-1)/2}.                       \tag{ABES.7}
```

The constant and exponent in BS.5 are correct.  The independent verifier
checks the volume prefactor and sine inequality in 156 dimension/accuracy
cases.

Applying ABES.7 to the distinct projected vectors from
`W_(2eta)` yields BS.10 because the number of projected vectors is at most
the number of witnesses.  Taking binary logarithms gives BS.14, including
the `-O(d)` term from the factor two.

## 3. Finite extraction and the constant four

At a net centre `v_j`, BS.15 and `||u_j||<=1` imply

```math
\|u_j-v_j\|^2\le2\epsilon,
\qquad
1-\|u_j\|^2\le2\epsilon.                            \tag{ABES.8}
```

For a target `v` within Euclidean distance `sqrt(epsilon)` of `v_j`,

```math
\|v-u_j\|\le(1+\sqrt2)\sqrt\epsilon.                \tag{ABES.9}
```

Using the exact polarization identity,

```math
1-\langle v,u_j\rangle
={1-\|u_j\|^2+\|v-u_j\|^2\over2}
\le {2+(1+\sqrt2)^2\over2}\epsilon
<4\epsilon.                                         \tag{ABES.10}
```

Thus the factor four is valid.  A standard sphere net has cardinality
`(C/sqrt(epsilon))^(d-1)` for `d>=2`, giving BS.16--BS.18.

### Required `d=1` repair

For `d=1`, the sphere is `S^0={-1,+1}`.  Exact pointwise witnesses can serve
both directions, but any library with support at least
`1-4epsilon>0` must contain both signs.  The right side of BS.16 is
`(C/sqrt(epsilon))^0=1`, regardless of `C`.  Hence BS.3 is literally false
in dimension one.  Add `d>=2`, as in BS.1, or include a leading factor/two-
point exception.  The independent verifier contains this explicit
counterexample.

With this repair, lower and upper witness-library orders match at the
sphere-covering scale.  This is a result about explicit witness libraries;
it does not preclude a shorter algebraic generator for those witnesses.

## 4. Spectral-gap/angular certificate

Because `V` is an eigenspace, it reduces `H`.  For a Boolean `x`, decompose
its normalized squared mass into `V` and `V^perp`.  The one-sided spectral
gap gives

```math
q(x)
\le\|u(x)\|^2+(1-\gamma)(1-\|u(x)\|^2),             \tag{ABES.11}
```

and therefore

```math
1-\|u(x)\|^2\le{1-q(x)\over\gamma}\le\delta/\gamma. \tag{ABES.12}
```

This proves BS.24 with no missing lower-spectrum assumption beyond
`||H||<=r`.

Suppose a selected witness has nonzero projection and
`<v,hat u> >=1-alpha`.  For `0<=alpha<=1` and `delta<=gamma`,

```math
\begin{aligned}
\langle v,u\rangle
&=\|u\|\langle v,\widehat u\rangle\\
&\ge\sqrt{1-\delta/\gamma}(1-\alpha)\\
&\ge1-\alpha-\delta/\gamma.                         \tag{ABES.13}
\end{aligned}
```

The positivity of `1-alpha` is used when multiplying the radial lower
bound.  If `alpha>1`, the displayed argument can reverse direction and the
certificate should instead be clipped/trivialized.  At `delta=gamma`, some
near-eigen witnesses may have zero projection; they cannot be normalized
and should simply be excluded from the angular list.  Any witness that
actually serves a nontrivial cover is nonzero.

Substitution into BS.12 proves BS.27.  After symmetrizing the library, the
support-function criterion is equivalent to

```math
(1-\epsilon)B_V\subseteq\operatorname{conv}\{u(x)\}, \tag{ABES.14}
```

so BS.28 is correct.

## 5. Exact-eigen count and the `1/32` floor

If `V` is the complete `+r` eigenspace and `q(x)=1`, equality in the
Rayleigh bound forces `x` to lie in `V`.  A `d`-dimensional subspace contains
at most `2^d` Boolean vertices: choose `d` coordinate restrictions forming
a basis of `V^*`; their `{+-1}^d` values inject `V cap {+-1}^n`.

An exact-eigen support cover with deficit `epsilon<=1/2` therefore obeys

```math
2^d\ge(2\epsilon)^{-(d-1)/2}.                       \tag{ABES.15}
```

Solving gives

```math
\epsilon\ge2^{-1-2d/(d-1)}.                         \tag{ABES.16}
```

The exponent is `3+2/(d-1)`, maximized at `d=2`; hence for every `d>=2`,

```math
2^{-1-2d/(d-1)}\ge2^{-5}=1/32.                      \tag{ABES.17}
```

If `epsilon>1/2`, the `1/32` conclusion is already trivial.  Thus BS.29 is
correct over its full intended range.  The pair-constant examples in the
independent verifier attain exactly `2^d` Boolean vertices, showing the
intersection count itself is sharp.

## 6. Scope relative to the full Boolean response

The abstract theorem declares all directions `v in S(V)` at one fixed
field strength `beta`.  An SA.3 endpoint channel with field
`z_epsilon` maps to this normalization only when its relevant field lies in
the active eigenspace, with

```math
v_\epsilon={z_\epsilon\over\|z_\epsilon\|},
\qquad
\beta_\epsilon={m\|z_\epsilon\|\over r\sqrt n}.      \tag{ABES.18}
```

The actual endpoint family has at most `2^p` such directions and generally
different `beta_epsilon`.  It need not cover the whole active sphere.
Therefore:

- BS.10 gives a lower bound for a declared continuum of active future
  fields, or for finite directions dense enough to inherit the cap bound;
- it does not automatically lower-bound state complexity for one bare
  finite-port maximum;
- BS.11/BS.27 remain valid sufficient certificates when restricted to the
  actual channel directions, with their individual strengths (using the
  worst relevant `beta` in a uniform bound).

For an absolute quadratic response, the positive outer channel uses the
`+r` active eigenspace.  The negative channel must be treated separately by
applying the theorem to `-H` and its `-r` eigenspace.  Fields with substantial
components outside the relevant edge eigenspace require an additional
decomposition; the identity `S_beta=1/2+beta` no longer applies verbatim.

Finally, BS.10 proves that many near-eigen Boolean witnesses must **exist**
and that an explicit witness library needs the corresponding index range.
It is not an unconditional bit lower bound for every possible response
representation: an algebraic generator or a direct support-function formula
could describe many witnesses compactly.  The phrase “information lower
bound” should be read in this witness-library sense unless a communication
or metric-entropy reduction is added.

The theorem is nonetheless a genuine intermediate result.  It identifies
exactly which Boolean support body would close the spherical gap for a rich
active query family, without storing the full `2^n` energy landscape.

## 7. Verifier assessment

The canonical verifier passes.  It checks the circle exponent, exact gap
decomposition on a two-dimensional pair-constant subspace, one extraction
geometry, and a random spectral-gap certificate.  It does not test cap
volume constants above `d=2`, the missing `d=1` extraction case, or the
exact-eigen `1/32` arithmetic.

The independent verifier adds:

- 156 cap-area constant checks through dimension 40;
- 400 exact gap decompositions in dimensions two and three;
- the extraction constant and explicit `d=1` counterexample;
- 1,664 spectral-gap checks, including `gamma>1`;
- the exact-eigen floor through dimension 50.

Run:

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_boolean_active_eigenspace_synchronization.py

./.venv/bin/python \
  extremal_information/experiments/verify_boolean_active_eigenspace_synchronization_independent_audit.py
```

Both pass.
