# Boolean synchronization in an active eigenspace

**Status.** Rigorous task-local information lower bound and matching-order
finite recovery condition.  A uniformly small Boolean-versus-spherical
trust gap over a `d`-dimensional active eigenspace requires a sphere-covering
number of near-eigen Boolean witnesses.  Conversely, a finite near-eigen
support cover is sufficient, and any pointwise recovery rule can be thinned
to a library of the same covering order.

This isolates a checkable intermediate state.  It is stronger than low
moments or one aggregate Gram matrix, but can be far smaller than the full
Boolean energy landscape when the active dimension is small.

## 1. Active trust queries

Let `H` be real symmetric with

```math
||H||_(2->2)<=r,
```

and let `V` be a `d`-dimensional subspace of the `+r` eigenspace.  Write
`P=P_V`.  For a Boolean witness `x`, define

```math
q(x)={x^THx\over rn}<=1,
\qquad
u(x)={Px\over\sqrt n}\in V,
\qquad ||u(x)||<=1.                                \tag{BS.1}
```

At field strength `beta>0`, the normalized spherical and Boolean responses
in direction `v in S(V)` are

```math
\begin{aligned}
S_beta(v)
&=\max_{||z||=\sqrt n}
 \left\{{z^THz\over2rn}+beta{v^Tz\over\sqrt n}\right\}
 ={1\over2}+beta,\\
B_beta(v)
&=\max_{x\in\{+-1\}^n}
 \left\{{q(x)\over2}+beta\langle v,u(x)\rangle\right\}.
                                                               \tag{BS.2}
\end{aligned}
```

The spherical identity is exact: the spectral and Cauchy bounds give the
upper value, and `z=sqrt(n)v` attains it.

For `delta>=0`, let

```math
\mathcal W_delta=\{x\in\{+-1\}^n:1-q(x)<=delta\}.  \tag{BS.3}
```

## 2. A support-covering lemma

### Lemma BS.1 (support deficit costs spherical entropy)

Let `U` be a finite subset of the Euclidean unit ball of `R^d`.  If

```math
\min_{v\in S^{d-1}}\max_{u\in U}\langle v,u\rangle
\ge1-epsilon,
\qquad 0<epsilon<=1/2,                              \tag{BS.4}
```

then, for `d>=2`,

```math
|U|\ge(2epsilon)^{-(d-1)/2}.                       \tag{BS.5}
```

Thus `|U|>=(c/sqrt(epsilon))^(d-1)` with an absolute `c`.

#### Proof

Whenever `u` serves a direction in BS.4, its inner product is positive.
Normalizing `u` can only increase that inner product.  The normalized
centres therefore cover the sphere by angular caps of radius

```math
theta=arccos(1-epsilon),
\qquad
sin(theta)=\sqrt{2epsilon-epsilon^2}<=\sqrt{2epsilon}. \tag{BS.6}
```

Because `epsilon<=1/2`, `theta<=pi/3`.  Orthogonal projection of one cap
onto its tangent `(d-1)`-ball has radius `sin(theta)` and surface Jacobian at
most `1/cos(theta)<=2`.  Hence its surface fraction is at most

```math
{2v_(d-1)\over d v_d}\sin(theta)^(d-1)
<=\sin(theta)^(d-1),                               \tag{BS.7}
```

where `v_k` is the unit `k`-ball volume; the displayed volume ratio is at
most one for `d>=2`.  Covering total surface measure one proves BS.5.
`square`

The exponent `d-1` is the intrinsic sphere dimension, and the square root
comes from angular deficit `1-cos(theta)=Theta(theta^2)`.

## 3. Trust-gap necessity and recovery sufficiency

Define the uniform active trust gap

```math
Gamma_beta=\sup_{v\in S(V)}(S_beta(v)-B_beta(v)).   \tag{BS.8}
```

### Theorem BS.2 (Boolean synchronization law)

The following two implications hold.

1. **Necessity.**  If `Gamma_beta<=eta`, then for every `v in S(V)` there is
   `x in mathcal W_(2eta)` with

   ```math
   \langle v,u(x)\rangle\ge1-{eta\over beta}.       \tag{BS.9}
   ```

   Consequently, when `0<eta/beta<=1/2`,

   ```math
   |\mathcal W_(2eta)|
   \ge\left({beta\over2eta}\right)^{(d-1)/2}.       \tag{BS.10}
   ```

2. **Sufficiency.**  If a finite library `mathcal L subset mathcal W_delta`
   satisfies

   ```math
   \min_{v\in S(V)}\max_{x\in\mathcal L}
       \langle v,u(x)\rangle\ge1-epsilon,           \tag{BS.11}
   ```

   then

   ```math
   Gamma_beta<=delta/2+beta epsilon.                \tag{BS.12}
   ```

#### Proof

For fixed `v`, take a Boolean optimizer in BS.2.  The exact gap decomposes
as

```math
S_beta(v)-B_beta(v)
={1-q(x)\over2}
 +beta(1-\langle v,u(x)\rangle).                   \tag{BS.13}
```

Both terms are nonnegative.  A gap at most `eta` therefore gives
`1-q(x)<=2eta` and BS.9.  Lemma BS.1 gives BS.10.

Conversely, for each `v` evaluate the Boolean response at the witness from
BS.11.  Its deficit in BS.13 is at most `delta/2+beta epsilon`, uniformly in
`v`. `square`

For an explicit witness-library carrier, BS.10 says that a uniformly
accurate active response requires at least

```math
{d-1\over2}\log_2(beta/eta)-O(d)                  \tag{BS.14}
```

bits merely to index the necessary near-eigen witnesses.

## 4. A matching-order finite extraction theorem

The sufficiency condition does not require retaining every near-eigen
Boolean state.

### Theorem BS.3 (pointwise recovery thins to a finite carrier)

Assume `d>=2`.  Suppose that for every `v in S(V)` there exists
`x(v) in mathcal W_delta`
with

```math
\langle v,u(x(v))\rangle\ge1-epsilon,
\qquad 0<epsilon<=1/4.                              \tag{BS.15}
```

Then there is a sublibrary `mathcal L` of these witnesses with

```math
|\mathcal L|\le(C/\sqrt epsilon)^(d-1)             \tag{BS.16}
```

and

```math
\min_{v\in S(V)}\max_{x\in\mathcal L}
 \langle v,u(x)\rangle\ge1-4epsilon.               \tag{BS.17}
```

Consequently this finite library certifies

```math
Gamma_beta<=delta/2+4beta epsilon.                 \tag{BS.18}
```

#### Proof

Take a Euclidean `sqrt(epsilon)`-net `{v_j}` of `S(V)` of cardinality at
most `(C/sqrt(epsilon))^(d-1)`; this is the standard spherical volumetric
upper bound.  Retain one witness `x(v_j)` for every net point and write
`u_j=u(x(v_j))`.

Equation BS.15 and `||u_j||<=1` give

```math
||u_j-v_j||^2
=||u_j||^2+1-2\langle v_j,u_j\rangle
<=2epsilon,                                        \tag{BS.19}
```

while `||u_j||>=1-epsilon` gives

```math
1-||u_j||^2<=2epsilon.                              \tag{BS.20}
```

For arbitrary `v`, choose `v_j` within `sqrt(epsilon)`.  Then

```math
||v-u_j||<=(1+\sqrt2)\sqrt epsilon.                \tag{BS.21}
```

Using

```math
1-\langle v,u_j\rangle
={1-||u_j||^2+||v-u_j||^2\over2},                  \tag{BS.22}
```

the right side is less than `4epsilon`.  Theorem BS.2 gives BS.18.
`square`

Thus the sphere-covering lower and upper library orders match up to absolute
constants and a constant factor in distortion.  The unproved part in any
specific model is the **Boolean realization** hypothesis BS.15, not the
compression of a realized family.

## 5. A checkable angular recovery certificate

Suppose in addition that the rest of the spectrum lies below the active edge:

```math
H|_(V^perp)<=r(1-gamma)I,
\qquad gamma>0.                                    \tag{BS.23}
```

Then every `x in mathcal W_delta` obeys

```math
1-||u(x)||^2<=delta/gamma.                          \tag{BS.24}
```

Indeed

```math
q(x)<=||u(x)||^2+(1-gamma)(1-||u(x)||^2).
```

Assume `delta<=gamma` and `0<=alpha<=1`.  If the nonzero normalized
projected directions

```math
\widehat u(x)=u(x)/||u(x)||
```

form an angular support cover

```math
\min_v\max_x\langle v,\widehat u(x)\rangle
\ge1-alpha,                                        \tag{BS.25}
```

then

```math
\langle v,u(x)\rangle
\ge\sqrt{1-delta/gamma}(1-alpha)
\ge1-alpha-delta/gamma.                            \tag{BS.26}
```

Hence a finite list of witness defects plus a spherical covering-radius
calculation certifies

```math
Gamma_beta
<=delta/2+beta(alpha+delta/gamma).                 \tag{BS.27}
```

Equivalently, BS.11 can be checked as the convex-geometric inradius condition

```math
(1-epsilon)B_V
\subseteq conv\{u(x):x\in\mathcal L\}.             \tag{BS.28}
```

after symmetrizing the library by `x -> -x`.  This is a finite support-body
certificate, not storage of all `2^n` energies.

## 6. Exact-eigen witnesses cannot give arbitrarily fine recovery

There is a useful rigidity check.  If `V` is the complete strict top
eigenspace, every zero-defect Boolean witness lies in `sqrt(n)V`.  A
`d`-dimensional subspace contains at most `2^d` Boolean cube vertices: choose
`d` coordinate functionals whose restrictions to `V` are independent; their
Boolean values inject the intersection into `{+-1}^d`.

Combining this fact with Lemma BS.1 shows that for `d>=2`, an exact-eigen
library satisfying BS.4 must have

```math
epsilon>=2^{-1-2d/(d-1)}>=1/32.                   \tag{BS.29}
```

Thus arbitrarily small uniform trust gap in a genuinely multidimensional
active space cannot come solely from exact Boolean eigenvectors.  It must use
a shell of near-eigen witnesses, or else the active dimension must collapse.

## 7. Scope and consequence

The theorem applies to one signed trust channel.  An absolute quadratic cap
has two channels; apply the same test to the dangerous active eigenspace of
each.  The lower bound declares a continuum of unit directions at one fixed
field strength.  An actual `p`-port SA.3 table exposes at most `2^p`
directions, with direction-dependent strength
`beta_epsilon=m||z_epsilon||/(r sqrt(n))`; the positive and negative
eigenspaces must be handled separately.  Thus BS.10 is not automatically an
unconditional bit lower bound for every algebraically generated finite
carrier.  It does not assert that the required Boolean synchronization
occurs for arbitrary sign matrices.  Instead it provides:

1. a sharp information price if synchronization does occur;
2. a finite, checkable recovery condition sufficient for uniform small
   spherical integrality gap;
3. an obstruction showing why a large active eigenspace cannot be represented
   by a few planted Boolean poles.

The relevant state is a near-eigen support body or an optimal-order witness
net.  It lies strictly between one global spectral summary and the entire
Boolean landscape.
