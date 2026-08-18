# Symmetry orbits are an exact posterior-retuning quotient

Status: **rigorous finite theorem**.  This note extracts the common algebra
behind the product-group and diffuse BSC examples.  A symmetry-invariant
latent channel never needs more than its latent orbit label to record the
*averaged* posterior retuning under a likelihood-dependent disorder tilt.
The statement is exact, but it controls neither the pointwise posterior nor
the canonical row-product gap without an additional directional theorem.

## 1. Equivariant channel

Let a finite group `G` act on finite sets `Z` and `B`.  Let `U` be a
`G`-invariant reference law on `B`, let `mu` be a `G`-invariant prior on
`Z`, and let `k_z` be channel densities relative to `U` satisfying

```math
k_{gz}(gb)=k_z(b)                                   \tag{OQ.1}
```

for every `g,z,b`.  Put

```math
p(b)=\sum_z\mu(z)k_z(b),
\qquad
\mu_b(z)={\mu(z)k_z(b)\over p(b)}.                 \tag{OQ.2}
```

Let `q` be any observation law whose density relative to `U` is a
measurable function of the scalar `p(b)`.  This includes every raw disorder
tilt `dq_s/dU proportional p^s`, positive or negative.  Define the averaged
posterior

```math
\bar\mu_q(z)=E_{b\sim q}\mu_b(z).                  \tag{OQ.3}
```

## 2. Exact quotient theorem

Let `pi:Z -> Z/G` be the orbit map.

**Theorem OQ.1 (orbit posterior quotient).**  Under (OQ.1)--(OQ.3),

```math
\boxed{
D(\bar\mu_q\Vert\mu)
=D(\pi_\#\bar\mu_q\Vert\pi_\#\mu).}              \tag{OQ.4}
```

More strongly, both `mu` and `bar mu_q` are uniform conditional on each
orbit.  Consequently every posterior cross-entropy or KL chain rule can be
evaluated entirely on the orbit alphabet.

*Proof.*  Invariance of `mu,U` and equivariance give

```math
p(gb)=\sum_z\mu(z)k_z(gb)
      =\sum_z\mu(g^{-1}z)k_{g^{-1}z}(b)=p(b).       \tag{OQ.5}
```

Thus `q` is also invariant.  Changing variables `b=gb'` in (OQ.3) gives

```math
\bar\mu_q(gz)
=\mu(gz)\sum_bq(b){k_{gz}(b)\over p(b)}
=\mu(z)\sum_{b'}q(b'){k_z(b')\over p(b')}
=\bar\mu_q(z).                                     \tag{OQ.6}
```

Every invariant probability measure is uniform on each finite orbit.
Writing `O=pi(z)`, therefore,

```math
\mu(z)={\mu_O(O)\over|O|},
\qquad
\bar\mu_q(z)={\bar\mu_O(O)\over|O|}.              \tag{OQ.7}
```

Substitution in the definition of relative entropy proves (OQ.4). `square`

There is an exact approximate-symmetry version.  For an arbitrary law
`nu` on `Z`, let

```math
\mathsf S_G\nu={1\over|G|}\sum_{g\in G}g_\#\nu.     \tag{OQ.8}
```

This symmetrization has the same orbit marginal as `nu` and is uniform
inside each orbit.  Whenever `mu` is `G`-invariant, the KL chain rule gives

```math
\boxed{
D(\nu\Vert\mu)
=D(\pi_\#\nu\Vert\pi_\#\mu)
 +D(\nu\Vert\mathsf S_G\nu).}                       \tag{OQ.9}
```

Thus the exact cost of replacing a posterior by its orbit quotient is not
an unspecified approximation error: it is its symmetry-breaking KL.  The
equivariant likelihood-dependent tilt has zero such cost by (OQ.6).

## 3. Complexity and three applications

If the action has `K_N` latent orbits, (OQ.4) is an exact `K_N`-state
retuning quotient.  It is subexponential precisely when
`log K_N=o(N)`.  This is a genuine information reduction even if individual
orbits are exponentially large.

1. **Uniform subgroup priors.**  In Theorems 37.63 and 37.65 the switching
   group acts transitively on latent support.  Thus `K_N=1` and (OQ.4)
   reduces to `bar mu=mu`.

2. **Diffuse BSC phase.**  In Theorem 37.66, row/column permutations and
   antipodal symmetry have `O(mn)` projective orbits, indexed by the two
   factor magnetization counts.  Equation (OQ.4) is exactly (37.255).

3. **Actual children with symmetry.**  For a child signing `A`, every graph
   automorphism of `A` acts equivariantly on its Gibbs sector and binary
   channel.  For a pair `(A,D)`, the product automorphism group therefore
   gives an exact orbit quotient of the actual latent retuning.  This is
   useful only when the number of joint factor orbits is subexponential;
   an arbitrary minimizing signing may have a trivial automorphism group.

## 4. Exact boundary

The theorem compresses the averaged posterior displacement
`D(bar mu||mu)`.  It does not say that the full inverse escort `q` is a
mixture of row products, nor does it bound `J-I^leftarrow`.  The subgroup
examples demonstrate this separation sharply: their orbit retuning is
zero while their coherent row-product gap is linear.  A successful actual-
child application must therefore combine an orbit (or approximate-orbit)
quotient with a theorem converting its displacement into a coherent product
direction, and must prove that the relevant orbit count is subexponential.
Storing the complete orbit partition when `K_N=e^{Theta(N)}` is not a
compression.

Equation (OQ.9) gives a precise approximate variant: an approximate symmetry
is useful exactly when it has both `log K_N=o(N)` and
`D(bar mu||S_G bar mu)=o(N)`.  Verifying either property for arbitrary
optimizing children remains open.
