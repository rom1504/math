# Independent audit of the partial-common-row obstruction

Status: **passed**.  This note independently checks
[`generic_partial_common_row_shell_lifetime_no_go.md`](generic_partial_common_row_shell_lifetime_no_go.md),
with emphasis on the added diffuse coordinates, the reverse product
projection, and posterior invariance.  The construction is correctly scoped
as a generic rank-one channel and not an actual-child sequence.

## 1. Entropy, collision, and cap spread

Put `ell=m-k`.  The latent map

```math
(\sigma,\xi_1,\ldots,\xi_\ell)
\longmapsto Q=(\sigma\mathbf1_k,\xi)\mathbf1_n^{\mathsf T}
```

is injective.  Its uniform support therefore has exactly `2^(ell+1)` words,
so both the maximum atom and collision probability equal `2^(-(ell+1))`.
There is no projective factor missing in PC.4: the two antipodal signed words
are distinct members of this support.

For any centre `uv^T`, the rank-one correlation factors as in PC.5.  If it
is at least `1-delta` in absolute value, then

```math
|\langle X,u\rangle|\ge(1-\delta)m.
```

The common coordinates contribute at most `k` in absolute value, hence the
sum of the `ell` independent coordinates has absolute value at least
`ell-delta m`.  For `delta<ell/m`, this tail event has probability

```math
{2\over2^\ell}
\sum_{j\le\lfloor\delta m/2\rfloor}{\ell\choose j}.  \tag{VPC.1}
```

Thus PC.6 is uniform over all rank-one centres.  If `ell/m->a` and
`delta<a`, the normalized radius is `delta/(2a)<1/2`, and the binomial
entropy estimate gives exactly the negative exponent in PC.7.  The bound is
not claimed to be the exact finite cap mass, but its exponential rate and
its stated range are correct.

## 2. Likelihood and inverse-escort factorization

Writing `S_i=sum_j B_ij`, direct averaging over the shared sign and the
independent tail signs gives

```math
p(B)
={\cosh(u\sum_{i\le k}S_i)\over(\cosh u)^{kn}}
 \prod_{i>k}{\cosh(uS_i)\over(\cosh u)^n}.           \tag{VPC.2}
```

Raising this density to `-lambda` and normalizing preserves the product
between the common block and every tail row.  Hence PC.9--PC.13 have the
correct constants and normalizers.

There is no hidden optimal-product flaw.  For an arbitrary full row product
`P=prod_i P_i`, the reference itself factors as
`q=q_C tensor prod_(i>k)q_i`, so

```math
D(P\Vert q)
=D\!\left(\bigotimes_{i\le k}P_i\middle\Vert q_C\right)
 +\sum_{i>k}D(P_i\Vert q_i).                         \tag{VPC.3}
```

The tail infimum is exactly zero, attained at `P_i=q_i`; therefore
`I^leftarrow(q)=I^leftarrow(q_C)`.  Likewise the canonical product uses the
same `q_i` on every tail row, giving `J(q,r)=J(q_C,r_C)`.  In particular,
the product used for the `O(sqrt N)` upper bound is fair only on the common
block and equals `q_i` on the tail.  Using a fully fair product would incur
an unnecessary linear tail cost, but the source does not make that mistake.

## 3. Asymptotics

Since `n/N->nu`, under one fair row

```math
uS_i\Longrightarrow N(0,\beta^2\nu).
```

The one-row inverse tilt therefore has the positive KL limit `d_0` in
PC.16.  There are `k=kappa N+o(N)` common rows.  The remaining global
common-block terms are `O(sqrt N)` by the same variance and central-window
bounds as Theorem 37.58, so

```math
J=\kappa d_0N+o(N).                                  \tag{VPC.4}
```

For the comparison product described above,

```math
D(U_{kn}\Vert q_C)
=\lambda E_U\log\cosh(uS_C)
 +\log E_U\cosh(uS_C)^{-\lambda}
\le\lambda u\sqrt{kn}=O(\sqrt N).                   \tag{VPC.5}
```

The sign of the normalizing term is correct because the integrand is at
most one.  This verifies PC.16--PC.19.

## 4. Averaged posterior invariance

The row-switching group `G` acts transitively on the latent support.  It
preserves `U` and the mixture likelihood, while

```math
k_X(T_gB)=k_{Xg}(B).
```

Every bridge law whose density is a function of `p(B)` is consequently
`G`-invariant.  Changing variables by `T_g` in

```math
\bar\mu(X)=\mu(X)E_w[k_X(B)/p(B)]
```

shows that the bracketed multiplier is constant on the entire support.
Normalization makes it one, proving `bar mu=mu`.  Pushforward and
conditional KL terms for every deterministic quotient therefore vanish,
as stated in PC.22--PC.23.

## 5. Verdict

All requested claims pass.  The example simultaneously has exponential
global atom/collision spread, exponentially small fixed-radius rank-one
caps, zero averaged-posterior retuning, and
`J-I^leftarrow=Theta(N)`.  Its limitation is also correctly identified: the
fixed right factor violates factorwise actual-child spread.  It therefore
rules out a generic shell-to-lifetime inequality but does not falsify a
future theorem using the full optimizer-specific factor structure.
