# Independent audit: full-range row-magnitude no-gain theorem

**Frozen source:**
`extremal_information/drafts/conference_row_magnitude_full_beta_no_gain.md`

**SHA-256:**
`eea295989cf782ac2e28bbd2904a35e318b243f9d92fc7b706e6c9e00e6abbd4`

**Verdict:** **PASS.**  The projected layer coupling, matrix concentration,
uniform-pressure transfer, convex rank-one restoration, constants, and
one-sided conclusion are correct.  In particular, the endpoint bridge need
not satisfy a high-temperature operator bound: convexity uses the covariance
bound only at the projected base point.  No repair is required.

## Gauge and projected population direction

Right multiplication by `D_v` sends the distinguished vector `v` to the
all-ones vector.  The corresponding spin change conjugates the right child
to `D_v A D_v`, which is again a symmetric conference signing.  Although
the displayed notation continues to use `A` on both diagonal blocks, this
does not change the uniform-pressure center: the audited conference proof
uses only the exact conference power identities, flatness, and symmetric
Bernoulli spectral law, all of which are invariant under switching.  It also
does not change any operator bound.

For `P=11^T/r`, a conditioned row has `S=<R,1>` and

```math
E[S^2\mid E_r]
\le p_r^{-1}E S^2\le r/p_0.
```

Independence of the rows consequently gives

```math
E\|B1\|_2
\le(E\|B1\|_2^2)^{1/2}
=\{rE[S^2\mid E_r]\}^{1/2}
\le r/\sqrt{p_0}.
```

Thus MF.7--MF.8 and the claimed uniformity in the magnitude set are exact.

## Projected layer coupling

Conditional on layer counts `K_i,K_i'`, the nesting coupling makes the
changed set a uniform `d_i=|K_i-K_i'|` subset, with a common change sign.
After projection, its row is therefore, up to sign,

```math
2(1_{T_i}-(d_i/r)1).
```

This verifies both

```math
\|D_i^\circ\|_2^2=4d_i(1-d_i/r)
```

and the exact covariance

```math
E[(D_i^\circ)^TD_i^\circ\mid d_i]
={4d_i(r-d_i)\over r(r-1)}(I-P).
```

The binomial tail before conditioning, division by `p_r>=p_0`, and a union
over `r` rows give

```math
P(\max_i d_i>2r^{3/4})\le e^{-c_{p_0}\sqrt r}
```

after changing constants.  On this event's complement, each independent
PSD row summand has norm `O(r^(3/4))`, and the norm of their summed
conditional expectation is also `O(r^(3/4))` because

```math
\sum_i {4d_i(r-d_i)\over r(r-1)}
\le {4\sum_i d_i\over r-1}=O(r^{3/4}).
```

Matrix Bernstein at threshold `C r^(3/4) log r` has exponent
`Omega(C log r)` after the dimension factor; taking `C` large gives the
stated `O(r^(-10))` failure.  Taking a square root yields

```math
\|D^\circ\|_{op}
=O_P(r^{3/8}\sqrt{\log r})=o_P(\sqrt r).
```

Projection is a Frobenius contraction, so the already audited layer-coupling
estimate also gives `E||D^circ||_F=O(r^(3/4))`.  No independence of entries
within a conditioned row is being assumed.

## Why both projected bridges have the correct pressure center

The iid bridge satisfies `||W||_op<=(2+o_P(1))sqrt(r)`, and
`||W(I-P)||_op<=||W||_op`.  Together with the preceding projected error,

```math
\|W^\circ\|_{op},\ \|B^\circ\|_{op}
\le(2+o_P(1))\sqrt r.
```

For every fixed `beta<sqrt(2)/6`, one can choose `delta>0` and
`kappa<1/2` with `beta(3+delta)/sqrt(2)<kappa`.  Hence both projected
parents lie in the same strict operator ball with probability tending to
one.  The polynomial `O(r^(-10))` projected-coupling failure and exponential
iid norm failure remain negligible after multiplication by the crude
`O(r^(3/2))` pressure cap.

On the common regular event, the archived stability theorem gives exactly

```math
|f(B^\circ)-f(W^\circ)|
\le {K_\kappa\beta\over\sqrt2}\|D^\circ\|_F.
```

Its expected right side is `O(r^(3/4))=o(r)`.  For the iid rank-one
projection,

```math
\|WP\|_*=\|W1\|_2/\sqrt r,
\qquad E\|W1\|_2\le r.
```

Using the nuclear form of the same stability estimate gives expected
pressure cost

```math
K_\kappa {\beta\over\sqrt{2r}}E\|WP\|_*=O(1),
```

so `W^circ` has the same leading pressure as the iid sign bridge.  This
proves MF.22 in probability and `L^1`; the argument does not merely compare
operator norms.

## Convex restoration of the rank-one component

Let `C=BP` and

```math
Y={\beta\over\sqrt{2r}}
\begin{pmatrix}0&C\\C^T&0\end{pmatrix}.
```

The symmetric dilation has nuclear norm `2||C||_*`, while

```math
\|C\|_*=\|B1\|_2/\sqrt r.
```

Therefore the high-temperature differential bound at the regular base
`X(B^circ)` is

```math
|g'(0)|
\le {K_\kappa\over2}\|Y\|_*
={K_\kappa\beta\over\sqrt2\,r}\|B1\|_2.
```

The log-cosh pressure is a log-sum-exp of affine functions of the
interaction, hence is globally convex.  Its supporting line at zero gives

```math
f(B)=g(1)
\ge g(0)+g'(0)
\ge f(B^\circ)
-{K_\kappa\beta\over\sqrt2\,r}\|B1\|_2.
```

This step uses no covariance or operator hypothesis at the full endpoint
`B`.  MF.8 makes the expected *unnormalized* loss `O_{p_0}(1)`, hence its
normalized contribution is `O(1/r)`.  On the exceptional base-irregular
event, `(h_beta-f(B)/r)_+<=h_beta` because pressure is nonnegative, so its
vanishing probability is already sufficient.

Combining this inequality with the `L^1` convergence of `f(B^circ)/r`
proves

```math
E[(h_\beta-f(B)/r)_+]\to0.
```

Markov then gives every fixed lower-deviation probability in MF.5.

## Scope

The theorem correctly claims only a one-sided no-gain result.  A population
spike along `v` can increase pressure, so two-sided convergence need not
hold.  The proof is uniform for constant one-row mass and one distinguished
magnitude direction; it does not cover row correlations, vanishing row
mass, or a growing exceptional subspace.  These limitations are accurately
stated.

## Corrections

None required.
