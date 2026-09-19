# Wave 41: quadratic-completion structure beyond matched erasure

## 1. General external Ising form

Fix `S`, put `T=S^c`, `k=|T|`, fix an orientation `sigma` and retained
spin `y`, and write

```math
X_y(z)=\sigma z^{\mathsf T}A[T]z
       +2\sigma z^{\mathsf T}A[T,S]y,
\qquad z\in\{\pm1\}^T.
```

The actual normalized completion partition and component likelihood are

```math
G_{\beta,S}(y)=\frac{K_{\beta,S}(y)}{2^k}
=\mathbb E_z e^{\beta X_y(z)},
\qquad
F_{\beta,S}=\log G_{\beta,S},
\qquad
b_S=e^{-F_{\beta,S}}=G_{\beta,S}^{-1}.
```

Thus actual components are not arbitrary positive functions.  Extending the
visible boundary to `theta in R^S`,

```math
\Lambda_S(\theta)
=\log\mathbb E_z\exp\{\beta\sigma z^{\mathsf T}A[T]z
 +2\beta\sigma z^{\mathsf T}A[T,S]\theta\},
```

gives the exact continuous identities

```math
\nabla\Lambda_S
=2\beta\sigma A[S,T]\mathbb E_\theta z,
\qquad
\nabla^2\Lambda_S
=4\beta^2A[S,T]\operatorname{Cov}_\theta(z)A[T,S]\succeq0.
```

At zero temperature, uniform-spin Walsh orthogonality gives

```math
F_{0,S}=\partial_\beta F_{0,S}=0,
\qquad
\left.\partial_\beta^2F_{\beta,S}(y)\right|_{\beta=0}
=2k(k-1)+4\lVert A[T,S]y\rVert_2^2.
```

The nonconstant part of this second derivative is therefore a degree-two
Walsh polynomial whose coefficient matrix is the positive semidefinite Gram
matrix `A[S,T]A[T,S]`.  These identities hold at every deletion density.
They do not by themselves exclude the OR profile below.

## 2. A single genuine completion component realizes OR exactly

At one deletion, `T={i}`.  Put

```math
h_i(d)=\sum_{j\ne i}a_{ij}d_{ij}.
```

Directly summing the two outside spins gives

```math
G_i(d)=\cosh(2\beta h_i(d)),
\qquad
b_i(d)=\operatorname{sech}(2\beta h_i(d)).
```

Gauge the three edges from the omitted spin to three retained spins positive.
With retained gauge `y_0=1`,

```math
b_i(y_1,y_2)=\operatorname{sech}(2\beta(1+y_1+y_2)).
```

It equals `sech(6 beta)` at `(y_1,y_2)=(1,1)` and `sech(2 beta)`
at the other three points.  Hence it is an exact one-low/three-high OR
component, with score gap

```math
L_\beta=\log\frac{\cosh(6\beta)}{\cosh(2\beta)}\sim4\beta.
```

For `G=1/b`, its two-bit four-point determinant is exactly

```math
G_{++}G_{--}-G_{+-}G_{-+}=\sinh^2(4\beta)
```

up to the sign of the product of the two star coefficients.  Thus neither a
four-point identity, continuous log-partition convexity, nor the Walsh
constraint can exclude OR at the level of one selector.  The needed
structure is simultaneous-selector and base-Gibbs coupling.

## 3. Exact simultaneous one-deletion coupling

For the actual one-deletion mixture,

```math
U_\beta(d)=\frac1n\sum_i b_i(d)
=\frac1n\sum_i\operatorname{sech}(2\beta h_i(d)).
```

On the `i`-edge `d -> d^i`, the parent Gibbs base odds and the shared
omitted-`i` component obey

```math
\omega_i(d)
=\log\frac{\nu_\beta(d^i)}{\nu_\beta(d)}
=-4\beta h_i(d),
\qquad
\boxed{b_i(d)=\operatorname{sech}(\omega_i(d)/2).}
```

This absolute component normalization is absent from the generic matched
erasure axioms.  Since `b_i/n <= U(d),U(d^i) <= 1`, the score jump
`chi_i=log U(d^i)-log U(d)` satisfies

```math
\boxed{
|\chi_i|
\le\log n+\log\cosh(|\omega_i|/2).
}
```

If the conditional log odds cross zero during `0<=t<=1`, then
`|chi_i|>=|omega_i|`.  Using `log cosh(a/2)<=a/2` gives the exact
anti-resonance consequence

```math
\boxed{
\text{crossing edge}
\quad\Longrightarrow\quad
|\omega_i|\le2\log n,
\qquad |\chi_i|\le2\log n.
}
```

In particular, the Wave 40 uniform-base OR wall has `omega_i=0` but
`|chi_i|=L -> infinity` at fixed dimension, and cannot be an actual
one-deletion quadratic endpoint.  The obstruction is not OR shape; it is the
failure to couple the OR component simultaneously to the Gibbs edge odds and
the other deletion components.

There is also a pointwise global bound.  For a full spin `x`, let
`h_i=x_i(Ax)_i`.  Then, using `sech u >= exp(-|u|)`, Jensen, and the Boolean
bilinear norm `B(A)<=2Q(A)`,

```math
\begin{aligned}
U_\beta(d)
&\ge\exp\left\{-\frac{2\beta}{n}\sum_i|h_i(d)|\right\},\\
\sum_i|h_i(d)|
&=z^{\mathsf T}Ax\le B(A)\le2Q(A),
\end{aligned}
```

where `z_i=sign(h_i)x_i`.  Hence

```math
\boxed{
e^{-4\beta Q(A)/n}\le U_\beta(d)\le1,
\qquad
\operatorname{osc}(\log U_\beta)\le\frac{4\beta Q(A)}n.
}
```

For an exact minimizer and project temperature
`beta=Theta(n^(-1/2+c))`, this is `O(n^c)`.  It is far stronger than the
unconstrained growing-`L` wall at one deletion, although inserting it into
the existing anti-evanescence estimate still costs `exp(O(n^c))` and does
not close fixed-density migration.

The same edge inequality also yields the uniform one-deletion curvature
bound `v_{i,s}<=4(log n)^2`: if `|omega|<=2log n`, use
`|chi|<=2log n`; otherwise there is no crossing and
`|omega+s chi|>=|omega|-|chi|`, after which
`1/(4cosh^2(u/2))<=e^{-|u|}`.  This is exact but summing it over all
coordinates is much too coarse for the project target.

## 4. Exact order-four audit

Exhaustion of all `2^6=64` order-four signings gives `q_4=8` and 48 exact
minimizers.  Every row field has magnitude one or three.  On every state of
every exact minimizer, at most two of the four fields have magnitude three.
Writing

```math
p=\operatorname{sech}(2\beta),
\qquad q=\operatorname{sech}(6\beta),
```

therefore gives

```math
\frac{\max_dU_\beta(d)}{\min_dU_\beta(d)}
\le\frac{2p}{p+q}<2.
```

The displayed exact minimizer `A_4` attains the two-large-field pattern and
has `osc(log U)<log 2` for every beta.  Nevertheless each individual
omitted-vertex component can be gauged to the exact OR component of §2.

For comparison, the all-positive order-four signing has `Q=12`, not `q_4`.
It has a state with all four fields of magnitude three and a balanced state
with all four fields of magnitude one, so the full-mixture ratio is

```math
\frac{\max U}{\min U}
=\frac{\operatorname{sech}(2\beta)}
       {\operatorname{sech}(6\beta)}
=e^{4\beta+O(1)}.
```

Even here the actual Gibbs base suppresses the Wave 40 transient.  The low
score state has parent energies `+-12`, while the high balanced states have
energies `+-4`; their orientation-marginal base weights are proportional to
`cosh(12 beta)` and `cosh(4 beta)`.  The base disadvantage of the high-score
states is `e^{-8 beta}`, twice the endpoint score advantage `e^{4 beta}`.
Thus they remain rare throughout `0<=t<=1`.  Quadrature gives rapidly
vanishing `K^2`, rather than the `Theta(beta)` uniform-base wall.

At `beta=1`, the checker reports

```text
exact A4:             K^2 = 1.39984e-10
nonminimal all-plus:  K^2 = 1.98407e-9
```

and both decrease further at `beta=2,4`.  These values are numerical; the
component identities, enumeration, cap values, row-field counts, and ratio
bounds are exact.

## 5. Scope and surviving target

The quadratic completion identities do exclude the full Wave 40 wall at one
deletion, and they identify precisely why: a single selector can have an OR
completion profile, but its absolute component scale, all other deletion
components, and the parent Gibbs odds come from the same signing.

This does **not** prove the fixed-density target (10.1096) or directly bound
the signed migration (10.1083).  For `|T|>1`,

```math
G_{\beta,S}(y)
=\mathbb E_z e^{\beta X_y(z)}
```

contains an interacting outside Ising model.  The simple identity
`b_i=sech(omega_i/2)` is replaced by the exponential-tilt/composition
relations (10.745) and (10.916), with no current pointwise lower bound tied
to one parent edge odds.  The continuous Hessian and zero-temperature Walsh
identities allow an exact OR component and do not control movement of context
mass.  A fixed-density successor needs either:

1. a simultaneous-selector lower bound on `U_beta` or on crossing-edge
   component mass derived from external-completion composition and exact
   minimality; or
2. a signed covariance identity that retains the cancellation in (10.1083)
   instead of passing through `K`.

The present result is therefore an exact one-deletion exclusion and a precise
fixed-density obstruction, not a convergence proof.
